import h5py
import numpy as np
import pandas as pd
from .utils import merge_chromosomes
from .exceptions import UnboundAttributeError
import phenograph
from collections import Counter
from sklearn.preprocessing import normalize
import os
import matplotlib

if os.environ.get("DISPLAY", "") == "":
    # print("no display found. Using non-interactive Agg backend")
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
from tqdm import tqdm


class SecondaryAnalysis:
    """
    Class to perform the SecondaryAnalysis for a single-cell DNA sample.
    """

    def __init__(self, h5_path, genes_path, all_genes_path, sample_name, output_path):
        """

        :param h5_path: Path to the HDF5 file created by cellranger-dna (10x Genomics)
        :param genes_path: Path to the genes of interest
        :param all_genes_path: Path to all of the genes listed in the human genome reference
        :param sample_name: Name of the sample, to be added to the output names
        :param output_path: Desired path for the output files
        """

        self.sample_name = sample_name
        self.output_path = output_path
        self.h5_path = h5_path
        self.genes_path = genes_path
        self.all_genes_path = all_genes_path

        paths = [output_path, output_path + "/filtering/", output_path + "/clustering/"]
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)

    def extract_genomic_info(self):
        """
        Outputs the chromosome stops and bin start stop positions
        :return:
        """
        
        h5f = h5py.File(self.h5_path, "r")

        n_cells = h5f["cell_barcodes"].value.shape[0]
        all_chromosomes = list(h5f["normalized_counts"].keys())

        number_chromosomes = sorted([int(x) for x in sorted(all_chromosomes)[:-2]])
        ordered_chromosomes = [str(x) for x in number_chromosomes] + sorted(
            all_chromosomes
        )[-2:]

        chr_lengths = []
        for ch in ordered_chromosomes:
            chr_lengths.append(h5f["normalized_counts"][ch][0:n_cells, :].shape[1])

        chr_ends = np.cumsum(chr_lengths)

        chr_stop_positions = [None for x in range(0, chr_ends[-1])]
        for idx, pos in enumerate(chr_ends):
            chr_stop_positions[pos - 1] = ordered_chromosomes[
                idx
            ]  # -1 because it is a size info

        bin_size = h5f["constants"]["bin_size"][()]
        normalized_counts = merge_chromosomes(h5f)
        n_bins = normalized_counts.shape[1]
        bin_ids = [x for x in range(0, n_bins)]
        bin_df = pd.DataFrame(bin_ids, columns=["bin_ids"])

        bin_df["start"] = bin_df["bin_ids"] * bin_size
        bin_df["end"] = bin_df["start"] + bin_size
        print(bin_df.head())

        chr_stops_arr = np.array(chr_stop_positions)

        df_chr_stops = pd.DataFrame(columns=["chr"])
        for idx, val in enumerate(chr_stops_arr):
            if val != None:
                # print((idx,val))
                df_chr_stops.loc[idx] = val

        output_path = os.path.join(self.output_path, "genomic_coordinates")

        bin_df.to_csv(
            os.path.join(output_path, self.sample_name) + "__bins_genome.tsv",
            sep="\t",
            index=False,
        )

        df_chr_stops.to_csv(
            os.path.join(output_path, self.sample_name) + "__chr_stops.tsv", sep="\t"
        )

        print("Output written to: " + output_path)

        h5f.close()
        

    def remove_tenx_genomics_artifacts(self, bins):
        """
        Filters out the technical noise produced by 10x genomics sequencing artifacts.
        :param bins: The list of bins corresponding to technical noise
        :return:
        """
        h5f = h5py.File(self.h5_path, "r")

        n_cells = h5f["cell_barcodes"].value.shape[0]
        all_chromosomes = list(h5f["normalized_counts"].keys())

        normalized_counts = merge_chromosomes(h5f)

        bin_size = h5f["constants"]["bin_size"][()]
        n_bins = normalized_counts.shape[1]
        bin_ids = [x for x in range(0, n_bins)]
        bin_df = pd.DataFrame(bin_ids, columns=["bin_ids"])

        bin_df["start"] = bin_df["bin_ids"] * bin_size
        bin_df["end"] = bin_df["start"] + bin_size
        print(bin_df.head())

        # exclude unmappable bins
        is_mappable = []
        for chr in all_chromosomes:
            is_mappable = np.concatenate([is_mappable,h5f['genome_tracks']['is_mappable'][chr][:]])

        unmappable = ~np.array(is_mappable, dtype=bool)

        print(f"unmappable bins len: {len(unmappable)}")
        print(f"unmappable bins sum: {sum(unmappable)}")

        # exclude 10x artifact bins
        artifact_bins = np.loadtxt(bins, delimiter="\t").astype(bool)

        print(f"artifact bins mask len: {len(artifact_bins)}")
        print(f"artifact bins mask sum: {sum(artifact_bins)}")

        assert artifact_bins.shape[0] == normalized_counts.shape[1]
        assert unmappable.shape[0] == artifact_bins.shape[0]
        to_filter_out = np.logical_or(unmappable, artifact_bins)
        print(f"combined filter len: {len(to_filter_out)}")
        print(f"combined filter sum: {sum(to_filter_out)}")

        print("normalized_counts matrix shape before & after filtering")
        print(normalized_counts.shape)
        normalized_counts = normalized_counts[:, ~to_filter_out]
        print(normalized_counts.shape)

        print("writing output...")

        output_path = os.path.join(self.output_path, "filtering")

        np.savetxt(
            os.path.join(output_path, self.sample_name) + "__excluded_bins.csv",
            to_filter_out,
            delimiter=','
        )

        np.savetxt(
            output_path + "/" + self.sample_name + "__filtered_counts_shape.txt",
            normalized_counts.shape
        )

        np.savetxt(
            output_path + "/" + self.sample_name + "__filtered_counts.csv",
            normalized_counts,
            delimiter=",",
        )

        print("Output written to: " + output_path)

        h5f.close()

    def apply_phenograph(self, normalised_regions_path, n_jobs=1):
        """
        Runs the phenograph clustering algorithm on the object and alters its fields
        :param n_jobs: The number of threads for clustering
        :param normalised_regions_path: The path to the normalised regions
        :return:
        """

        print("loading the normalised regions...")
        normalised_regions = np.loadtxt(normalised_regions_path, delimiter=',')
        print(f"shape of normalised regions: {normalised_regions.shape}")

        n_cells = normalised_regions.shape[0]

        print(f"n_cells: {str(n_cells)}")
        n_neighbours = int(n_cells / 10)
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        communities, graph, Q = phenograph.cluster(
            data=normalised_regions, k=n_neighbours, n_jobs=n_jobs, jaccard=True
        )

        # computing the distance matrix from graph
        arr = graph.toarray()
        arr_full = arr + arr.T
        np.fill_diagonal(arr_full, 1)
        dist = (arr_full - arr_full.max()) * (-1)
        np.fill_diagonal(dist, 0)

        print(f"shape of the distance matrix: {dist.shape}")

        dist_fname = (
            os.path.join(self.output_path, "clustering", self.sample_name) + "__phenograph_distance.csv"
        )
        np.savetxt(fname=dist_fname, X=dist, delimiter=",")

        print(f"Communities: {communities}")
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[
            ["cell_barcode", "cluster"]
        ]

        output_path = os.path.join(self.output_path, "clustering", self.sample_name)
        print(f"output path: {output_path}")
        communities_df.to_csv(
            output_path + "__clusters_phenograph_assignment.tsv",
            sep="\t",
            index=False
        )

        # write the modularity score, for stability
        with open(output_path + "__clustering_score.txt", "w") as f:
            f.write(str(Q))

    def add_filtered_bins_back(self, unique_cnvs_path, bin_mask_path):
        """
        Adds the filtered bins back to the inferred cnvs
        :param unique_cnvs_path: path to the file containing unique copy number profiles
        :param bin_mask_path: path to the excluded bins file
        :return:
        """
        unique_cnvs = np.loadtxt(unique_cnvs_path, delimiter=',')
        print(f"unique_cnvs shape: {unique_cnvs.shape}")
        bin_mask = np.loadtxt(bin_mask_path, delimiter=',')
        print(f"bin_mask shape: {bin_mask.shape}")

        cnvs_mat = []
        cnvs_counter = 0

        for bin_idx, bin_val in enumerate(bin_mask):
            c_row = []
            if bin_val:
                for c_id in range(unique_cnvs.shape[0]):        
                    c_row.append(None)
            else:
                for c_id in range(unique_cnvs.shape[0]):
                    c_row.append(unique_cnvs[c_id][cnvs_counter])
                cnvs_counter += 1
            cnvs_mat.append(c_row)

        cnvs_arr =  np.array(cnvs_mat, dtype=float).T
        print(f"cnvs_arr shape: {cnvs_arr.shape}")

        print("writing the inferred cnvs...")
        np.savetxt(
            os.path.join(self.output_path, "inferred_cnvs", self.sample_name) + "__inferred_cnvs.csv",
            cnvs_arr,
            delimiter=",",
            fmt="%s"
        )

    def plot_clusters(self, chr_stops_path, cnvs_arr_path):
        """
        Plots the copy number values
        for each cluster across the chromosome
        :param chr_stops_path: Path to the file containing the chromosome stop positions
        :param cnvs_array_path: path to the file containing unique copy number profiles (including missing values)
        :return:
        """

        cnvs_arr = np.loadtxt(cnvs_arr_path, delimiter=',')

        chr_stops_df = pd.read_csv(chr_stops_path, sep='\t', index_col=0)
        chr_stops_df.columns = ["pos"]
        

        # the max val to be used in the plot
        # max_val = np.nanmax(cluster_means.values)
        max_val = 8  # max CN value

        output_path = os.path.join(self.output_path, "inferred_cnvs")

        cmap = matplotlib.cm.get_cmap("Dark2")

        # use the formula below to get the distinct colors
        # color = cmap(float(i)/N)
        print("saving copy number figures by cluster.")
        for i, row in enumerate(cnvs_arr):
            plt.figure(figsize=(20, 6))
            ax = plt.scatter(
                y=cnvs_arr[i],
                x=range(len(row)),
                label="cluster id: " + str(i),
                color=cmap(
                    float(i + 1) / len(cnvs_arr)
                ),
                s=1,
            )
            plt.axis([None, None, -0.2, max_val])  # to make the axises same
            plt.legend(loc="upper left")
            plt.xticks([], [])
            plt.text(0, -0.5, "0")
            for index, row in chr_stops_df.iterrows():
                plt.text(index, -0.5, "chr " + row["pos"], rotation=90)
            plt.savefig(
                os.path.join(output_path, self.sample_name) + "__cluster_profile_"
                + str(i)
                + ".png"
            )
            plt.close()

        print("Saving overlapping copy number profile figures by cluster.")
        plt.figure(figsize=(20, 6))
        for i, row in enumerate(cnvs_arr):
            ax = plt.scatter(
                y=cnvs_arr[i],
                x=range(len(row)),
                label="cluster id: " + str(i),
                color=cmap(
                    float(i + 1) / len(cnvs_arr)
                ),
                alpha=0.8,
                s=1,
            )
            plt.axis([None, None, -0.2, max_val])  # to make the axises same
            plt.legend(loc="upper left")
            plt.xticks([], [])
            plt.text(0, -0.5, "0")
            for index, row in chr_stops_df.iterrows():
                plt.text(index, -0.5, "chr " + row["pos"], rotation=90)
        plt.savefig(
            os.path.join(output_path, self.sample_name)+"__cluster_profile_overlapping.png"
        )
        plt.close()

    def get_gene_cluster_cn(self, genes, cell_assignment, with_chr_names=True):
        """
        Creates and returns the dataframe of copy numbers, genes by cluster ids
        :param genes: The input list of genes to be specified
        :param cell_assignment: Dataframe containing cell_barcode and cell cluster id
        :param with_chr_names: Boolean variable to denote whether to print the chromosome names or not
        :return: CN dataframe of genes by clusters
        """

        cnv_data = h5py.File(self.h5_path)

        n_cells = cell_assignment.shape[0]
        print("n_cells: ", n_cells)

        bin_size = cnv_data["constants"]["bin_size"][()]
        print("bin_size: ", bin_size)

        cluster_ids = sorted(list(set(cell_assignment.cluster)))
        print("cluster_ids: ", cluster_ids)

        gene_cn_df = pd.DataFrame(index=cluster_ids)
        # for each gene
        for index, row in tqdm(genes.iterrows()):
            start_bin = int(row["Gene start (bp)"] / bin_size)
            stop_bin = int(row["Gene end (bp)"] / bin_size)
            chromosome = str(row["Chromosome/scaffold name"])
            gene_name = row["Gene name"]
            # print(start_bin, stop_bin)
            median_copy_numbers = []
            for c_id in cluster_ids:
                # get all the cells that belong to that cluster
                cells = cell_assignment[cell_assignment.cluster == c_id][
                    "cell_barcode"
                ].values.tolist()
                cn_states = cnv_data["cnvs"][chromosome][:n_cells][
                    cells, start_bin : stop_bin + 1
                ]

                cn_states = cn_states.astype("float")
                cn_states[cn_states < 0] = np.nan

                min_cn_cell_bin = np.nanmin(
                    cn_states, axis=1
                )  # all the bins within the gene, min due to biology
                median_cn_cell = np.nanmedian(
                    min_cn_cell_bin
                )  # median value across all cells
                median_copy_numbers.append(median_cn_cell)

            # print(mean_copy_numbers)
            if with_chr_names:
                gene_cn_df[chromosome + "/" + gene_name] = median_copy_numbers
            else:
                gene_cn_df[gene_name] = median_copy_numbers

        return gene_cn_df

    def create_cn_cluster_h5(self):
        """
        Creates the HDF5 for copy number values per cluster
        :return:
        """

        if self.communities_df is None:
            raise UnboundAttributeError(
                "The object attribute, namely communities_df is not set"
            )

        all_genes = pd.read_csv(self.all_genes_path, sep="\t")
        cnv_data = h5py.File(self.h5_path)

        all_gene_cn_df = self.get_gene_cluster_cn(
            all_genes, self.communities_df, with_chr_names=False
        )
        output_path = self.output_path + "/clustering/"

        cn_cluster_h5 = h5py.File(
            output_path + "/" + self.sample_name + "__cn_cluster.h5", "w"
        )
        gene_attributes = cn_cluster_h5.create_group("gene_attrs")

        column_values = np.array(all_gene_cn_df.columns.values, dtype="S16")

        gene_attributes.create_dataset("gene_names", data=column_values)
        cn_cluster_h5.create_dataset("matrix", data=all_gene_cn_df.values)

        cnv_data.close()
        cn_cluster_h5.close()
