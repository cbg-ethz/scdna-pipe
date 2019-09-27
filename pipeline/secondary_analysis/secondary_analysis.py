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

    def remove_tenx_genomics_artifacts(self, bins):
        """
        Filters out the technical noise produced by 10x genomics sequencing artifacts.
        :param bins: The list of bins corresponding to technical noise
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

        cnvs = merge_chromosomes(h5f, key="cnvs")
        normalized_counts = merge_chromosomes(h5f)

        bin_size = h5f["constants"]["bin_size"][()]
        n_bins = normalized_counts.shape[1]
        bin_ids = [x for x in range(0, n_bins)]
        bin_df = pd.DataFrame(bin_ids, columns=["bin_ids"])

        bin_df["start"] = bin_df["bin_ids"] * bin_size
        bin_df["end"] = bin_df["start"] + bin_size
        print(bin_df.head())

        # exclude 10x artifact bins
        artifact_bins = np.loadtxt(bins, delimiter="\t").astype(bool)

        assert artifact_bins.shape[0] == normalized_counts.shape[1]

        print("artifact bins mask len")
        print(len(artifact_bins))
        print("artifact bins mask sum")
        print(sum(artifact_bins))

        print("normalized_counts matrix shape before & after filtering")
        print(normalized_counts.shape)
        normalized_counts = normalized_counts[:, ~artifact_bins]
        print(normalized_counts.shape)

        print("cnvs matrix shape before & after filtering")
        print(cnvs.shape)
        cnvs = cnvs[:, ~artifact_bins]
        print(cnvs.shape)

        print("bin_df shape before & after filtering")
        print(bin_df.shape)
        bin_df = bin_df[~artifact_bins]
        print(bin_df.shape)

        print("filtering chromosome stop positions")
        filtered_chr_stops = np.array(chr_stop_positions)[~artifact_bins]

        df_chr_stops = pd.DataFrame(columns=["chr"])
        for idx, val in enumerate(filtered_chr_stops):
            if val != None:
                # print((idx,val))
                df_chr_stops.loc[idx] = val

        cnvs = cnvs.astype("float")
        cnvs[cnvs < 0] = None

        print("writing output...")

        output_path = self.output_path + "/filtering/"

        np.savetxt(
            output_path + "/" + self.sample_name + "__filtered_counts_shape.txt",
            normalized_counts.shape
        )

        np.savetxt(
            output_path + "/" + self.sample_name + "__filtered_counts.csv",
            normalized_counts,
            delimiter=",",
        )

        np.savetxt(
            output_path + "/" + self.sample_name + "__filtered_cnvs.csv",
            cnvs,
            delimiter=",",
        )

        bin_df.to_csv(
            output_path + "/" + self.sample_name + "__bins_genome.tsv",
            sep="\t",
            index=False,
        )

        df_chr_stops.to_csv(
            output_path + "/" + self.sample_name + "__chr_stops.tsv", sep="\t"
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

    def plot_clusters(self, chr_stops_path, unique_cnvs_path):
        """
        Plots the copy number values
        for each cluster across the chromosome
        :param chr_stops_path: Path to the file containing the chromosome stop positions
        :param unique_cnvs_path: path to the file containing unique copy number profiles
        :return:
        """

        unique_cnvs = np.loadtxt(unique_cnvs_path, delimiter=',')
        # reshape if 1D array
        if (unique_cnvs.ndim == 1):
            unique_cnvs = np.reshape(unique_cnvs,(-1,unique_cnvs.shape[0]))
        cluster_ids = range(unique_cnvs.shape[0])

        chr_stops_df = pd.read_csv(chr_stops_path, sep='\t', index_col=0)
        chr_stops_df.columns = ["pos"]
        

        # the max val to be used in the plot
        # max_val = np.nanmax(cluster_means.values)
        max_val = 8  # max CN value

        output_path = os.path.join(self.output_path, "tree_learning")

        cmap = matplotlib.cm.get_cmap("Dark2")

        # use the formula below to get the distinct colors
        # color = cmap(float(i)/N)
        print("saving copy number figures by cluster.")
        for i, row in enumerate(unique_cnvs):
            plt.figure(figsize=(20, 6))
            ax = plt.scatter(
                y=unique_cnvs[i],
                x=range(len(row)),
                label="cluster id: " + str(i),
                color=cmap(
                    float(i + 1) / len(unique_cnvs)
                ),
                s=1,
            )
            plt.axis([None, None, -0.2, max_val])  # to make the axises same
            plt.legend(loc="upper left")
            plt.xticks([], [])
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
        for i, row in enumerate(unique_cnvs):
            ax = plt.scatter(
                y=row + (i + 1) / 10,
                x=range(len(row)),
                label="cluster id: " + str(i),
                color=cmap(
                    float(i + 1) / len(unique_cnvs)
                ),
                alpha=0.8,
                s=1,
            )
            plt.axis([None, None, -0.2, max_val])  # to make the axises same
            plt.legend(loc="upper left")
            plt.xticks([], [])
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

    def plot_heatmap(self):
        """
        Creates the heapmap of CN values per gene per cluster
        :return:
        """

        if self.communities_df is None:
            raise UnboundAttributeError(
                "The object attribute, namely communities_df is not set"
            )

        genes = pd.read_csv(self.genes_path, sep="\t")
        cnv_data = h5py.File(self.h5_path)

        gene_cn_df = self.get_gene_cluster_cn(genes, self.communities_df)
        gene_cn_df = gene_cn_df.T
        output_path = self.output_path + "/clustering/"

        gene_cn_df.to_csv(
            output_path + "/" + self.sample_name + "__cn_gene_cluster.tsv", sep="\t"
        )

        for i in range(1, 23):
            chr_df = gene_cn_df[gene_cn_df.index.str.startswith(f"{str(i)}/")]

            if chr_df.empty:
                continue

            chr_df.index = chr_df.index.to_series().str.split("/").str[1:].str.join("/")
            chr_df = chr_df.sort_index(axis=0)

            figure_width = chr_df.shape[0] / 2 + 1.5
            plt.figure(figsize=(8, figure_width))
            heatmap = sns.heatmap(
                chr_df,
                annot=True,
                cmap="bwr",
                vmin=0,
                vmax=4,
                xticklabels=True,
                yticklabels=True,
                cbar_kws={"ticks": [0, 1, 2, 3, 4]},
            )
            heatmap.set_title(f"Chromosome {str(i)}")
            heatmap.set_facecolor("#656565")
            # heatmap.get_legend().remove()
            heatmap = heatmap.get_figure()
            # heatmap.set_size_inches(19.7, 20.27)
            heatmap.savefig(
                f"{output_path}/{self.sample_name}__cn_genes_clusters_chr{str(i)}_heatmap.png"
            )
            plt.close()

        cnv_data.close()
