import h5py
import numpy as np
import pandas as pd
from utils import merge_chromosomes
import phenograph
from collections import Counter
from sklearn.preprocessing import normalize
import os
import matplotlib
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns


class SecondaryAnalysis:

    def __init__(self, h5_path, genes_path, sample_name, output_path):
        self.filtered_normalized_counts = None
        self.filtered_cnvs = None
        self.chr_stops = None
        self.bin_positions = None
        self.clustering_distance = None
        self.sample_name = sample_name
        self.output_path = output_path
        self.h5_path = h5_path
        self.genes_path = genes_path

        paths = [output_path, output_path+'/filtering/', output_path + '/clustering/']
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)

    def remove_tenx_genomics_artifacts(self, bins):

        h5f = h5py.File(self.h5_path, 'r')

        n_cells = h5f['cell_barcodes'].value.shape[0]
        all_chromosomes = list(h5f['normalized_counts'].keys())

        number_chromosomes = sorted([int(x) for x in sorted(all_chromosomes)[:-2]])
        ordered_chromosomes = [str(x) for x in number_chromosomes] + sorted(all_chromosomes)[-2:]

        chr_lengths = []
        for ch in ordered_chromosomes:
            chr_lengths.append(h5f['normalized_counts'][ch][0:n_cells, :].shape[1])

        chr_ends = np.cumsum(chr_lengths)

        chr_stop_positions = [None for x in range(0, chr_ends[-1])]
        for idx, pos in enumerate(chr_ends):
            chr_stop_positions[pos - 1] = ordered_chromosomes[idx]  # -1 because it is a size info

        cnvs = merge_chromosomes(h5f, key='cnvs')
        normalized_counts = merge_chromosomes(h5f)

        bin_size = h5f["constants"]["bin_size"][()]
        n_bins = normalized_counts.shape[1]
        bin_ids = [x for x in range(0, n_bins)]
        bin_df = pd.DataFrame(bin_ids, columns=["bin_ids"])

        bin_df["start"] = bin_df["bin_ids"] * bin_size
        bin_df["end"] = bin_df["start"] + bin_size
        print(bin_df.head())

        # exclude 10x artifact bins
        artifact_bins = np.loadtxt(bins, delimiter='\t').astype(bool)

        assert (artifact_bins.shape[0] == normalized_counts.shape[1])

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
            if (val != None):
                # print((idx,val))
                df_chr_stops.loc[idx] = val

        '''
        cnvs
        Denote the NaN values with None
        Make the imputed values non-negative
        '''

        # cnvs==-127 means imputed to 0
        cnvs[cnvs == -127] = 0
        cnvs[cnvs == -128] = 129
        cnvs = np.abs(cnvs)
        cnvs = cnvs.astype('float')
        cnvs[cnvs == 129] = None

        self.filtered_cnvs = cnvs
        self.chr_stops = df_chr_stops
        self.filtered_normalized_counts = normalized_counts
        self.bin_positions = bin_df

        print("writing output...")
        output_path = self.output_path + '/filtering/'
        np.savetxt(output_path + '/' + self.sample_name + "__filtered_counts.tsv", normalized_counts,
                   delimiter='\t')

        np.savetxt(output_path + '/' + self.sample_name + "__filtered_cnvs.tsv", cnvs, delimiter='\t')

        bin_df.to_csv(output_path + '/' + self.sample_name + "__bins_genome.tsv", sep='\t', index=False)

        df_chr_stops.to_csv(output_path + '/' + self.sample_name + "__chr_stops.tsv", sep='\t')

        print("Output written to: " + output_path)

        h5f.close()

    def apply_phenograph(self, n_jobs=16):

        # points in the knn neighbourhood are weighted by the distance
        weight = 'distance'

        filtered_counts = self.filtered_normalized_counts
        normalized_filtered_counts = normalize(filtered_counts, axis=1, norm='l1')

        cnvs = self.filtered_cnvs

        n_cells = normalized_filtered_counts.shape[0]

        print("n_cells: " + str(n_cells))
        n_neighbours = int(n_cells / 10)
        print("n_neighbours: " + str(n_neighbours))
        communities, graph, Q = phenograph.cluster(data=normalized_filtered_counts, k=n_neighbours, n_jobs=n_jobs,
                                                   jaccard=True)

        # computing the distance matrix from graph
        arr = graph.toarray()
        arr_full = arr + arr.T
        np.fill_diagonal(arr_full, 1)
        dist = (arr_full - arr_full.max()) * (-1)
        np.fill_diagonal(dist, 0)

        print("shape of the distance matrix:")
        print(dist.shape)

        # write dist to file
        # later use dist for all of the plots
        self.clustering_distance = dist
        # dist_fname = args.output_path + '/' + args.sample_name + "_phenograph_distance.csv"
        # np.savetxt(fname=dist_fname, X=dist, delimiter=',')

        print(communities)  # one of the outputs

        communities_df = pd.DataFrame(communities, columns=['cluster'])
        communities_df['cell_barcode'] = communities_df.index
        communities_df = communities_df[['cell_barcode', 'cluster']]  # order the columns
        communities_df.head()

        output_path = self.output_path + '/clustering/'
        communities_df.to_csv(output_path + '/' + self.sample_name + "__clusters_phenograph_assignment.tsv",
                              sep='\t', index=False)

        # write the modularity score, for stability
        f = open(output_path + '/' + self.sample_name + "__clustering_score.txt", 'w')
        f.write(str(Q))
        f.close()

        cells_by_cluster = []

        community_dict = dict((Counter(communities)))

        community_ids = sorted(list(community_dict))

        with open(output_path + '/' + self.sample_name + "__cluster_sizes.txt", 'w') as community_dict_file:
            community_dict_file.write(str(community_dict))

        for cluster in community_ids:
            cells_by_cluster.append(filtered_counts[communities == cluster])

        avg_clusters = [m.mean(0) for m in cells_by_cluster]

        avg_clusters_df = pd.DataFrame(avg_clusters)

        avg_clusters_df['cluster_ids'] = community_ids  # add the community_ids

        avg_clusters_df.to_csv(output_path + '/' + self.sample_name + "__clusters_phenograph_count_profiles.tsv",
                               sep='\t', index=False, header=True)

        cnvs_per_cluster = []
        for cluster in community_ids:
            cnvs_per_cluster.append(cnvs[communities == cluster])
        cn_avg_clusters = [np.nanmean(c, axis=0) for c in cnvs_per_cluster]

        cn_avg_clusters_df = pd.DataFrame(cn_avg_clusters)
        cn_avg_clusters_df['cluster_ids'] = community_ids

        cn_avg_clusters_df.to_csv(output_path + '/' + self.sample_name + "__clusters_phenograph_cn_profiles.tsv",
                                  sep='\t', index=False, header=True)

        self.plot_clusters(cluster_means=cn_avg_clusters_df, dist=dist, communities=communities)
        self.plot_heatmap(communities_df)
    def plot_clusters(self, cluster_means, dist, communities):

        # the max val to be used in the plot
        # max_val = np.nanmax(cluster_means.values)
        max_val = 12  # max CN value

        output_path = self.output_path + '/clustering/'

        cmap = matplotlib.cm.get_cmap('Dark2')
        tsne = TSNE(n_components=2, perplexity=30, metric='precomputed').fit_transform(dist)
        df_tsne = pd.DataFrame(tsne)
        df_tsne['cluster'] = communities
        df_tsne['color'] = (df_tsne['cluster'] + 1) / len(cluster_means.index)  # +1 because outliers are -1
        ax = df_tsne.plot(kind='scatter', x=0, y=1, c=cmap(df_tsne['color']), figsize=(10, 8), colorbar=False,
                          grid=True, title='Phenograph Clusters on CNV Data')
        fig = ax.get_figure()

        fig.savefig(output_path + '/' + self.sample_name + "__tsne_output.png")

        chr_stops_df = self.chr_stops
        chr_stops_df.columns = ["pos"]

        # use the formula below to get the distinct colors
        # color = cmap(float(i)/N)
        for i, cluster_idx in enumerate(cluster_means.index):
            plt.figure(figsize=(20, 6))
            ax = plt.plot(cluster_means.iloc[i].values, label="cluster id: " + str(cluster_idx),
                          color=cmap(float(cluster_idx + 1) / len(cluster_means.index)))
            plt.axis([None, None, 0, max_val])  # to make the axises same
            plt.legend(loc='upper left')
            plt.xticks([], [])
            for index, row in chr_stops_df.iterrows():
                plt.text(index, -0.5, "chr " + row['pos'], rotation=90)
            plt.savefig(output_path + '/' + self.sample_name + "__cluster_profile_" + str(cluster_idx) + ".png")

        plt.figure(figsize=(20, 6))
        for i, cluster_idx in enumerate(cluster_means.index):
            ax = plt.plot(cluster_means.iloc[i].values, label="cluster id: " + str(cluster_idx),
                          color=cmap(float(cluster_idx + 1) / len(cluster_means.index)), alpha=0.6)
            plt.axis([None, None, 0, max_val])  # to make the axises same
            plt.legend(loc='upper left')
            plt.xticks([], [])
            for index, row in chr_stops_df.iterrows():
                plt.text(index, -0.5, "chr " + row['pos'], rotation=90)
        plt.savefig(output_path + '/' + self.sample_name + "__cluster_profile_overlapping.png")

    def plot_heatmap(self, cell_assignment):

        genes = pd.read_csv(self.genes_path, sep='\t')
        cnv_data = h5py.File(self.h5_path)

        n_cells = cell_assignment.shape[0]
        print("n_cells: ", n_cells)

        bin_size = cnv_data['constants']['bin_size'][()]
        print("bin_size: ", bin_size)

        cluster_ids = sorted(list(set(cell_assignment.cluster)))
        print("cluster_ids: ", cluster_ids)

        gene_cn_df = pd.DataFrame(index=cluster_ids)
        # for each gene
        for index, row in genes.iterrows():
            start_bin = int(row['Gene start (bp)'] / bin_size)
            stop_bin = int(row['Gene end (bp)'] / bin_size)
            chromosome = str(row['Chromosome/scaffold name'])
            gene_name = row['Gene name']
            # print(start_bin, stop_bin)
            mean_copy_numbers = []
            for c_id in cluster_ids:
                # get all the cells that belong to that cluster
                cells = cell_assignment[cell_assignment.cluster == c_id]['cell_barcode'].values.tolist()
                cn_states = cnv_data['cnvs'][chromosome][:n_cells][cells, start_bin:stop_bin + 1]

                # -127 means imputed to 0
                cn_states[cn_states == -127] = 0
                cn_states[cn_states == -128] = 129
                cn_states = np.abs(cn_states)
                cn_states = cn_states.astype('float')
                cn_states[cn_states == 129] = np.nan

                min_cn_cell_bin = np.nanmin(cn_states, axis=1)
                avg_cn_cell = np.nanmean(min_cn_cell_bin)
                mean_copy_numbers.append(avg_cn_cell)

            # print(mean_copy_numbers)
            gene_cn_df[chromosome + '/' + gene_name] = mean_copy_numbers

        output_path = self.output_path + '/clustering/'

        gene_cn_df.T.to_csv(output_path + '/' + self.sample_name + '__cn_gene_cluster.tsv', sep='\t')
        heatmap = sns.heatmap(gene_cn_df.T, annot=True, cmap='bwr', vmin=0, vmax=4, xticklabels=True, yticklabels=True,
                              cbar_kws={"ticks": [0, 1, 2, 3, 4]})
        heatmap.get_legend().remove()
        heatmap = heatmap.get_figure()
        heatmap.set_size_inches(19.7, 20.27)
        heatmap.savefig(output_path + '/' + self.sample_name + "__cn_genes_clusters_heatmap.png")

        cnv_data.close()

            
