from scgenpy.plotting.plotting import *
from scgenpy.process_cnvs.process_cnvs import *
from scgenpy.preprocessing.utils import sort_chromosomes

rule plots:
    input:
        sorted_inferred_cnvs_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_cnvs_bins.png",
        lib_sizes = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clone_lib_sizes.png",
        cluster_tree_genes_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.png",
        overlapping_cluster_plot = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_profile_overlapping.png",
        heatmap_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__heatmap_cnvs.png"
    run:
        print("Getting plots")

rule plot_cnv_matrix:
    input:
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins_final.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions_final.csv",
        candidate_segmented_regions = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions_candidate.txt",
        final_segmented_regions = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions_final.txt",
        excluded_bins = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt",
        inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_final.csv",
        is_outlier = os.path.join(analysis_path, "filtering", analysis_prefix) + "_is_outlier.txt",
        all_cells_normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__all_cells_normalised_bins_final.csv",
        cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_final.csv",
    output:
        sorted_normalised_counts_heatmap_candidate_bp = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_normalised_counts_bins_candidate_bps.png",
        sorted_normalised_counts_heatmap_final_bp = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_normalised_counts_bins_final_bps.png",
        sorted_normalised_counts_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_normalised_counts_bins.png",
        sorted_inferred_cnvs_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_cnvs_bins.png"
    benchmark:
        "benchmarks/plot_cnv_matrix.tsv"
    run:
        from scipy.cluster.hierarchy import ward, leaves_list
        from scipy.spatial.distance import pdist

        is_outlier = np.loadtxt(input.is_outlier, delimiter=',').astype(bool).ravel()
        print("loading the all_cells normalised bins...")
        all_cells_normalised_bins = np.loadtxt(input.all_cells_normalised_bins, delimiter=',')
        print(f"shape of all_cells_normalised bins: {all_cells_normalised_bins.shape}")

        print("loading the normalised bins...")
        normalised_bins = np.loadtxt(input.normalised_bins, delimiter=',')
        print(f"shape of normalised bins: {normalised_bins.shape}")

        print("loading the inferred cnvs...")
        inferred_cnvs = np.loadtxt(input.inferred_cnvs, delimiter=',')
        print(f"shape of inferred cnvs: {inferred_cnvs.shape}")

        print("loading the cell labels...")
        cell_labels = np.loadtxt(input.cell_labels, delimiter=',')

        # Annotate with chromosome coordinates
        chr_indicator = pd.read_csv(input.bin_chr_indicator, header=None).values.ravel()

        # Remove excluded bins
        excluded_bins = pd.read_csv(input.excluded_bins, header=None).values.ravel().astype(bool)
        chr_indicator_filtered = np.array(chr_indicator)[~excluded_bins]
        chromosome_list = utils.sort_chromosomes(np.unique(chr_indicator))
        chr_stops_filtered_bins = [np.where(chr_indicator_filtered==chr)[0][-1] for chr in chromosome_list]
        chr_stops_filtered_bins_dict = dict(zip(chromosome_list, chr_stops_filtered_bins))

        # Now plot data with the inferred CNV information (i.e. cells ordered and annotated by inferred clone)
        clustered_cells = np.argsort(cell_labels)
        clustered_labels = cell_labels[clustered_cells]
        clustered_inferred_cnvs = inferred_cnvs[clustered_cells]
        clustered_normalised_bins = normalised_bins[clustered_cells]

        # Add the outlier cells
        outlier_cnvs = np.nan * np.ones((np.count_nonzero(is_outlier), clustered_inferred_cnvs.shape[1]))
        clustered_inferred_cnvs_full = np.vstack([outlier_cnvs, clustered_inferred_cnvs])

        outlier_normbins = all_cells_normalised_bins[is_outlier]
        clustered_normalised_bins_full = np.vstack([outlier_normbins, clustered_normalised_bins])

        outlier_labels = np.array(['-'] * np.count_nonzero(is_outlier))
        clustered_labels = clustered_labels.astype(int).astype(str)
        clustered_labels_full = np.hstack([outlier_labels, clustered_labels]).tolist()

        candidate_bps = np.loadtxt(input.candidate_segmented_regions)
        print(f"number of candidate breakpoints: {len(candidate_bps)}")

        final_bps = np.loadtxt(input.final_segmented_regions)
        print(f"number of final breakpoints: {len(final_bps)}")

        print("plotting the heatmap of normalised counts ordered by clone...")
        scicone.plotting.plot_matrix(clustered_normalised_bins_full, mode='data', cbar_title='  Normalized  \ncounts', vmax=2,
                            chr_stops_dict=chr_stops_filtered_bins_dict,
                            labels=clustered_labels_full,
                            figsize=(24,8), output_path=output.sorted_normalised_counts_heatmap, dpi=100)

        print("plotting the heatmap of normalised counts ordered by clone with candidate breakpoints...")
        scicone.plotting.plot_matrix(clustered_normalised_bins_full, mode='data', cbar_title='  Normalized  \ncounts', vmax=2,
                            chr_stops_dict=chr_stops_filtered_bins_dict,
                            labels=clustered_labels_full,
                            bps=candidate_bps,
                            figsize=(24,8), output_path=output.sorted_normalised_counts_heatmap_candidate_bp, dpi=100)

        print("plotting the heatmap of normalised counts ordered by clone with final breakpoints...")
        scicone.plotting.plot_matrix(clustered_normalised_bins_full, mode='data', cbar_title='  Normalized  \ncounts', vmax=2,
                            chr_stops_dict=chr_stops_filtered_bins_dict,
                            labels=clustered_labels_full,
                            bps=final_bps,
                            figsize=(24,8), output_path=output.sorted_normalised_counts_heatmap_final_bp, dpi=100)

        print("plotting the heatmap of inferred cnvs ordered by clone...")
        scicone.plotting.plot_matrix(clustered_inferred_cnvs_full, mode='cnv', cbar_title='Copy number\nvalues', vmax=4,
                            chr_stops_dict=chr_stops_filtered_bins_dict,
                            labels=clustered_labels_full,
                            figsize=(24,8), output_path=output.sorted_inferred_cnvs_heatmap, dpi=100)

rule plot_clone_lib_sizes:
    input:
        filtered_counts = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_final.csv",
    output:
        lib_sizes = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clone_lib_sizes.png",
    run:
        filtered_counts = np.loadtxt(input.filtered_counts, delimiter=',')
        cell_node_labels = np.loadtxt(input.cell_labels, delimiter=',')

        lib_sizes = np.sum(filtered_counts, axis=1).reshape(-1,1)
        # Sort cells by clone for visualization
        lib_sizes, labels = scicone.utils.cluster_clones(lib_sizes, cell_node_labels, within_clone=False)

        fig = plt.figure(figsize=(24,8), dpi=100)
        for i, label in enumerate(np.unique(labels)):
            xx = np.arange(np.where(label==labels)[0][0], np.where(label==labels)[0][-1]+1)
            plt.scatter(xx, lib_sizes[label==labels], label=str(int(label)), color=scicone.constants.LABEL_CPAL_HEX[int(label)],
                     s=80, alpha=0.6)
            plt.hlines(np.mean(lib_sizes[label==labels]), xx[0], xx[-1]+1,
                        color=scicone.constants.LABEL_CPAL_HEX[int(label)], linewidth=5, linestyle='-')
        leg = plt.legend(fontsize=24)
        leg.set_title('Subclone',prop={'size':24})
        plt.ylabel('Total read counts',fontsize=24)
        plt.xlabel('Cells',fontsize=24)
        ax = plt.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(22)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(22)
        plt.savefig(output.lib_sizes, bbox_inches="tight")


rule visualise_trees:
    params:
        max_genes_per_line = config['plotting']['trees']['max_genes_per_line'],
        gender = config['gender']
    input:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final.txt",
        cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_final.tsv",
        cluster_tree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final.json",
        bin_gene_region_df_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv",
        genes_to_highlight_path = disease_genes_path
    output:
        cluster_tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.graphviz",
        cluster_tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.png",
        cluster_tree_genes_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.graphviz",
        cluster_tree_genes_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.png"
    run:
        # Node labels
        with open(input.cluster_tree_json) as json_file:
            cluster_tree = json.load(json_file)
        node_labels = dict()
        for node in cluster_tree:
            if cluster_tree[node]['label'] != "":
                node_labels[node] = cluster_tree[node]['label']

        # Node sizes
        cell_node_assignments = np.loadtxt(input.cell_node_assignments, delimiter='\t')

        # Remove outliers
        cell_node_assignments = cell_node_assignments[cell_node_assignments != '-']
        keys, values = [arr.tolist() for arr in np.unique(cell_node_assignments, return_counts=True)]
        node_sizes = dict(zip(np.array(keys).astype(int).astype(str).tolist(), values))

        # Tree with region labels
        regions_tree = tree_to_graphviz(input.cluster_tree, node_sizes=node_sizes, node_labels=node_labels, gender=params.gender)
        with open(output.cluster_tree_graphviz, "w") as file:
            for line in regions_tree:
                file.write(f"{line}\n")
        plot_tree_graphviz(output.cluster_tree_graphviz, output.cluster_tree_figure)

        # Tree with gene labels
        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df_path, index_col=0)
        genes_to_highlight = pd.read_csv(input.genes_to_highlight_path, sep="\t", header=None)
        genes_to_highlight.columns = ["gene"]
        genes_to_highlight = genes_to_highlight.gene.values.tolist()
        genes_tree = tree_to_graphviz(input.cluster_tree, node_sizes=node_sizes, gene_labels=True, bin_gene_region_df=bin_gene_region_df,
                                genes_to_highlight=genes_to_highlight, max_genes_per_line=params.max_genes_per_line, node_labels=node_labels,
                                gender=params.gender)
        with open(output.cluster_tree_genes_graphviz, "w") as file:
            for line in genes_tree:
                file.write(f"{line}\n")
        plot_tree_graphviz(output.cluster_tree_genes_graphviz, output.cluster_tree_genes_figure)

rule plot_cluster_cnvs:
    params:
        offset_sizes = config['plotting']['profiles']['offset_sizes'],
        s = config['plotting']['profiles']['s'],
        max_amp_val = config['inference']['copy_number_limit'],
        ploidy = config['inference']['ploidy']
    input:
        inferred_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__unique_cnvs_final.csv",
        chr_stops = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv"
    output:
        overlapping_cluster_plot = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_profile_overlapping.png"
    benchmark:
        "benchmarks/plot_cluster_cnvs.tsv"
    run:
        cnvs_arr = np.loadtxt(input.inferred_cnvs, delimiter=',')
        chr_stops = pd.read_csv(input.chr_stops, sep="\t")

        ymax = params.ploidy + params.max_amp_val + 1

        output_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix)
        plot_cluster_cnvs(cnvs_arr, chr_stops, output_path=output_path, ymax=ymax, offset_sizes=params.offset_sizes, s=params.s)


try:
    max_genes = config['plotting']['heatmaps']['max_genes']
except KeyError:
    max_genes = 20

rule plot_heatmaps:
    input:
        general_gene_lists_path = os.path.join(gene_lists_path, 'general'),
        disease_genes_path = disease_genes_path,
        bin_gene_region_df_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv"
    output:
        gene_cn_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_gene_df.csv",
        cn_gene_df_roche_gene_list = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_gene_df_roche_gene_list.csv",
        heatmap_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__heatmap_cnvs.png"
    benchmark:
        "benchmarks/plot_heatmap.tsv"
    run:
        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df_path, index_col=0, low_memory=False)

        # General lists
        for dirpath,_,filenames in os.walk(input.general_gene_lists_path):
           for f in filenames:
               print(f'Creating heatmap for {f}')
               gene_list_path = os.path.abspath(os.path.join(dirpath, f))
               genes_list = pd.read_csv(gene_list_path, sep="\t", header=None)
               genes_list.columns = ["gene"]
               genes_list = genes_list.gene.values.tolist()

               gene_cn_df_imputed = get_gene_cn_df(genes_list, bin_gene_region_df, impute=True)
               gene_cn_df = gene_cn_df_imputed.drop(['is_imputed'], axis=1)
               n_genes = gene_cn_df.shape[0]
               n_cols = int(np.ceil(n_genes/max_genes))

               gene_list_name = os.path.splitext(f)[0]
               output_png_file_name = os.path.splitext(output.heatmap_cnvs)[0] + '_' + gene_list_name + '.png'

               fig = plt.figure(figsize=(4*n_cols,7), dpi=100)
               outer = gridspec.GridSpec(1, n_cols, wspace=.5, hspace=0.2)
               for i in range(n_cols):
                   plot_heatmap(gene_cn_df.iloc[max_genes*i:max_genes*(i+1),:], title="", colorbar=False, outer=outer[i], fig=fig, max=max_genes,
                                   final=(i==n_cols-1))
               plt.suptitle("Copy number value of gene per subclone", fontsize=20)
               plt.savefig(output_png_file_name, bbox_inches = 'tight',)

               output_csv_file_name = os.path.splitext(output.gene_cn_df)[0] + '_' + gene_list_name + '.csv'
               gene_cn_df_imputed.to_csv(
                   output_csv_file_name
               )

        # Disease specific
        disease_genes_list = pd.read_csv(disease_genes_path, sep="\t", header=None)
        disease_genes_list.columns = ["gene"]
        disease_genes_list = disease_genes_list.gene.values.tolist()

        gene_cn_df_imputed = get_gene_cn_df(disease_genes_list, bin_gene_region_df, impute=True)
        gene_cn_df = gene_cn_df_imputed.drop(['is_imputed'], axis=1)
        n_genes = gene_cn_df.shape[0]
        n_cols = int(np.ceil(n_genes/max_genes))

        fig = plt.figure(figsize=(4*n_cols,7), dpi=100)
        outer = gridspec.GridSpec(1, n_cols, wspace=.5, hspace=0.2)
        for i in range(n_cols):
            plot_heatmap(gene_cn_df.iloc[max_genes*i:max_genes*(i+1),:], title="", colorbar=False, outer=outer[i], fig=fig, max=max_genes,
                            final=(i==n_cols-1))
        plt.suptitle("Copy number value of gene per subclone", fontsize=20)
        plt.savefig(output.heatmap_cnvs, bbox_inches = 'tight',)

        gene_cn_df_imputed.to_csv(
            output.gene_cn_df
        )
