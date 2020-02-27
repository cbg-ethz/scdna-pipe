from scgenpy.plotting.plotting import *
from scgenpy.process_cnvs.process_cnvs import *
from scgenpy.preprocessing.utils import sort_chromosomes

rule plot_cnv_matrix:
    input:
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv",
        segmented_regions = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        excluded_bins = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt",
        inferred_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__unique_cluster_tree_cnvs.csv"
    output:
        inferred_cnvs_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_cnvs_bins.png",
        sorted_normalised_counts_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_normalised_counts_bins.png",
        sorted_inferred_cnvs_heatmap = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_tree_sorted_cnvs_bins.png"
    benchmark:
        "benchmark/plot_cnv_matrix.tsv"
    run:
        from scipy.cluster.hierarchy import ward, leaves_list
        from scipy.spatial.distance import pdist

        print("loading the normalised bins...")
        normalised_bins = np.loadtxt(input.normalised_bins, delimiter=',')
        print(f"shape of normalised bins: {normalised_bins.shape}")

        print("loading the normalised regions...")
        normalised_regions = np.loadtxt(input.normalised_regions, delimiter=',')
        print(f"shape of normalised regions: {normalised_regions.shape}")

        Z = ward(pdist(normalised_regions))
        hclust_index = leaves_list(Z)

        print("loading the inferred cnvs...")
        inferred_cnvs = np.loadtxt(input.inferred_cnvs, delimiter=',')
        print(f"shape of inferred cnvs: {inferred_cnvs.shape}")

        # Annotate with chromosome coordinates
        chr_indicator = pd.read_csv(input.bin_chr_indicator, header=None).values.ravel()

        # Remove excluded bins
        excluded_bins = pd.read_csv(input.excluded_bins, header=None).values.ravel().astype(bool)
        chr_indicator_filtered = np.array(chr_indicator)[~excluded_bins]
        chromosome_list = utils.sort_chromosomes(np.unique(chr_indicator))
        chr_stops_filtered_bins = [np.where(chr_indicator_filtered==chr)[0][-1] for chr in chromosome_list]
        chr_stops_filtered_bins_dict = dict(zip(chromosome_list, chr_stops_filtered_bins))

        print("plotting the heatmap of clustered inferred cnvs ordered according to input...")
        plot_bins(inferred_cnvs[hclust_index], chr_stops_filtered_bins_dict, vlines=True,
                    cbar_title='Copy number\nstates', vmax=None,
                    figsize=(24,8), output_path=output.inferred_cnvs_heatmap, dpi=300)

        # Now plot data with the inferred CNV information (i.e. cells ordered and annotated by inferred clone)
        clustered_inferred_cnvs, clustered_normalised_bins, clustered_labels = cluster_clones(inferred_cnvs, normalised_bins, normalised_regions)

        print("plotting the heatmap of normalised counts ordered by clone...")
        plot_bins(clustered_normalised_bins, chr_stops_filtered_bins_dict, vlines=True,
                    cbar_title='Normalised\ncounts', vmax=2, annotations=clustered_labels,
                    figsize=(24,8), output_path=output.sorted_normalised_counts_heatmap, dpi=300)

        print("plotting the heatmap of inferred cnvs ordered by clone...")
        plot_bins(clustered_inferred_cnvs, chr_stops_filtered_bins_dict, vlines=True,
                    cbar_title='Copy number\nstates', vmax=None, annotations=clustered_labels,
                    figsize=(24,8), output_path=output.sorted_inferred_cnvs_heatmap, dpi=300)


rule plot_data_matrix:
    input:
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv",
        segmented_regions = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        excluded_bins = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt",
    output:
        normalised_bins_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins.png",
        normalised_bins_clustered_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered.png",
        normalised_bins_clustered_bps_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered_bps.png"
    benchmark:
        "benchmark/plot_breakpoints.tsv"
    run:
        from scipy.cluster.hierarchy import ward, leaves_list
        from scipy.spatial.distance import pdist

        print("loading the normalised bins...")
        normalised_bins = np.loadtxt(input.normalised_bins, delimiter=',')
        print(f"shape of normalised bins: {normalised_bins.shape}")

        print("loading the normalised regions...")
        normalised_regions = np.loadtxt(input.normalised_regions, delimiter=',')
        print(f"shape of normalised regions: {normalised_regions.shape}")

        bps = np.loadtxt(input.segmented_regions)
        print(f"number of breakpoints: {len(bps)}")

        # Annotate with chromosome coordinates
        chr_indicator = pd.read_csv(input.bin_chr_indicator, header=None).values.ravel()
        excluded_bins = pd.read_csv(input.excluded_bins, header=None).values.ravel().astype(bool)
        chr_indicator_filtered = np.array(chr_indicator)[~excluded_bins]
        chromosome_list = utils.sort_chromosomes(np.unique(chr_indicator))
        chr_stops_filtered_bins = [np.where(chr_indicator_filtered==chr)[0][-1] for chr in chromosome_list]
        chr_stops_filtered_bins_dict = dict(zip(chromosome_list, chr_stops_filtered_bins))

        Z = ward(pdist(normalised_regions))
        hclust_index = leaves_list(Z)

        print("plotting the heatmap ...")
        plot_bins(normalised_bins, chr_stops_filtered_bins_dict, vlines=True,
                    cbar_title='Normalised\ncounts',
                    figsize=(24,8), output_path=output.normalised_bins_heatmap)

        print("plotting the clustered heatmap ...")
        plot_bins(normalised_bins[hclust_index], chr_stops_filtered_bins_dict, vlines=True,
                    cbar_title='Normalised\ncounts',
                    figsize=(24,8), output_path=output.normalised_bins_clustered_heatmap)

        print("plotting the clustered heatmap with breakpoints ...")
        plot_bins(normalised_bins[hclust_index], chr_stops_filtered_bins_dict, vlines=False, bps=bps,
                    cbar_title='Normalised\ncounts',
                    figsize=(24,8), output_path=output.normalised_bins_clustered_bps_heatmap)

rule visualise_trees:
    params:
        highlight_color = config['plotting']['trees']['highlight_color'],
        max_genes_per_line = config['plotting']['trees']['max_genes_per_line']
    input:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids.tsv",
        bin_gene_region_df_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv",
        genes_to_highlight_path = disease_genes_path
    output:
        cluster_tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.graphviz",
        cluster_tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.png",
        cluster_tree_genes_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.graphviz",
        cluster_tree_genes_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.png"
    run:
        # Node sizes
        cell_node_assignments = np.loadtxt(input.cell_node_assignments, delimiter='\t')
        keys, values = [arr.tolist() for arr in np.unique(cell_node_assignments[:,1], return_counts=True)]
        node_sizes = dict(zip(np.array(keys).astype(int).astype(str).tolist(), values))

        # Tree with region labels
        regions_tree = tree_to_graphviz(input.cluster_tree, node_sizes=node_sizes)
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
                                genes_to_highlight=genes_to_highlight, highlight_color=params.highlight_color,
                                max_genes_per_line=params.max_genes_per_line)
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
        inferred_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__inferred_cnvs.csv",
        chr_stops = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv"
    output:
        overlapping_cluster_plot = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_profile_overlapping.png"
    benchmark:
        "benchmark/plot_cluster_cnvs.tsv"
    run:
        cnvs_arr = np.loadtxt(input.inferred_cnvs, delimiter=',')
        chr_stops = pd.read_csv(input.chr_stops, sep="\t")

        ymax = params.ploidy + params.max_amp_val + 1

        output_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix)
        plot_cluster_cnvs(cnvs_arr, chr_stops, output_path=output_path, ymax=ymax, offset_sizes=params.offset_sizes, s=params.s)

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
        "benchmark/plot_heatmap.tsv"
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

               gene_list_name = os.path.splitext(f)[0]
               output_png_file_name = os.path.splitext(output.heatmap_cnvs)[0] + '_' + gene_list_name + '.png'
               plot_heatmap(gene_cn_df_imputed, output_path=output_png_file_name)

               output_csv_file_name = os.path.splitext(output.gene_cn_df)[0] + '_' + gene_list_name + '.csv'
               gene_cn_df_imputed.to_csv(
                   output_csv_file_name
               )

        # Disease specific
        disease_genes_list = pd.read_csv(disease_genes_path, sep="\t", header=None)
        disease_genes_list.columns = ["gene"]
        disease_genes_list = disease_genes_list.gene.values.tolist()

        gene_cn_df_imputed = get_gene_cn_df(disease_genes_list, bin_gene_region_df, impute=True)

        plot_heatmap(gene_cn_df_imputed, output_path=output.heatmap_cnvs)

        gene_cn_df_imputed.to_csv(
            output.gene_cn_df
        )
