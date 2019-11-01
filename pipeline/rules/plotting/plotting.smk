from plotting import *

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
        cell_node_assignments = pd.read_csv(input.cell_node_assignments, sep='\t', index_col=0, header=None)
        keys, values = [arr.tolist() for arr in np.unique(cell_node_assignments[1].values, return_counts=True)]
        node_sizes = dict(zip(np.array(keys).astype(str).tolist(), values))

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
        add_offsets = config['plotting']['profiles']['add_offsets'],
        s = config['plotting']['profiles']['s']
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

        add_offsets=False
        if params.add_offsets=='True':
            add_offsets=True

        output_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix)
        plot_cluster_cnvs(cnvs_arr, chr_stops, output_path=output_path, add_offsets=add_offsets, s=params.s)

rule plot_heatmap:
    input:
        disease_genes_path = disease_genes_path,
        bin_gene_region_df_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv"
    output:
        gene_cn_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_gene_df.csv",
        heatmap_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__heatmap_cnvs.png"
    benchmark:
        "benchmark/plot_heatmap.tsv"
    run:
        disease_genes_list = pd.read_csv(disease_genes_path, sep="\t", header=None)
        disease_genes_list.columns = ["gene"]
        disease_genes_list = disease_genes_list.gene.values.tolist()
        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df_path, index_col=0, low_memory=False)

        gene_cn_df_imputed = get_gene_cn_df(disease_genes_list, bin_gene_region_df, impute=True)
        gene_cn_df = get_gene_cn_df(disease_genes_list, bin_gene_region_df, impute=False)

        is_imputed = np.array(np.isnan(gene_cn_df))

        plot_heatmap(gene_cn_df_imputed, is_imputed=is_imputed, output_path=output.heatmap_cnvs)

        gene_cn_df_imputed.to_csv(
            output.gene_cn_df
        )

rule create_cluster_plots_list:
    params:
        analysis_prefix = analysis_prefix,
        analysis_path = analysis_path
    input:
        secondary_analysis_done = "secondary_analysis_done.txt"
    output:
        cluster_plots_list = os.path.join(analysis_path, "clustering", analysis_prefix ) + "__cluster_profile_files.txt"
    shell:
        "cd {params.analysis_path}/clustering; ls | egrep '{params.analysis_prefix}__cluster_profile_..?\.png' > {params.analysis_prefix}__cluster_profile_files.txt"

rule create_heatmap_plots_list:
    params:
        analysis_prefix = analysis_prefix,
        analysis_path = analysis_path
    input:
        secondary_analysis_done = "secondary_analysis_done.txt"
    output:
        heatmap_plots_list = os.path.join(analysis_path, "clustering", analysis_prefix ) + "__cn_genes_clusters_files.txt"
    shell:
        "cd {params.analysis_path}/clustering; ls | egrep '{params.analysis_prefix}__cn_genes_clusters_chr..?_heatmap\.png' > {params.analysis_prefix}__cn_genes_clusters_files.txt"
