rule visualise_trees:
    params:
        highlight_color = config['secondary_analysis']['highlight_color'],
        max_genes_per_line = config['secondary_analysis']['max_genes_per_line']
    input:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        bin_gene_region_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv",
        genes_to_highlight_path = melanoma_genes
    output:
        cluster_tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.graphviz",
        cluster_tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.png",
        cluster_tree_genes_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.graphviz",
        cluster_tree_genes_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_genes.png"
    run:
        # Tree with region labels
        tree_as_list = tree_to_graphviz(input.cluster_tree)

        for line in tree_as_list:
            print(f"{line}\n")

        with open(output.cluster_tree_graphviz, "w") as file:
            for line in tree_as_list:
                file.write(f"{line}\n")

        try:
            cmd_output = subprocess.run(["dot", "-Tpng", f"{output.cluster_tree_graphviz}", "-o", f"{output.cluster_tree_figure}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        # Tree with gene labels
        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df, index_col=0)
        genes_to_highlight = pd.read_csv(input.genes_to_highlight_path, sep="\t", header=None)
        genes_to_highlight.columns = ["gene"]
        genes_to_highlight = genes_to_highlight.gene.values.tolist()
        tree_genes_as_list = tree_to_graphviz(input.cluster_tree, gene_labels=True, bin_gene_region_df=bin_gene_region_df,
                                genes_to_highlight=genes_to_highlight, highlight_color=params.highlight_color, max_genes_per_line=params.max_genes_per_line)

        for line in tree_genes_as_list:
            print(f"{line}\n")

        with open(output.cluster_tree_genes_graphviz, "w") as file:
            for line in tree_genes_as_list:
                file.write(f"{line}\n")

        try:
            cmd_output = subprocess.run(["dot", "-Tpng:cairo", f"{output.cluster_tree_genes_graphviz}", "-o", f"{output.cluster_tree_genes_figure}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

rule plot_cluster_cnvs:
    params:
    input:
        inferred_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__inferred_cnvs.csv",
        chr_stops = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv"
    output:
        overlapping_cluster_plot = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cluster_profile_overlapping.png"
    benchmark:
        "benchmark/plot_cluster_cnvs.tsv"
    run:
        sa.plot_clusters(input.chr_stops, input.inferred_cnvs)

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
        bin_gene_region_df = pd.read_csv(bin_gene_region_df_path, index_col=0, low_memory=False)

        gene_cn_df = get_gene_cn_df(disease_genes_list, bin_gene_region_df)

        gene_cn_df.to_csv(
            output.gene_cn_df
        )

        figure_width = gene_cn_df.shape[0] / 2 + 1.5
        plt.figure(figsize=(8, figure_width))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        heatmap = sns.heatmap(
            gene_cn_df,
            annot=True,
            cmap=cmap,
            vmin=0,
            vmax=4,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"ticks": [0, 1, 2, 3, 4]},
        )
        heatmap.set_title("Copy number values of genes per cluster")
        heatmap.set_facecolor("#656565")
        heatmap = heatmap.get_figure()
        heatmap.savefig(
            output.heatmap_cnvs
        )
        plt.close()

rule plot_roche_heatmap:
    params:
        roche_genes_path = roche_genes_path,
        bin_size = bin_size
    input:
        cn_cluster_h5 = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_cluster.h5"
    output:
        roche_gene_cn_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__roche_cn_gene_df.csv",
        roche_heatmap_cnvs = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__roche_heatmap_cnvs.png"
    benchmark:
        "benchmark/plot_roche_heatmap.tsv"
    run:
        gene_list = pd.read_csv(params.roche_genes_path, sep="\t", header=None)
        gene_list.columns = ["gene"]
        gene_list = gene_list.gene.values.tolist()

        cn_cluster_h5 = h5py.File(input.cn_cluster_h5, "r")
        cnvs = cn_cluster_h5['matrix']
        gene_names = cn_cluster_h5['gene_attrs']['gene_names'][()].astype(str)

        gene_cn_df = pd.DataFrame(columns=gene_list)
        for gene in gene_list:
            gene_cn_df[gene] = cnvs[np.where(gene_names==gene)[0][0]]
        gene_cn_df = gene_cn_df.T

        gene_cn_df.to_csv(
            output.roche_gene_cn_df
        )

        figure_width = gene_cn_df.shape[0] / 2 + 1.5
        plt.figure(figsize=(8, figure_width))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        heatmap = sns.heatmap(
            gene_cn_df,
            annot=True,
            cmap=cmap,
            vmin=0,
            vmax=4,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"ticks": [0, 1, 2, 3, 4]},
        )
        heatmap.set_title("Copy number values of genes (in Roche's list) per cluster")
        heatmap.set_facecolor("#656565")
        heatmap = heatmap.get_figure()
        heatmap.savefig(
            output.roche_heatmap_cnvs
        )
        plt.close()
