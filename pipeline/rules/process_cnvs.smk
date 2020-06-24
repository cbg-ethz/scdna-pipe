from scgenpy.process_cnvs.process_cnvs import *

rule create_bin_gene_region_df:
    input:
        gene_coordinates_path = gene_coordinates_path,
        chr_stops_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv",
        excluded_bins_path = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        region_stops_path = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        inferred_cnvs_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__unique_cnvs.csv",
        general_main_gene_list_path = general_main_gene_list_path,
        disease_genes_path = disease_genes_path
    params:
        bin_size = bin_size
    output:
        bin_gene_region_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv"
    benchmark:
        "benchmarks/create_bin_gene_region_df.tsv"
    run:
        genes = pd.read_csv(input.gene_coordinates_path, sep="\t")
        chr_stops = pd.read_csv(input.chr_stops_path, sep="\t", index_col=0)
        excluded_bins = pd.read_csv(input.excluded_bins_path, header=None)
        region_stops = pd.read_csv(input.region_stops_path, header=None)
        region_stops.columns = ["bin"]
        inferred_cnvs = np.loadtxt(input.inferred_cnvs_path, delimiter=',')

        general_main_gene_list = pd.read_csv(input.general_main_gene_list_path, sep="\t", header=None)
        general_main_gene_list.columns = ["gene"]
        general_main_gene_list = general_main_gene_list.gene.values.tolist()

        merged_lists = general_main_gene_list
        try:
            disease_genes_path = os.path.join(gene_lists_path, 'disease_specific', f"{config['disease']}_genes.txt")
            disease_genes_list = pd.read_csv(input.disease_genes_path, sep="\t", header=None)
            disease_genes_list.columns = ["gene"]
            disease_genes_list = disease_genes_list.gene.values.tolist()
            merged_lists = list(set(general_main_gene_list + disease_genes_list))
        except KeyError:
            pass

        df = get_bin_gene_region_df(params.bin_size, genes, chr_stops, region_stops, excluded_bins, cnvs=inferred_cnvs, priority_genes=merged_lists)
        df.to_csv(output.bin_gene_region_df)

rule create_cn_cluster_h5:
    input:
        bin_gene_region_df_path = rules.create_bin_gene_region_df.output.bin_gene_region_df,
        gene_coordinates_path = gene_coordinates_path
    output:
        cn_cluster_h5 = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_cluster.h5"
    benchmark:
        "benchmarks/create_cn_cluster_h5"
    run:
        all_genes = pd.read_csv(input.gene_coordinates_path, sep="\t", index_col=0)
        all_genes_list = all_genes['gene_name'].values.tolist()

        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df_path, index_col=0, low_memory=False)

        gene_cn_df = get_gene_cn_df(all_genes_list, bin_gene_region_df, impute=False)

        cn_cluster_h5 = h5py.File(output.cn_cluster_h5, "w")
        gene_attributes = cn_cluster_h5.create_group("gene_attrs")
        gene_names = np.array(gene_cn_df.index.values, dtype="S16")

        gene_attributes.create_dataset("gene_names", data=gene_names)
        cn_cluster_h5.create_dataset("matrix", data=gene_cn_df.values)
        cn_cluster_h5.close()

rule identify_diploid:
    input:
        inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs.csv",
        excluded_bins_path = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt",
    output:
        is_diploid = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__is_diploid.txt",
    run:
        inferred_cnvs = np.loadtxt(input.inferred_cnvs, delimiter=',')
        excluded_bins = np.loadtxt(input.excluded_bins_path)
        excluded_bins = excluded_bins.astype(bool)
        bin_chr_indicator = np.loadtxt(input.bin_chr_indicator_path)
        bin_chr_indicator = bin_chr_indicator[~excluded_bins].ravel()
        x_start = np.where(bin_chr_indicator=='X')[0]

        inferred_cnvs = inferred_cnvs[:, :x_start]
        n_bins = x_start
        n_cells = inferred_cnvs.shape[0]

        is_diploid = np.array(np.sum(inferred_cnvs, axis=1) / n_bins > 0.98)
        is_diploid = is_diploid.reshape(n_cells, 1)

        np.savetxt(output.is_diploid, is_diploid)

rule assess_cell_fusions:
    input:
        input_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        segmented_counts = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_neutral_states.txt",
        ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs.csv",
        ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids.tsv",
    params:
        inference_binary = os.path.join(scicone_path, "inference"),
        scripts_dir = scripts_dir
    output:
        output_file = os.path.join(analysis_path, "cell_fusions", analysis_prefix) + "__cell_fusions_summary.txt"
    shell:
         "python {params.scripts_dir}/assess_cell_fusions.py -t {input.input_tree}  -d {input.segmented_counts} -r {input.segmented_region_sizes} -n {input.segmented_neutral_states} -i {input.ctree_inferred_cnvs} -c {input.ctree_cell_node_assignments} -b {params.inference_binary} -o {output.output_file}"
