from process_cnvs import *

rule create_bin_gene_region_df:
    input:
        gene_coordinates_path = gene_coordinates_path,
        chr_stops_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv",
        excluded_bins_path = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        region_stops_path = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        inferred_cnvs_path = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__inferred_cnvs.csv",
        priority_genes_path = priority_genes_path
    params:
        bin_size = bin_size
    output:
        bin_gene_region_df = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__bin_gene_region_df.csv"
    benchmark:
        "benchmark/create_bin_gene_region_df.tsv"
    run:
        genes = pd.read_csv(input.gene_coordinates_path, sep="\t")
        chr_stops = pd.read_csv(input.chr_stops_path, sep="\t", index_col=1)
        excluded_bins = pd.read_csv(input.excluded_bins_path, header=None)
        region_stops = pd.read_csv(input.region_stops_path, header=None)
        region_stops.columns = ["bin"]
        priority_genes = pd.read_csv(input.priority_genes_path, sep="\t", header=None)
        priority_genes.columns = ["gene"]
        priority_genes = priority_genes.gene.values.tolist()
        inferred_cnvs = np.loadtxt(input.inferred_cnvs_path, delimiter=',')

        df = get_bin_gene_region_df(params.bin_size, genes, chr_stops, region_stops, excluded_bins, cnvs=inferred_cnvs, priority_genes=priority_genes)
        df.to_csv(output.bin_gene_region_df)

rule create_cn_cluster_h5:
    input:
        bin_gene_region_df_path = rules.create_bin_gene_region_df.output.bin_gene_region_df,
        gene_coordinates_path = gene_coordinates_path
    output:
        cn_cluster_h5 = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "__cn_cluster.h5"
    benchmark:
        "benchmark/create_cn_cluster_h5"
    run:
        all_genes = pd.read_csv(input.gene_coordinates_path, sep="\t", index_col=0)
        all_genes_list = all_genes['Gene name'].values.tolist()

        bin_gene_region_df = pd.read_csv(input.bin_gene_region_df_path, index_col=0, low_memory=False)

        gene_cn_df = get_gene_cn_df(all_genes_list, bin_gene_region_df, impute=False)

        cn_cluster_h5 = h5py.File(output.cn_cluster_h5, "w")
        gene_attributes = cn_cluster_h5.create_group("gene_attrs")
        gene_names = np.array(gene_cn_df.index.values, dtype="S16")

        gene_attributes.create_dataset("gene_names", data=gene_names)
        cn_cluster_h5.create_dataset("matrix", data=gene_cn_df.values)
        cn_cluster_h5.close()
