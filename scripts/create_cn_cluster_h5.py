import argparse
from scgenpy.process_cnvs.process_cnvs import *
import h5py
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-d", "--bin_gene_region_df_path", required=True,
)
parser.add_argument(
    "-a",
    "--annotations",
    required=False,
    default="/cluster/work/bewi/members/pedrof/tupro_code/scdna-pipe/required_files/genes_of_interest/ensembl_hg19_annotations.tsv",
)

args = parser.parse_args()

bin_gene_region_df_path = args.bin_gene_region_df_path
gene_coordinates_path = args.annotations
output_file = bin_gene_region_df_path.split("__")[0] + "__cn_cluster.h5"

all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)
all_genes_list = all_genes["gene_name"].values.tolist()

bin_gene_region_df = pd.read_csv(bin_gene_region_df_path, index_col=0, low_memory=False)

gene_cn_df = get_gene_cn_df(all_genes_list, bin_gene_region_df, impute=False)

cn_cluster_h5 = h5py.File(output_file, "w")
gene_attributes = cn_cluster_h5.create_group("gene_attrs")
gene_names = np.array(gene_cn_df.index.values, dtype="S16")

gene_attributes.create_dataset("gene_names", data=gene_names)
cn_cluster_h5.create_dataset("matrix", data=gene_cn_df.values)
cn_cluster_h5.close()
