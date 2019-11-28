import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dnapipe
from tqdm import tqdm

sample_name = "MOTAMUH-T"

# Generate strings to be added to tree file
bin_size = 20000
gene_coordinates_path = f"/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/emsembl_hg19_genes_simplified.tsv"
chr_stops_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/genomic_coordinates/{sample_name}_scD_Ar1v1.9__chr_stops.tsv"
bin_mask_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/filtering/{sample_name}_scD_Ar1v1.9__excluded_bins.csv"
region_stops_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/breakpoint_detection/{sample_name}_scD_Ar1v1.9_segmented_regions.txt"
region_stops = pd.read_csv(region_stops_path, header=None)
region_stops.columns = ["bin"]
chr_stops = pd.read_csv(chr_stops_path, sep="\t", index_col=1)
bin_mask = pd.read_csv(bin_mask_path, header=None)
all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)

bin_gene_region_df = dnapipe.process_cnvs.get_bin_gene_region_df(
    bin_size, all_genes, chr_stops, region_stops, bin_mask, priority_genes=None
)

# Get first region with chromosome X
def add_genome_str(is_male, bin_gene_region_df):
    regions_with_x = bin_gene_region_df["region"][
        bin_gene_region_df["chr"] == "X"
    ].unique()
    regions_with_y = bin_gene_region_df["region"][
        bin_gene_region_df["chr"] == "Y"
    ].unique()

    autosomal_last_region = [region for region in regions_with_x if region is not None][
        0
    ] - 1
    x_last_region = [region for region in regions_with_y if region is not None][0] - 1

    n_regions = region_stops.shape[0]
    string = []
    for i in range(autosomal_last_region + 1):
        string.append("{}:2".format(i))

    if is_male:  # add 1 copy of X and 1 copy of Y
        for i in range(autosomal_last_region + 1, n_regions + 1):
            string.append("{}:1".format(i))
    else:  # add 2 copies of X and 0 copies of Y
        for i in range(autosomal_last_region + 1, x_last_region + 1):
            string.append("{}:2".format(i))

    string = "[" + ",".join(string) + "]"
    return string


string = add_genome_str(True, bin_gene_region_df)
print(string) # add this event vector to fusion nodes in txt file of tree

# Fusion nodes have ID > 1000
assignments_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/_cell_node_ids.tsv"
assignments = np.loadtxt(assignments_path, delimiter="\t")
unique_assignments, unique_assignments_idx, tree_cluster_sizes = np.unique(
    assignments[:, 1], return_index=True, return_counts=True
)  # unique IDs

unique_assignments
tree_cluster_sizes

# Number of cells assigned to fusion clusters
np.count_nonzero(assignments > 1000)

inferred_cnvs_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/3nodes_201regions__inferred_cnvs.csv"
inferred_cnvs = np.loadtxt(inferred_cnvs_path, delimiter=",")
unique_cnvs, unique_cnv_idx = np.unique(inferred_cnvs, axis=0, return_index=True)

unique_assignments_idx, unique_cnv_idx

node_ids = assignments[:, 1][unique_assignments_idx]
cnvs = inferred_cnvs[unique_assignments_idx]
cnvs.shape

# Sort clones by distance to diploid profile
dist_to_diploid = []
diploid_profile = np.ones([cnvs.shape[1]]) * 2
for c_id in range(cnvs.shape[0]):
    dist_to_diploid.append(np.linalg.norm(cnvs[c_id] - diploid_profile))
order = np.argsort(dist_to_diploid)
cnvs = cnvs[order]
node_ids = node_ids[order]

node_ids

# Get genes
def add_filtered_bins_back(unique_cnvs, bin_mask):
    """
    Adds the filtered bins back to the inferred cnvs
    :param unique_cnvs_path: path to the file containing unique copy number profiles
    :param bin_mask_path: path to the excluded bins file
    :return:
    """
    cnvs_mat = []
    cnvs_counter = 0

    for bin_idx, bin_val in enumerate(bin_mask[0]):
        c_row = []
        if bin_val:
            for c_id in range(unique_cnvs.shape[0]):
                c_row.append(None)
        else:
            for c_id in range(unique_cnvs.shape[0]):
                c_row.append(unique_cnvs[c_id][cnvs_counter])
            cnvs_counter += 1
        cnvs_mat.append(c_row)

    cnvs_arr = np.array(cnvs_mat, dtype=float).T
    print(f"cnvs_arr shape: {cnvs_arr.shape}")

    return cnvs_arr


cnvs_arr = add_filtered_bins_back(cnvs, bin_mask)

bin_size = 20000
gene_coordinates_path = f"/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/emsembl_hg19_genes_simplified.tsv"
chr_stops_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/genomic_coordinates/{sample_name}_scD_Ar1v1.9__chr_stops.tsv"
region_stops_path = f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/breakpoint_detection/{sample_name}_scD_Ar1v1.9_segmented_regions.txt"
region_stops = pd.read_csv(region_stops_path, header=None)
region_stops.columns = ["bin"]
chr_stops = pd.read_csv(chr_stops_path, sep="\t", index_col=1)
all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)
bin_gene_region_df = dnapipe.process_cnvs.get_bin_gene_region_df(
    bin_size,
    all_genes,
    chr_stops,
    region_stops,
    bin_mask,
    cnvs=cnvs_arr,
    priority_genes=None,
)

general_main_gene_list = pd.read_csv(
    f"/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/general/dna_long_gene_list.txt",
    sep="\t",
    header=None,
)
general_main_gene_list.columns = ["gene"]
general_main_gene_list = general_main_gene_list.gene.values.tolist()
disease_genes_list = pd.read_csv(
    f"/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/disease_specific/melanoma_genes.txt",
    sep="\t",
    header=None,
)
disease_genes_list.columns = ["gene"]
disease_genes_list = disease_genes_list.gene.values.tolist()
merged_lists = list(set(general_main_gene_list + disease_genes_list))

gene_cn_df = dnapipe.process_cnvs.get_gene_cn_df(
    merged_lists, bin_gene_region_df, impute=True
)

np.unique(merged_lists).shape

# Label subclones for heatmap display
node_ids > 1000
np.where(node_ids > 1000)[0]
gene_cn_df.columns = node_ids.astype(int).astype(str)

# Save copy number values of aberrant subclones
gene_cn_df[gene_cn_df.columns[node_ids > 1000]].to_csv(
    f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/{sample_name}_fusion_cnvs.csv"
)

n_subclones = cnvs_arr.shape[0]
offsets = [0.0] * n_subclones
if n_subclones > 1:
    offset = 0.1
    offsets = [offset * c for c in range(n_subclones)]
    offsets = offsets - np.mean(offsets)

# Plot clusters
fig = plt.figure(figsize=(10, 4))
for c in range(n_subclones):
<<<<<<< HEAD
    plt.scatter(range(cnvs_arr.shape[1]), cnvs_arr[c] + offsets[c], label='{} ({} cells)'.format(c, tree_cluster_sizes[c]), s=1)
plt.legend(bbox_to_anchor=[1, 1], markerscale=5.)
# plt.savefig(f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/figure.png')
=======
    plt.scatter(
        range(cnvs_arr.shape[1]),
        cnvs_arr[c] + offsets[c],
        label="{} ({} cells)".format(c, tree_cluster_sizes[c]),
        s=1,
    )
plt.legend(bbox_to_anchor=[1, 1], markerscale=5.0)
plt.savefig(
    f"/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/figure.png"
)
>>>>>>> 6a2c4ab33564c95cbb7e0a0005bcef0b36d7f336
plt.show()

sample_name = 'MADIBUG-T'
# Create heatmap
sample_names = ['MYNESYB-T', 'MOTAMUH-T', 'MEKOBAB-T', 'MOVAZYQ-T', 'MOQAVIJ-T', 'MODUDOL-T', 'MOBUBOT-T', 'MISYPUP-T', 'MEXUXEH-T', 'MEVIXYV-T', 'MEHYLOB-T', 'MANOFYB-T', 'MAJOFIJ-T', 'MAHACEB-T', 'MADIBUG-T']

for sample_name in sample_names:
    # Generate strings to be added to tree file
    bin_size = 20000
    gene_coordinates_path = f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/emsembl_hg19_genes_simplified.tsv'
    chr_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/genomic_coordinates/{sample_name}_scD_Ar1v1.9__chr_stops.tsv'
    bin_mask_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/filtering/{sample_name}_scD_Ar1v1.9__excluded_bins.csv'
    region_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/breakpoint_detection/{sample_name}_scD_Ar1v1.9_segmented_regions.txt'
    inferred_cnvs_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/inferred_cnvs/{sample_name}_scD_Ar1v1.9__inferred_cnvs.csv'
    inferred_cnvs = np.loadtxt(inferred_cnvs_path, delimiter=',')
    region_stops = pd.read_csv(region_stops_path, header=None)
    region_stops.columns = ["bin"]
    chr_stops = pd.read_csv(chr_stops_path, sep="\t", index_col=1)
    bin_mask = pd.read_csv(bin_mask_path, header=None)
    all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)

    bin_gene_region_df = dnapipe.process_cnvs.get_bin_gene_region_df(bin_size, all_genes, chr_stops, region_stops, bin_mask, cnvs=inferred_cnvs, priority_genes=None)

    dnapipe.process_cnvs.get_region_with_gene('HLA-A', bin_gene_region_df)

    roche_gene_list = pd.read_csv(f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/general/roche_gene_list.txt', sep="\t", header=None)
    roche_gene_list.columns = ["gene"]
    roche_gene_list = roche_gene_list.gene.values.tolist()

    gene_cn_df = get_gene_cn_df(roche_gene_list, bin_gene_region_df, impute=True)

    tree_cluster_sizes_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/tree_learning/{sample_name}_scD_Ar1v1.9__tree_cluster_sizes.csv'
    tree_cluster_sizes = pd.read_csv(tree_cluster_sizes_path, header=None)[0].values.tolist()
    gene_cn_df.columns = [f'subclone {i} ({tree_cluster_sizes[i]} cells)' for i, name in enumerate(gene_cn_df.columns.astype(str).values.tolist()[:-1])] + ['is_imputed']

    gene_cn_df.to_csv(f'/cluster/work/bewi/members/pedrof/for_petra/{sample_name}_cnvs.csv')

dnapipe.process_cnvs.get_region_with_gene('HLA-C', bin_gene_region_df)
bin_gene_region_df.iloc[:50]
def get_surrounding_regions(gene, bin_gene_region_df):
    """
        Finds the regions to the left and to the right of a gene's bins
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: Tuple containing (left_region, right_region)
    """
    bin_gene_region_df = bin_gene_region_df.copy(deep=True)
    bin_gene_region_df['gene'] = bin_gene_region_df['gene'].astype(str).apply(lambda x: x.split(',')).apply(lambda x: set(x))
    bins = bin_gene_region_df[bin_gene_region_df['gene'].apply(lambda x: gene in x)].index.values

    left_bin = bins[0]
    right_bin = bins[-1]

    left_regions = bin_gene_region_df['region'].values[:left_bin]
    left_regions = np.array([x for x in left_regions if x is not None])
    left_region = left_regions[~np.isnan(left_regions)][-1]

    right_regions = bin_gene_region_df['region'].values[right_bin+1:]
    right_regions = np.array([x for x in right_regions if x is not None])
    right_region = right_regions[~np.isnan(right_regions)][0]

    return (int(left_region), int(right_region))

def get_gene_cn_df(gene_list, bin_gene_region_df, impute=False):
    """
        Creates and returns the dataframe of copy numbers, genes by cluster ids
        :param gene_list: the input list of genes to be specified
        :param bin_gene_region_df: pd.DataFrame with cnv per bin, region and gene
        :param impute: replace NaN values with median of left and rightmost regions
        :return: CN dataframe of genes by subclone
    """
    n_subclones = np.count_nonzero(['cnv' in column for column in bin_gene_region_df.columns])
    cluster_ids = range(n_subclones)
    gene_cn_df = pd.DataFrame(index=cluster_ids)

    df = bin_gene_region_df.copy(deep=True)
    df['gene'] = df['gene'].astype(str).apply(lambda x: x.split(',')).apply(lambda x: set(x))

    is_imputed = np.empty(len(gene_list))

    # for each gene
    i = -1
    for gene in tqdm(gene_list):
        i += 1
        gene_cn_per_cluster = []
        is_imputed[i] = False

        for c_id in cluster_ids:
            bins = df[df['gene'].apply(lambda x: gene in x)].index.values
            median_cn = np.nanmedian(
                bin_gene_region_df['cnv_{}'.format(c_id)][bins].values
            )
            # If NaN, impute with median value of regions surrounding it
            if np.isnan(median_cn) and impute:
                if len(bins) > 0:
                    # Get regions surrounding gene
                    left_region, right_region = get_surrounding_regions(gene, bin_gene_region_df)

                    # get CNV values of region surronding gene
                    left_cn = bin_gene_region_df['cnv_{}'.format(c_id)][bin_gene_region_df['region']==left_region].iloc[0]
                    right_cn = bin_gene_region_df['cnv_{}'.format(c_id)][bin_gene_region_df['region']==right_region].iloc[0]

                    median_cn = np.nanmedian([left_cn, right_cn])

                    is_imputed[i] = True
                else:
                    print(f'Gene {gene} does not exist.')

            if not np.isnan(median_cn):
                if median_cn > 2:
                    median_cn = int(np.floor(median_cn))
                elif median_cn < 2:
                    median_cn = int(np.ceil(median_cn))

            gene_cn_per_cluster.append(median_cn)

        gene_cn_df[gene] = gene_cn_per_cluster

    print("Transposing the dataframe...")
    gene_cn_df = gene_cn_df.T
    if impute:
        gene_cn_df['is_imputed'] = is_imputed.tolist()
    # gene_cn_df = gene_cn_df.rename(columns = {'two':'new_name'})
    print("Sorting the genes...")
    gene_cn_df.sort_index(inplace=True)

    return gene_cn_df

gene_cn_df = get_gene_cn_df(roche_gene_list, bin_gene_region_df, impute=False)

gene_cn_df


import glob, os
import pandas as pd

from pandas import DataFrame, ExcelWriter

writer = ExcelWriter("/Users/pedrof/projects/for_petra/new.xlsx")

for filename in glob.glob("/Users/pedrof/projects/for_petra/*.csv"):
    df_csv = pd.read_csv(filename)

    (_, f_name) = os.path.split(filename)
    (f_shortname, _) = os.path.splitext(f_name)

    df_csv.to_excel(writer, f_shortname, index=False)

writer.save()
