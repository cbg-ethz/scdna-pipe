import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dnapipe
from tqdm import tqdm

sample_name = 'MOTAMUH-T'

# Generate strings to be added to tree file
bin_size = 20000
gene_coordinates_path = f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/emsembl_hg19_genes_simplified.tsv'
chr_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/genomic_coordinates/{sample_name}_scD_Ar1v1.9__chr_stops.tsv'
bin_mask_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/filtering/{sample_name}_scD_Ar1v1.9__excluded_bins.csv'
region_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/breakpoint_detection/{sample_name}_scD_Ar1v1.9_segmented_regions.txt'
region_stops = pd.read_csv(region_stops_path, header=None)
region_stops.columns = ["bin"]
chr_stops = pd.read_csv(chr_stops_path, sep="\t", index_col=1)
bin_mask = pd.read_csv(bin_mask_path, header=None)
all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)

bin_gene_region_df = dnapipe.process_cnvs.get_bin_gene_region_df(bin_size, all_genes, chr_stops, region_stops, bin_mask, priority_genes=None)

# Get first region with chromosome X
def add_genome_str(is_male, bin_gene_region_df):
    regions_with_x = bin_gene_region_df['region'][bin_gene_region_df['chr']=='X'].unique()
    regions_with_y = bin_gene_region_df['region'][bin_gene_region_df['chr']=='Y'].unique()

    autosomal_last_region = [region for region in regions_with_x if region is not None][0]-1
    x_last_region = [region for region in regions_with_y if region is not None][0]-1

    n_regions = region_stops.shape[0]
    string = []
    for i in range(autosomal_last_region+1):
        string.append('{}:2'.format(i))

    if is_male: # add 1 copy of X and 1 copy of Y
        for i in range(autosomal_last_region+1,n_regions+1):
            string.append('{}:1'.format(i))
    else: # add 2 copies of X and 0 copies of Y
        for i in range(autosomal_last_region+1,x_last_region+1):
            string.append('{}:2'.format(i))

    string = '[' + ','.join(string) + ']'
    return string

string = add_genome_str(True, bin_gene_region_df)
print(string)

# Fusion nodes have ID > 1000
assignments_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/_cell_node_ids.tsv'
assignments = np.loadtxt(assignments_path, delimiter='\t')
unique_assignments, unique_assignments_idx, tree_cluster_sizes = np.unique(assignments[:, 1], return_index=True, return_counts=True) # unique IDs

unique_assignments
tree_cluster_sizes

# Number of cells assigned to fusion clusters
np.count_nonzero(assignments > 1000)

inferred_cnvs_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/3nodes_201regions__inferred_cnvs.csv'
inferred_cnvs = np.loadtxt(inferred_cnvs_path, delimiter=',')
unique_cnvs, unique_cnv_idx = np.unique(inferred_cnvs, axis=0, return_index=True)

unique_assignments_idx, unique_cnv_idx

node_ids = assignments[:, 1][unique_assignments_idx]
cnvs = inferred_cnvs[unique_assignments_idx]
cnvs.shape

# Sort clones by distance to diploid profile
dist_to_diploid = []
diploid_profile = np.ones([cnvs.shape[1]]) * 2
for c_id in range(cnvs.shape[0]):
    dist_to_diploid.append(np.linalg.norm(cnvs[c_id]-diploid_profile))
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
gene_coordinates_path = f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/emsembl_hg19_genes_simplified.tsv'
chr_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/genomic_coordinates/{sample_name}_scD_Ar1v1.9__chr_stops.tsv'
region_stops_path = f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/breakpoint_detection/{sample_name}_scD_Ar1v1.9_segmented_regions.txt'
region_stops = pd.read_csv(region_stops_path, header=None)
region_stops.columns = ["bin"]
chr_stops = pd.read_csv(chr_stops_path, sep="\t", index_col=1)
all_genes = pd.read_csv(gene_coordinates_path, sep="\t", index_col=0)
bin_gene_region_df = dnapipe.process_cnvs.get_bin_gene_region_df(bin_size, all_genes, chr_stops, region_stops, bin_mask, cnvs=cnvs_arr, priority_genes=None)

general_main_gene_list = pd.read_csv(f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/general/dna_long_gene_list.txt', sep="\t", header=None)
general_main_gene_list.columns = ["gene"]
general_main_gene_list = general_main_gene_list.gene.values.tolist()
disease_genes_list = pd.read_csv(f'/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/disease_specific/melanoma_genes.txt', sep="\t", header=None)
disease_genes_list.columns = ["gene"]
disease_genes_list = disease_genes_list.gene.values.tolist()
merged_lists = list(set(general_main_gene_list + disease_genes_list))

gene_cn_df = dnapipe.process_cnvs.get_gene_cn_df(merged_lists, bin_gene_region_df, impute=True)

np.unique(merged_lists).shape

# Label subclones for heatmap display
node_ids > 1000
np.where(node_ids > 1000)[0]
gene_cn_df.columns = node_ids.astype(int).astype(str)

# Save copy number values of aberrant subclones
gene_cn_df[gene_cn_df.columns[node_ids > 1000]].to_csv(f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/{sample_name}_fusion_cnvs.csv')

n_subclones = cnvs_arr.shape[0]
offsets = [0.] * n_subclones
if n_subclones > 1:
    offset = 0.1
    offsets = [offset*c for c in range(n_subclones)]
    offsets = offsets - np.mean(offsets)

# Plot clusters
fig = plt.figure(figsize=(10, 4))
for c in range(n_subclones):
    plt.scatter(range(cnvs_arr.shape[1]), cnvs_arr[c] + offsets[c], label='{} ({} cells)'.format(c, tree_cluster_sizes[c]), s=1)
plt.legend(bbox_to_anchor=[1, 1], markerscale=5.)
plt.savefig(f'/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/{sample_name}/singlecell_dna/analysis/cell_fusions/figure.png')
plt.show()
