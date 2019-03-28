import numpy as np
import pandas as pd
import h5py
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c","--cell_assignment",required=True, help="the cluster assignments of cells")
parser.add_argument("-o","--output_path",required=False, default="./", help="path to the output")
parser.add_argument("-s", "--sample_name",required=False, default="", help="name of the sample")
parser.add_argument("-g", "--genes", required=True, help="genes of interest and their positions")
parser.add_argument("-h5", "--cnv_hdf5", required=True, help="cnv_data.h5 file")


args = parser.parse_args()


cell_assignment = pd.read_csv(args.cell_assignment,sep='\t')
melanoma_genes = pd.read_csv(args.genes,sep='\t')
cnv_data = h5py.File(args.cnv_hdf5)

n_cells = cell_assignment.shape[0]
print("n_cells: ", n_cells)

bin_size = cnv_data['constants']['bin_size'][()]
print("bin_size: ", bin_size)

cluster_ids = sorted(list(set(cell_assignment.cluster)))
print("cluster_ids: ", cluster_ids)

gene_cn_df = pd.DataFrame(index=cluster_ids)
# for each gene
for index, row in melanoma_genes.iterrows():
    start_bin = int(row['Gene start (bp)']/bin_size)
    stop_bin = int(row['Gene end (bp)']/bin_size)
    chromosome = str(row['Chromosome/scaffold name'])
    gene_name = row['Gene name']
    #print(start_bin, stop_bin)
    mean_copy_numbers = []
    for c_id in cluster_ids:
        # get all the cells that belong to that cluster
        cells = cell_assignment[cell_assignment.cluster == c_id]['cell_barcode'].values.tolist()
        cn_states = cnv_data['cnvs'][chromosome][:n_cells][cells,start_bin:stop_bin+1]

        # -127 means imputed to 0
        cn_states[cn_states == -127] = 0
        cn_states[cn_states == -128] = 129
        cn_states = np.abs(cn_states)
        cn_states = cn_states.astype('float')
        cn_states[cn_states == 129] = np.nan

        min_cn_cell_bin = np.nanmin(cn_states,axis=1)
        avg_cn_cell = np.nanmean(min_cn_cell_bin)
        mean_copy_numbers.append(avg_cn_cell)

    #print(mean_copy_numbers)
    gene_cn_df[gene_name] = mean_copy_numbers

gene_cn_df.T.to_csv(args.output_path + '/' + args.sample_name + '__cn_gene_cluster.tsv', sep='\t')
heatmap = sns.heatmap(gene_cn_df.T, annot=True, cmap='RdBu',vmin=0,vmax=4).get_figure()
heatmap.savefig(args.output_path + '/' + args.sample_name + "__cn_genes_clusters_heatmap.png")
