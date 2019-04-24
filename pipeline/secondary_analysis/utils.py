import numpy as np
import h5py

def merge_chromosomes(h5, key='normalized_counts'):

    n_cells = h5['cell_barcodes'][:].shape[0]
    all_chromosomes = list(h5[key].keys())
    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(h5[key][chr][:][0:n_cells,:]) # select only the cells, not cell groups

    cell_all_chrs = np.concatenate(cnv_matrices, axis=1)
    return cell_all_chrs
