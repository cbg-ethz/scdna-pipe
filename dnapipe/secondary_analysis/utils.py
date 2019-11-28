import numpy as np
import h5py
import re


def merge_chromosomes(h5, key="normalized_counts"):

    n_cells = h5["cell_barcodes"][:].shape[0]
    all_chromosomes = list(h5[key].keys())
    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(
            h5[key][chr][:][0:n_cells, :]
        )  # select only the cells, not cell groups

    cell_all_chrs = np.concatenate(cnv_matrices, axis=1)
    return cell_all_chrs


def get_tree_scores(tree_paths):
    """
        Creates the list of tree scores by parsing the tree files
        :param trees_path: a list of tree paths
        :return: the list of tree scores
    """
    tree_scores = []
    for tree_path in tree_paths:
        with open(tree_path) as f:
            list_tree_file = list(f)

        for line in list_tree_file:
            if line.startswith("Tree score:"):
                score = line.rstrip("\n").lstrip("Tree score:").lstrip(" ")
                tree_scores.append(float(score))

    return tree_scores


def rename_fastq(s_name):
    """
        renames the merged fastqs according to the bcl2fastq naming convention
        Sample input name: MERGED_BSSE_QGF_123456_ZXVN2SHG5_1_QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz
    """
    split_name = s_name.split("_")
    new_name = "_".join(split_name[6:7] + split_name[-4:])
    return new_name
