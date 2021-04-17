import numpy as np
import h5py
import re
import difflib


def sort_chromosomes(chromosome_list):
    """
    Sorts a list of unordered chromosome names
    :param chromosome_list: list of unordered characters denoting chromosomes '1', '2', ..., 'X', 'Y'
    """
    # Replace X and Y with 23 and 24
    sorted_chromosome_list = np.array(chromosome_list)
    sorted_chromosome_list[np.where(sorted_chromosome_list == "X")[0]] = '23'
    sorted_chromosome_list[np.where(sorted_chromosome_list == "Y")[0]] = '24'

    # Sort
    sorted_chromosome_list = sorted([int(x) for x in sorted_chromosome_list])
    sorted_chromosome_list = [str(x) for x in sorted_chromosome_list]

    # Convert back to string
    sorted_chromosome_list[sorted_chromosome_list.index("23")] = "X"
    sorted_chromosome_list[sorted_chromosome_list.index("24")] = "Y"

    return sorted_chromosome_list


def merge_chromosomes(h5, key="normalized_counts"):

    n_cells = h5["cell_barcodes"][:].shape[0]
    all_chromosomes = list(h5[key].keys())

    number_chromosomes = sorted([re.findall('[0-9XY]+', x)[0] for x in all_chromosomes])
    all_chromosomes = sort_chromosomes(number_chromosomes)

    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(
            h5[key]["chr" + chr][:][0:n_cells, :]
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


def rename_fastq(s_name, sample_name=None):
    """
        renames the merged fastqs according to the bcl2fastq naming convention
        Sample input name: MERGED_BSSE_QGF_123456_ZXVN2SHG5_1_QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz
    """
    split_name = s_name.split("_")
    if sample_name is not None:
        new_name = "_".join(
            difflib.get_close_matches(sample_name, split_name, 1) + split_name[-4:]
        )
    else:
        new_name = "_".join(split_name[6:7] + split_name[-4:])

    return new_name


rename_fastq(
    "MERGED_BSSE_QGF_123456_ZXVN2SHG5_1__QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz",
    sample_name="QWEERTY_-TES2",
)


l = "MERGED_BSSE_QGF_123456_ZXVN2SHG5_1__QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz".split(
    "_"
)
sample_name = "QWEERTY-T"
import difflib

difflib.get_close_matches("QWEERTY-T", l, 1)

l[6:7]
l[-4:]
l
