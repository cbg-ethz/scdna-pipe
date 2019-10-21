import glob
import os
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from secondary_analysis import SecondaryAnalysis
from collections import Counter
import re
sns.set()


bin_size = config['bin_size']
analysis_path = config['analysis_path']

symlinks_path = os.path.join(analysis_path, '..', 'symlinks')
sym_raw_path = os.path.join(symlinks_path, "raw")
cellranger_path = os.path.join(analysis_path, "cellranger")
sample_name = config['sample_name']
cr_sample_name = sample_name[:-3] # e.g. MHELAVELA_S2 becomes MHELAVELA
all_genes_path = config['secondary_analysis']['all_genes_path']

analysis_prefix = config['analysis_prefix']

try:
    cluster_tree_rep = config["inference"]["cluster_trees"]["n_reps"]
except KeyError:
    cluster_tree_rep = 10

try:
    full_tree_rep = config["inference"]["full_trees"]["n_reps"]
except KeyError:
    full_tree_rep = 10

sa = SecondaryAnalysis(
    sample_name=analysis_prefix,
    output_path=analysis_path,
    h5_path=None,
    genes_path=None,
    all_genes_path=all_genes_path,
)

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
    '''
        renames the merged fastqs according to the bcl2fastq naming convention
        Sample input name: MERGED_BSSE_QGF_123456_ZXVN2SHG5_1_QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz
    '''
    split_name = s_name.split('_')
    new_name = '_'.join(split_name[6:7]+split_name[-4:])
    return new_name

def tree_to_graphviz(tree_path):
    """
        reads the file containing trees converts it to graphviz format
        :param tree_path: path to the tree file.
        :return: string object containing the graphviz formatted tree
    """
    with open(tree_path) as f:
        list_tree_file = list(f)

    graphviz_header = ["digraph { \n", "node [style=filled,color=\"#D4C0D6\"]"
                "edge [arrowhead=none, color=\"#602A86\"]"]
    
    graphviz_labels = []
    graphviz_links = []

    graphviz_labels.append("0[label=\"Neutral\"]") # root

    for line in list_tree_file: 
        if line.startswith("node 0:"):
            continue
        elif line.startswith("node"):
            comma_splits = line.split(",")
            
            comma_first = re.split(" |:",comma_splits[0])
            node_id = comma_first[1]
            p_id = comma_first[4]
            comma_rest = comma_splits[1:]
            comma_rest[0] = comma_rest[0].lstrip('[')
            comma_rest[-1] = comma_rest[-1].rstrip(']\n')
            merged_labels = [] 
            [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
            k_end = k_begin
            for term in comma_rest[1:]: # events vector
                [k,v] = (int(x) for x in term.split(":"))
                if k==k_end+1 and v==previous_v:
                    k_end = k # update the end
                else:
                    if k_begin == k_end:
                        merged_labels.append(f"{previous_v:+}R{k_begin}")
                    else:
                        merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
                    k_begin = k_end = k          
                previous_v = v
            # print the last one
            if k_begin == k_end:
                merged_labels.append(f"{previous_v:+}R{k_begin}")
            else:
                merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
            
            str_merged_labels = " ".join(f"{x}\n" if i%10 == 0 and i>0 else str(x) for i, x in enumerate(merged_labels))
            graphviz_labels.append(f"{node_id}[label=\"{str_merged_labels}\"]")
            graphviz_links.append(f"{p_id} -> {node_id}")

    return graphviz_header+graphviz_labels+graphviz_links+["}"]

# import rules
include: os.path.join(workflow.basedir, "rules", "breakpoint_detection.smk")
include: os.path.join(workflow.basedir, "rules", "tree_learning.smk")

onstart:
    print(f"Workflow main directory: {workflow.basedir}")

rule all:
    input:
        segmented_regions = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt",

        normalised_bins_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins.png",
        normalised_bins_clustered_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered.png",
        normalised_bins_clustered_bps_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered_bps.png",

        segmented_counts = os.path.join(analysis_path,\
        "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",

        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv",

        clustering_score = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clustering_score.txt",
        clusters_phenograph_assignment = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clusters_phenograph_assignment.tsv",

        avg_counts = os.path.join(analysis_path,\
                "clustering", analysis_prefix) + "_avg_counts.csv",

        empty_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__empty_tree.txt",
        empty_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__empty_tree_cnvs.csv",

        cluster_tree_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree.txt",\
             repeat_id=[x for x in range(0,cluster_tree_rep)]),

        cluster_tree_inferred_cnvs_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree_cnvs.csv",\
            repeat_id=[x for x in range(0,cluster_tree_rep)]),

        robustness_results = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_cluster_tree_robustness.txt",

        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cnvs.csv",

        unique_cnv_profiles = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__unique_cluster_tree_cnvs.csv",
        tree_cluster_sizes =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__tree_cluster_sizes.csv",

        cluster_tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.graphviz",
        cluster_tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.png",

        nu_on_cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree.txt",
        nu_on_cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree_cnvs.csv",

        full_tree_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__full_tree.txt",\
             repeat_id=[x for x in range(0,full_tree_rep)]),

        full_tree_inferred_cnvs_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__full_tree_cnvs.csv",\
            repeat_id=[x for x in range(0,full_tree_rep)])

    output:
        
    run:
        print("echo rule all")

rule visualise_trees:
    input:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt"
    output:
        cluster_tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.graphviz",
        cluster_tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.png"
    run:
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
