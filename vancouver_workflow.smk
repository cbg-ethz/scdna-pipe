from scgenpy import *
from scgenpy.preprocessing.utils import *

import glob
import os
import h5py
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from collections import Counter
import re
import warnings
sns.set_style("ticks")

bin_size = config['bin_size']
fastqs_path = config['fastqs_path']
moved_fastqs_path = os.path.join(fastqs_path, "merged", "tricked")
analysis_path = config['analysis_path']

analysis_prefix = config['analysis_prefix']
seq_prefix = config["sequencing_prefix"]

sample_name = config['sample_name']
gene_lists_path = config['genes_path']
gene_coordinates_path = os.path.join(gene_lists_path, 'ensembl_hg19_annotations.tsv')
general_main_gene_list_path = os.path.join(gene_lists_path, 'general', f"{config['general_main_gene_list']}")

try:
    tree_rep = config["inference"]["cluster_trees"]["n_reps"]
except KeyError:
    tree_rep = 10

tree_outputs = ["cluster_tree", "full_tree"]

sa = SecondaryAnalysis(
    sample_name=analysis_prefix,
    output_path=analysis_path,
    h5_path=None
)

# import rules
include: os.path.join(workflow.basedir, "rules", "tree_learning.smk")
include: os.path.join(workflow.basedir, "rules", "breakpoint_detection.smk")
include: os.path.join(workflow.basedir, "rules", "process_cnvs.smk")
include: os.path.join(workflow.basedir, "rules", "plotting.smk")

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
             repeat_id=[x for x in range(0,tree_rep)]),

        cluster_tree_inferred_cnvs_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree_cnvs.csv",\
            repeat_id=[x for x in range(0,tree_rep)]),

        best_trees = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.txt", tree_name=tree_outputs),
        best_tree_inferred_cnvs = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}_cnvs.csv", tree_name=tree_outputs),

        unique_cnv_profiles = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__unique_cluster_tree_cnvs.csv",
        tree_cluster_sizes =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__tree_cluster_sizes.csv",

        tree_graphviz = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.graphviz", tree_name=tree_outputs),
        tree_figure =  expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.png", tree_name=tree_outputs),

        nu_on_cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree.txt",
        nu_on_cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree_cnvs.csv",

        full_tree_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__full_tree.txt",\
             repeat_id=[x for x in range(0,tree_rep)]),

        full_tree_inferred_cnvs_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__full_tree_cnvs.csv",\
            repeat_id=[x for x in range(0,tree_rep)]),

        robustness_results = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{tree_name}_robustness.txt", tree_name=tree_outputs),
        full_tree_inferred_cnvs_png = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__full_tree_cnvs.png"
    run:
        print("echo rule all")
