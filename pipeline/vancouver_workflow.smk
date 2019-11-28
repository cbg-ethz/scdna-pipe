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
from secondary_analysis.utils import *
from collections import Counter
import re
import warnings
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
    tree_rep = config["inference"]["cluster_trees"]["n_reps"]
except KeyError:
    tree_rep = 10

tree_outputs = ["cluster_tree", "full_tree"]

sa = SecondaryAnalysis(
    sample_name=analysis_prefix,
    output_path=analysis_path,
    h5_path=None,
    genes_path=None,
    all_genes_path=all_genes_path,
)


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

    output:
        
    run:
        print("echo rule all")

rule visualise_trees:
    input:
        tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.txt"
    output:
        tree_graphviz = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.graphviz",
        tree_figure =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__{tree_name}.png"
    run:
        tree_as_list = tree_to_graphviz(input.tree)
        
        for line in tree_as_list:
            print(f"{line}\n")
        
        with open(output.tree_graphviz, "w") as file:
            for line in tree_as_list:
                file.write(f"{line}\n")

        try:
            cmd_output = subprocess.run(["dot", "-Tpng", f"{output.tree_graphviz}", "-o", f"{output.tree_figure}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

rule plot_inferred_cnvs:
    input:
        full_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__full_tree_cnvs.csv"
    output:
        full_tree_inferred_cnvs_png = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__full_tree_cnvs.png"
    benchmark:
        "benchmark/plot_inferred_cnvs.tsv"
    run:
        cnvs = np.loadtxt(input.full_tree_inferred_cnvs, delimiter=',', dtype=int)
        cmap = sns.color_palette("RdBu_r", 5)
        print("plotting")
        plt.figure(figsize=(24, 8))
        ax = sns.heatmap(cnvs, cmap=cmap)
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
        colorbar.set_ticklabels(['0', '1', '2', '3', '4'])
        plt.savefig(output.full_tree_inferred_cnvs_png)
        plt.close()