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
        phenograph_distance = os.path.join(analysis_path, "clustering", analysis_prefix) + "__phenograph_distance.csv",
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

rule detect_breakpoints:
    params:
        binary = config["breakpoint_detection"]["bin"],
        window_size = config["breakpoint_detection"]["window_size"],
        verbosity = config["breakpoint_detection"]["verbosity"],
        threshold = config["breakpoint_detection"]["threshold"],
        bp_detection_path = os.path.join(analysis_path, "breakpoint_detection"),
        posfix = analysis_prefix 
    input:
        d_matrix_file = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        matrix_shape = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts_shape.txt"
    output:
        segmented_regions = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        segmented_region_sizes = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt"
    run:
        try:
            os.makedirs(params.bp_detection_path)
        except FileExistsError:
            print("breakpoint detection directory already exists.")
        input_shape = np.loadtxt(input.matrix_shape)
        (n_cells, n_bins) = [int(input_shape[i]) for i in range(2)]

        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.d_matrix_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={params.window_size}", f"--postfix={params.posfix}",\
                f"--verbosity={params.verbosity}", f"--threshold={params.threshold}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")
        
        os.rename(params.posfix+"_segmented_regions.txt", output.segmented_regions)
        os.rename(params.posfix+"_segmented_region_sizes.txt", output.segmented_region_sizes)

rule plot_breakpoints:
    input:
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv",
        segmented_regions = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt"
    output:
        normalised_bins_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins.png",
        normalised_bins_clustered_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered.png",
        normalised_bins_clustered_bps_heatmap = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_normalised_bins_clustered_bps.png"
    benchmark:
        "benchmark/plot_breakpoints.tsv"
    run:
        from scipy.cluster.hierarchy import ward, leaves_list
        from scipy.spatial.distance import pdist

        print("loading the normalised bins...")
        normalised_bins = np.loadtxt(input.normalised_bins, delimiter=',')
        print(f"shape of normalised bins: {normalised_bins.shape}")

        print("loading the normalised regions...")
        normalised_regions = np.loadtxt(input.normalised_regions, delimiter=',')
        print(f"shape of normalised regions: {normalised_regions.shape}")

        bps = np.loadtxt(input.segmented_regions)
        print(f"number of breakpoints: {len(bps)}")

        mean_normalised_bins = normalised_bins.mean()
        mean_normalised_regions = normalised_regions.mean()

        Z = ward(pdist(normalised_regions))
        hclust_index = leaves_list(Z)
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        print("plotting the heatmap ...")
        plt.figure(figsize=(24, 8))
        ax = sns.heatmap(normalised_bins, vmax=2, cmap=cmap)
        plt.savefig(output.normalised_bins_heatmap)
        plt.close()

        print("plotting the clustered heatmap ...")
        plt.figure(figsize=(24, 8))
        ax = sns.heatmap(normalised_bins[hclust_index], vmax=2, cmap=cmap)
        plt.savefig(output.normalised_bins_clustered_heatmap)
        plt.close()

        print("plotting the clustered heatmap with breakpoints ...")
        plt.figure(figsize=(24, 8))
        ax = sns.heatmap(normalised_bins[hclust_index], vmax=2, cmap=cmap)
        ax = ax.vlines(bps, *ax.get_xlim(), colors='b', linestyles='dashed')
        plt.savefig(output.normalised_bins_clustered_bps_heatmap)
        plt.close()


rule segment_regions:
    input:
        filtered_counts = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        segmented_region_sizes = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt"
    output:
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt"
    benchmark:
        "benchmark/segment_regions.tsv"
    run:
        print("loading the filtered counts...")
        filtered_counts = np.loadtxt(input.filtered_counts, delimiter=',')
        n_cells = filtered_counts.shape[0]
        region_sizes = np.loadtxt(input.segmented_region_sizes)
        n_regions = len(region_sizes)
        sum_region_sizes = np.sum(region_sizes)
        condensed_mat = np.zeros((n_cells, n_regions))

        print("segmenting the bins...")
        for i in tqdm(range(n_cells)):
            region_id = 0
            region_count = 0
            # import ipdb; ipdb.set_trace() # debugging starts here
            for j in range(int(sum_region_sizes)):
                to_add = filtered_counts[i][j]
                condensed_mat[i][region_id] += to_add
                region_count += 1
                if region_count == region_sizes[region_id]:
                    region_id += 1
                    region_count = 0

        if not np.allclose(condensed_mat.sum(axis=1), filtered_counts.sum(axis=1)):
            raise AssertionError(
                "not all values of the sums before & after "
                "segmentation are close")

        print("saving the segmented regions...")
        np.savetxt(
            output.segmented_counts,
            condensed_mat,
            delimiter=",",
        )

        np.savetxt(
            output.segmented_counts_shape,
            condensed_mat.shape
        )


rule normalise_counts:
    input:
        filtered_counts = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv"
    output:
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv"
    benchmark:
        "benchmark/normalise_counts.tsv"    
    run:
        from sklearn.preprocessing import normalize

        print("loading the filtered counts...")
        filtered_counts = np.loadtxt(input.filtered_counts, delimiter=',')
        print("normalising the bins...")
        normalized_filtered_bins = normalize(filtered_counts, axis=1, norm="l1")
        # normalise but make the row sum equal to n_bins
        normalized_filtered_bins *= normalized_filtered_bins.shape[1]
        print(f"shape of normalised bins: {normalized_filtered_bins.shape}")
        print("saving the normalised bins...")
        np.savetxt(
            output.normalised_bins,
            normalized_filtered_bins,
            delimiter=",",
        )
        
        print("loading the segmented counts...")
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        print("normalising the regions...")
        normalized_filtered_regions = normalize(segmented_counts, axis=1, norm="l1")
        print(f"shape of normalised regions: {normalized_filtered_regions.shape}")
        print("saving the normalised bins...")
        np.savetxt(
            output.normalised_regions,
            normalized_filtered_regions,
            delimiter=",",
        )

rule clustering:
    input:
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv"
    output:
        clustering_score = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clustering_score.txt",
        phenograph_distance = os.path.join(analysis_path, "clustering", analysis_prefix) + "__phenograph_distance.csv",
        clusters_phenograph_assignment = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clusters_phenograph_assignment.tsv"
    benchmark:
        "benchmark/clustering.tsv"
    run:
        sa.apply_phenograph(input.normalised_regions)


rule create_averaged_region_matrix:
    """
    Creates the averaged regions by using cluster assignments
    """
    input:
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        clusters_phenograph_assignment = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clusters_phenograph_assignment.tsv"
    output:
        avg_counts = os.path.join(analysis_path,\
                "clustering", analysis_prefix) + "_avg_counts.csv"
    benchmark:
        "benchmark/create_averaged_region_matrix.tsv"
    run:
        print("loading the segmented counts...")
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        print("normalising the regions...")

        phenograph_assignments = pd.read_csv(input.clusters_phenograph_assignment, sep='\t')
        communities = phenograph_assignments.cluster.values

        community_dict = dict((Counter(communities)))
        community_ids = sorted(list(community_dict))

        cluster_sizes = [v for (k,v) in sorted(community_dict.items())]
        print(f"cluster sizes: {cluster_sizes}")

        cells_by_cluster = []
        for cluster in community_ids:
            cells_by_cluster.append(segmented_counts[communities == cluster])

        avg_clusters = [m.mean(0) for m in cells_by_cluster]
        avg_clusters_df = pd.DataFrame(avg_clusters)
        print(f"shape of average clusters: {avg_clusters_df.shape}")

        replicated_df = pd.DataFrame(np.repeat(avg_clusters_df.values,cluster_sizes,axis=0))

        np.savetxt(output.avg_counts, replicated_df.values, delimiter=",")

rule learn_empty_tree:
    params:
        binary = config["inference"]["bin"],
        ploidy = config["inference"]["ploidy"],
        verbosity = config["inference"]["verbosity"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        n_iters = config["inference"]["learn_nu"]["n_iters"],
        n_nodes = config["inference"]["learn_nu"]["n_nodes"],
        move_probs = config["inference"]["learn_nu"]["move_probs"],
        seed = config["inference"]["seed"],
        posfix = "empty_tree"
    input:
        avg_counts = os.path.join(analysis_path,\
                "clustering", analysis_prefix) + "_avg_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",
        segmented_region_sizes = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt"
    output:
        empty_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__empty_tree.txt",
        empty_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__empty_tree_cnvs.csv"
    benchmark:
        "benchmark/learn_empty_tree.tsv"
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]
        
        move_probs_str = ",".join(str(p) for p in params.move_probs)

        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.avg_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--move_probs={move_probs_str}", f"--seed={params.seed}", f"--region_sizes_file={input.segmented_region_sizes}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.empty_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.empty_tree_inferred_cnvs)

rule learn_cluster_trees:
    params:
        binary = config["inference"]["bin"],
        ploidy = config["inference"]["ploidy"],
        verbosity = config["inference"]["verbosity"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        n_iters = config["inference"]["cluster_trees"]["n_iters"],
        n_nodes = config["inference"]["cluster_trees"]["n_nodes"],
        move_probs = config["inference"]["cluster_trees"]["move_probs"],
        posfix = "cluster_trees" + "_{cluster_tree_rep}"
    input:
        avg_counts = os.path.join(analysis_path,\
                "clustering", analysis_prefix) + "_avg_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",
        segmented_region_sizes = os.path.join(analysis_path,\
             "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt",
        empty_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__empty_tree.txt"
    output:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{cluster_tree_rep}" + "__cluster_tree.txt",
        cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{cluster_tree_rep}" + "__cluster_tree_cnvs.csv"
    benchmark:
        "benchmark/learn_cluster_trees.tsv"
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]
        
        move_probs_str = ",".join(str(p) for p in params.move_probs)

        with open(input.empty_tree) as file:
            for l in file:
                l_parts = l.split(':')
                if l_parts[0] == 'Nu':
                    nu = l_parts[1].strip()
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.avg_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--move_probs={move_probs_str}", f"--seed={wildcards.cluster_tree_rep}", f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.cluster_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.cluster_tree_inferred_cnvs)

rule cluster_tree_robustness:
    params:
        robustness_thr = config["inference"]["robustness_thr"]
    input:
        cluster_tree_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree.txt",\
             repeat_id=[x for x in range(0,cluster_tree_rep)])
    output:
        robustness_results = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_cluster_tree_robustness.txt"
    benchmark:
        "benchmark/cluster_tree_robustness.tsv"
    run:
        tree_scores = get_tree_scores(input.cluster_tree_with_rep)
        max_score = max(tree_scores)
        score_dispersions = [abs(max_score-x) for x in tree_scores]

        is_robust = [x<=1 for x in score_dispersions]
        robustness_ratio = sum(is_robust) / len(is_robust)
        print(f"Robustness ratio: {robustness_ratio}")
        
        with open(output.robustness_results, "w") as f:
            f.write(f"Scores differences from the max: {score_dispersions}\n")
            f.write(f"Robustness vector of trees: {is_robust}\n")
            f.write(f"Robustness ratio: {robustness_ratio}\n")    

        if robustness_ratio < params.robustness_thr:
            raise Exception("The trees found are not robust, you may want to change the configurations") 

rule pick_max_cluster_tree:
    input:
        cluster_tree_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree.txt",\
             repeat_id=[x for x in range(0,cluster_tree_rep)]),
        cluster_tree_inferred_cnvs_with_rep = expand(os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{repeat_id}" + "__cluster_tree_cnvs.csv",\
            repeat_id=[x for x in range(0,cluster_tree_rep)])
    output:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cnvs.csv"
    benchmark:
        "benchmark/pick_max_cluster_tree.tsv"
    run:
        import operator

        trees_sorted = sorted(input.cluster_tree_with_rep)
        trees_inferred_cnvs_sorted = sorted(input.cluster_tree_inferred_cnvs_with_rep)

        tree_scores = get_tree_scores(trees_sorted)
        max_index, max_score = max(enumerate(tree_scores), key=operator.itemgetter(1))

        os.symlink(trees_sorted[max_index], output.cluster_tree)
        os.symlink(trees_inferred_cnvs_sorted[max_index], output.cluster_tree_inferred_cnvs)

rule learn_nu_on_cluster_tree:
    params:
        binary = config["inference"]["bin"],
        ploidy = config["inference"]["ploidy"],
        verbosity = config["inference"]["verbosity"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        n_iters = config["inference"]["learn_nu_cluster_trees"]["n_iters"],
        move_probs = config["inference"]["learn_nu_cluster_trees"]["move_probs"],
        n_nodes = 0, # needed only for the output naming
        seed = config["inference"]["seed"],
        posfix = "nu_on_cluster_tree"
    input:
        cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree.txt",
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt"
    output:
        nu_on_cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree.txt",
        nu_on_cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree_cnvs.csv"
    benchmark:
        "benchmark/learn_nu_on_cluster_tree.tsv"
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]
        
        move_probs_str = ",".join(str(p) for p in params.move_probs)

        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.segmented_counts}", f"--tree_file={input.cluster_tree}",\
             f"--n_regions={n_regions}",f"--n_nodes={params.n_nodes}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}",\
                f"--move_probs={move_probs_str}", f"--seed={params.seed}", f"--region_sizes_file={input.segmented_region_sizes}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.nu_on_cluster_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.nu_on_cluster_tree_inferred_cnvs)

rule learn_full_trees:
    params:
        binary = config["inference"]["bin"],
        ploidy = config["inference"]["ploidy"],
        verbosity = config["inference"]["verbosity"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        n_iters = config["inference"]["full_trees"]["n_iters"],
        n_nodes = config["inference"]["full_trees"]["n_nodes"],
        move_probs = config["inference"]["full_trees"]["move_probs"],
        posfix = "full_trees" + "_{full_tree_rep}"
    input:
        nu_on_cluster_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__nu_on_cluster_tree.txt",
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt"
    output:
        full_tree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{full_tree_rep}" + "__full_tree.txt",
        full_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "_{full_tree_rep}" + "__full_tree_cnvs.csv"
    benchmark:
        "benchmark/learn_full_trees.tsv"
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]
        
        move_probs_str = ",".join(str(p) for p in params.move_probs)

        with open(input.nu_on_cluster_tree) as file:
            for l in file:
                l_parts = l.split(':')
                if l_parts[0] == 'Nu':
                    nu = l_parts[1].strip()
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.segmented_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--tree_file={input.nu_on_cluster_tree}",\
                f"--move_probs={move_probs_str}", f"--seed={wildcards.full_tree_rep}", f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.full_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.full_tree_inferred_cnvs)

rule cell_assignment:
    input:
        cluster_tree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cnvs.csv"
    output:
        unique_cnv_profiles = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__unique_cluster_tree_cnvs.csv",
        tree_cluster_sizes =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__tree_cluster_sizes.csv"
    benchmark:
        "benchmark/cell_assignments.tsv"
    run:
        inferred_cnvs = np.loadtxt(input.cluster_tree_inferred_cnvs, delimiter=',')
        unique_cnvs, tree_cluster_sizes = np.unique(inferred_cnvs, axis=0, return_counts=True)

        print("saving the unique cnv profiles...")
        np.savetxt(
            output.unique_cnv_profiles,
            unique_cnvs,
            delimiter=",",
            fmt='%d'
        )
        print("saving the tree cluster sizes...")
        np.savetxt(
            output.tree_cluster_sizes,
            tree_cluster_sizes,
            delimiter=",",
            fmt='%d'
        )

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
