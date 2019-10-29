import os
import itertools
import pandas as pd
import numpy as np
from collections import Counter
from sklearn.metrics import roc_curve, auc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from tqdm import tqdm as tqdm
import re
import phenograph
sns.set(rc={'figure.figsize':(15.7,8.27)})

'''
parameters
'''

'''
Trees of size n=10, 20 and 40
n_regions = n, 2n and 4n
10,000 bins and 500 cells
20,000, 40,000 and 80,000 reads per cell
'''

n_nodes = config["simulate"]["n_nodes"] # values: [10,20,30]
n_regions = [n_nodes,2*n_nodes,4*n_nodes]
n_bins = config["simulate"]["n_bins"]
n_reads = config["simulate"]["n_reads"]
nu = config["simulate"]["nu"]


try:
    all_n_tps = config["simulate"]["n_reps"]
except KeyError:
    all_n_tps = 100 

try:
    n_inference_reps = config["inference"]["full_trees"]["n_reps"]
except KeyError:
    n_inference_reps = 10

try:
    tree_rep = config["inference"]["cluster_trees"]["n_reps"]
except KeyError:
    tree_rep = 10

n_cells = config["simulate"]["n_cells"]
n_iters = config["inference"]["full_trees"]["n_iters"]  # int(1000000*n_nodes/10)

sim_output_file_exts = ['d_mat.txt','ground_truth.txt','region_sizes.txt', 'tree.txt'] #, 'inferred_cnvs.txt', 'tree_inferred.txt', 'HMMcopy_inferred.txt','inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']
bp_output_file_exts = ['segmented_regions.txt', 'segmented_region_sizes.txt', 'sp_vec.csv', 'log_posterior_vec.csv', 'expected_k_vec.csv', 'all_bps_comparison.csv', 'aic_vec.csv']
# trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt'] # , 'inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

SIM_OUTPUT= config["simulate"]["output"]
BP_OUTPUT = config["breakpoint_detection"]["output"]
PHENO_OUTPUT = config["phenograph"]["output"]
TREES_OUTPUT = config["inference"]["output"]
sim_prefix=config["simulate"]["prefix"]


rule all:
    input:
        cluster_tree_inferred_cnvs = expand(f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_cluster_tree_cnvs.csv'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)], tree_rep_id=[x for x in range(0,tree_rep)]),

        breakpoints = expand(os.path.join(f"{BP_OUTPUT}_{sim_prefix}_trees", str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads', '{rep_id}' +'_' + '{bp_output_ext}')\
        ,bp_output_ext=bp_output_file_exts, regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)]),

        hmmcopy_inferred_cnvs = expand(SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_HMMcopy_inferred.txt'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)]),

        phenograph_assignments = expand(f'{PHENO_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)])

        
    run:
        print("rule all")

rule run_sim:
    params:
        sim_bin = config["simulate"]["bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        all_n_tps = all_n_tps,
        n_iters = n_iters,
        nu = nu
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ground_truth.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        tree = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree.txt'
    shell:
        "{params.sim_bin} --n_regions {wildcards.regions} --n_reads {wildcards.reads} --n_iters {params.n_iters} --n_cells {params.n_cells} --n_bins {params.n_bins} --n_nodes \
        {params.n_nodes} --verbosity 0 --ploidy 2 --postfix {wildcards.rep_id} --nu {params.nu} --seed {wildcards.rep_id}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.csv {output.d_mat}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_ground_truth.csv {output.ground_truth}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_region_sizes.txt {output.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_tree.txt {output.tree}"

rule detect_breakpoints:
    params:
        binary = config["breakpoint_detection"]["bin"],
        window_size = config["breakpoint_detection"]["window_size"],
        verbosity = config["breakpoint_detection"]["verbosity"],
        threshold = config["breakpoint_detection"]["threshold"],
        bp_limit = config["breakpoint_detection"]["bp_limit"],
        posfix = sim_prefix
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        segmented_regions = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_regions.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt',
        
        sp_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_sp_vec.csv',
        log_posterior_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_log_posterior_vec.csv',
        expected_k_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_expected_k_vec.csv',
        all_bps_comparison = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_all_bps_comparison.csv',
        aic_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_aic_vec.csv'
    run:
        
        # n_cells and n_bins are global
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.d_matrix_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={params.window_size}", f"--postfix={n_nodes}nodes_{wildcards.regions}regions_"
                f"{wildcards.reads}reads_{wildcards.rep_id}",\
                f"--verbosity={params.verbosity}", f"--threshold={params.threshold}", f"--bp_limit={params.bp_limit}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        for filename in bp_output_file_exts:
            os.rename(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{filename}", os.path.join(f"{BP_OUTPUT}_{sim_prefix}_trees", str(n_nodes) + 'nodes_' + wildcards.regions +'regions_'+ wildcards.reads + 'reads', wildcards.rep_id +'_' + filename))

rule segment_regions:
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' +\
         '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt"
    run:
        filtered_counts = np.loadtxt(input.d_matrix_file, delimiter=',')
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
        np.savetxt(output.segmented_counts, condensed_mat, delimiter=",")
        np.savetxt(output.segmented_counts_shape, condensed_mat.shape)

rule phenograph_clustering_trees:
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv"
    output:
        clusters_phenograph_assignment = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'
    threads:
        config["phenograph"]["threads"]
    run:
        normalised_regions = np.loadtxt(input.segmented_counts, delimiter=",")
        n_cells = normalised_regions.shape[0]
        print(f"n_cells: {str(n_cells)}")
        n_neighbours = int(n_cells / 10)
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        communities, graph, Q = phenograph.cluster(data=normalised_regions, k=n_neighbours, n_jobs=threads, jaccard=True)

        print(f"Communities: {communities}")
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]

        communities_df.to_csv(output.clusters_phenograph_assignment, sep="\t", index=False)

rule create_averaged_region_matrix:
    """
    Creates the averaged regions by using cluster assignments
    """
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        clusters_phenograph_assignment = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'
    output:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv'
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
        posfix = f"empty_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}rep"
    input:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv',
        segmented_counts_shape = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) +\
         'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        empty_tree = f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_empty_tree.txt',
        empty_tree_inferred_cnvs = f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}__empty_tree_cnvs.csv'
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

        if params.verbosity > 0:
            debug_info_out_path = os.path.join(analysis_path, "tree_learning", "debug_info", "empty_tree")
            debug_info_with_ap = os.path.join(debug_info_out_path, analysis_prefix)
            if not os.path.exists(debug_info_out_path):
                os.makedirs(debug_info_out_path)
            os.rename(f"{params.posfix}_cell_node_ids.tsv", f"{debug_info_with_ap}__cell_node_ids.tsv")
            os.rename(f"{params.posfix}_cell_region_cnvs.csv", f"{debug_info_with_ap}__cell_region_cnvs.csv")
            os.rename(f"{params.posfix}_markov_chain.csv", f"{debug_info_with_ap}__markov_chain.csv")
            os.rename(f"{params.posfix}_rel_markov_chain.csv", f"{debug_info_with_ap}__rel_markov_chain.csv")
            os.rename(f"{params.posfix}_acceptance_ratio.csv", f"{debug_info_with_ap}__acceptance_ratio.csv")
            os.rename(f"{params.posfix}_gamma_values.csv", f"{debug_info_with_ap}__gamma_values.csv")

rule learn_cluster_trees:
    params:
        binary = config["inference"]["bin"],
        ploidy = config["inference"]["ploidy"],
        verbosity = config["inference"]["verbosity"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        n_iters = config["inference"]["cluster_trees"]["n_iters"],
        n_nodes = config["inference"]["cluster_trees"]["n_nodes"],
        move_probs = config["inference"]["cluster_trees"]["move_probs"],
        posfix = f"cluster_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}rep_{tree_rep_id}",
        alpha = config["inference"]["cluster_trees"]["alpha"]
    input:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv',
        segmented_counts_shape = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) +\
         'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt',
        empty_tree = f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_empty_tree.txt'
    output:
        cluster_tree = f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_cluster_tree_cnvs.csv'
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
                f"--move_probs={move_probs_str}", f"--seed={wildcards.tree_rep}", f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}",
                f"--alpha={params.alpha}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.cluster_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.cluster_tree_inferred_cnvs)

        if params.verbosity > 0:
            debug_info_out_path = os.path.join(analysis_path, "tree_learning", "debug_info", "cluster_tree")
            debug_info_with_ap = os.path.join(debug_info_out_path, analysis_prefix)
            if not os.path.exists(debug_info_out_path):
                os.makedirs(debug_info_out_path)
            os.rename(f"{params.posfix}_cell_node_ids.tsv", f"{debug_info_with_ap}__{wildcards.tree_rep}_cell_node_ids.tsv")
            os.rename(f"{params.posfix}_cell_region_cnvs.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_cell_region_cnvs.csv")
            os.rename(f"{params.posfix}_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_markov_chain.csv")
            os.rename(f"{params.posfix}_rel_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_rel_markov_chain.csv")
            os.rename(f"{params.posfix}_acceptance_ratio.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_acceptance_ratio.csv")
            os.rename(f"{params.posfix}_gamma_values.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_gamma_values.csv")


rule plot_breakpoints_thresholds:
    params:
        n_bins = n_bins,
        n_nodes = n_nodes
    input:
        gt_files = expand(os.path.join(f"{SIM_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{{regions}}'+'regions_'+ '{{reads}}'+'reads', '{rep_id}' +'_' + 'ground_truth.txt')\
        , rep_id=[x for x in range(0,all_n_tps)]),
        bps_files = expand(os.path.join(f"{BP_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{{regions}}'+'regions_'+ '{{reads}}'+'reads', '{rep_id}' +'_' + 'all_bps_comparison.csv')\
        , rep_id=[x for x in range(0,all_n_tps)])
    output:
        mean_tprs = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'mean_tpr_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        mean_fprs = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'mean_fpr_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        sum_tps = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'sum_tps_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        sum_fps = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'sum_fps_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png')
    run:
        all_bins = range(0,params.n_bins)
        all_tpr, all_fpr, all_n_tps, all_n_fps, all_bps_tables = ([] for i in range(5))

        threshold_coeffs = np.linspace(1.0,16.0, 50) # 50 thresholds btw 1 and 16

        for gt_file, bps_file in tqdm(zip(input.gt_files, input.bps_files)):
            bps = pd.read_csv(bps_file, header=None)
            bps.columns = ['idx','log_sp','stdev']
            bps['ranking'] = bps['log_sp'] / bps['stdev']
            # bps = bps.sort_values('ranking',ascending=False)
            bps = bps.dropna()
            
            all_bps_tables.append(bps)
            
            # get the ground truth
            cell_genotypes = pd.read_csv(gt_file, sep=',' ,header=None)
            cell_genotypes = cell_genotypes[cell_genotypes.columns[:-1]] # remove the last (only NaN) column
            cell_bps = cell_genotypes.diff(periods=1, axis=1) # apply diff to detect breakpoints
            cell_bps = cell_bps.fillna(value=0.0) # diff makes the 1st row NaN, make it zero
            cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
            grouped_cell_bps = cell_bps.sum(axis=0) # count the non-zeroes
            ground_truth = grouped_cell_bps[grouped_cell_bps > 0] # if they occur in at least 1 cell
            ground_truth = ground_truth.index.tolist()
            # end of ground truth    

            # correcting for the bps 1-2 bins nearby
            for index, row in bps.iterrows():
                idx_val = bps.loc[index, 'idx']
                for gt in ground_truth:
                    if (abs(idx_val - gt) <=2 and idx_val != gt):
                        print('correcting ' + str(idx_val) + '->' + str(gt))
                        bps.loc[index,'idx'] = gt
            
            tpr_values, fpr_values, n_tps, n_fps = [], [], [], []
            for thr in threshold_coeffs:
                predicted_positives = []
                predicted_negatives = []
                for index, row in bps.iterrows():
                    if row['ranking'] > thr:
                        predicted_positives.append(row['idx'])
                    else:
                        break 
                        
                #import ipdb; ipdb.set_trace()
                predicted_negatives = [i for i in all_bins if i not in predicted_positives]

                true_positives = [i for i in predicted_positives if i in ground_truth]
                false_positives = [i for i in predicted_positives if i not in ground_truth]

                true_negatives = [i for i in predicted_negatives if i not in ground_truth]
                false_negatives = [i for i in predicted_negatives if i in ground_truth]

                # import ipdb; ipdb.set_trace()
                assert(len(ground_truth) == (len(true_positives) + len(false_negatives)))
                tpr = len(true_positives) / len(ground_truth) # len(ground_truth)
                fpr = len(false_positives) / (params.n_bins - len(ground_truth)) # (len(false_positives) + len(true_negatives))
                tpr_values.append(tpr)
                fpr_values.append(fpr)
                n_tps.append(len(true_positives))
                n_fps.append(len(false_positives))
                
            
            all_tpr.append(tpr_values)
            all_fpr.append(fpr_values)
            all_n_tps.append(n_tps)
            all_n_fps.append(n_fps)


        # Average TRP, FPR, number of predicted positives for each threshold value over all runs
        print("Plotting...")
        tpr_df = pd.DataFrame(all_tpr)
        tpr_df.columns = threshold_coeffs
        mean_tpr = tpr_df.mean()
        ax = sns.lineplot(x=mean_tpr.index,y=mean_tpr.values).set_title('Average TPR values')
        plt.savefig(output.mean_tprs)
        plt.clf()

        fpr_df = pd.DataFrame(all_fpr)
        fpr_df.columns = threshold_coeffs
        mean_fpr = fpr_df.mean()
        ax = sns.lineplot(x=mean_fpr.index,y=mean_fpr.values).set_title('Average FPR values')
        plt.savefig(output.mean_fprs)
        plt.clf()

        n_tps_df = pd.DataFrame(all_n_tps)
        n_tps_df.columns = threshold_coeffs
        sum_tps = n_tps_df.sum()
        ax = sns.lineplot(x=sum_tps.index,y=sum_tps.values).set_title('Sum TP values')
        plt.savefig(output.sum_tps)
        plt.clf()

        n_fps_df = pd.DataFrame(all_n_fps)
        n_fps_df.columns = threshold_coeffs
        sum_fps = n_fps_df.sum()
        ax = sns.lineplot(x=sum_fps.index,y=sum_fps.values).set_title('Sum FP values')
        plt.savefig(output.sum_fps)
        plt.clf()

rule hmm_copy_inference:
    params:
        script = config["hmm_copy"]["script"],
        n_nodes = n_nodes,
        scratch = config["hmm_copy"]["scratch"],
        mem = config["hmm_copy"]["mem"],
        time = config["hmm_copy"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_HMMcopy_inferred.txt
        hmmcopy_inferred_cnvs = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_HMMcopy_inferred.txt'
    conda:
        config["hmm_copy"]["conda_env"]
    shell:
        " Rscript {params.script} {input.d_mat}" 

rule phenograph_clustering:
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        clusters_phenograph_assignment = f'{PHENO_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'
    threads:
        config["phenograph"]["threads"]
    run:
        normalised_regions = np.loadtxt(input.d_mat, delimiter=",")
        n_cells = normalised_regions.shape[0]
        print(f"n_cells: {str(n_cells)}")
        n_neighbours = int(n_cells / 10)
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        communities, graph, Q = phenograph.cluster(data=normalised_regions, k=n_neighbours, n_jobs=threads, jaccard=True)

        print(f"Communities: {communities}")
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]

        communities_df.to_csv(output.clusters_phenograph_assignment, sep="\t", index=False)