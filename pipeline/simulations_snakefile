import os
import itertools

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
    n_repetitions = config["simulate"]["n_reps"]
except KeyError:
    n_repetitions = 100 

try:
    n_inference_reps = config["inference"]["full_trees"]["n_reps"]
except KeyError:
    n_inference_reps = 10

n_cells = config["simulate"]["n_cells"]
n_iters = config["inference"]["full_trees"]["n_iters"]  # int(1000000*n_nodes/10)

sim_output_file_exts = ['d_mat.txt','ground_truth.txt','region_sizes.txt', 'tree.txt'] #, 'inferred_cnvs.txt', 'tree_inferred.txt', 'HMMcopy_inferred.txt','inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']
bp_output_file_exts = ['segmented_regions.txt', 'segmented_region_sizes.txt', 'sp_vec.csv', 'log_posterior_vec.csv', 'expected_k_vec.csv', 'all_bps_comparison.csv', 'aic_vec.csv']
# trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt'] # , 'inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

SIM_OUTPUT= config["simulate"]["output"]
BP_OUTPUT = config["breakpoint_detection"]["output"]
sim_prefix=config["simulate"]["prefix"]


rule all:
    input:
        breakpoints = expand(os.path.join(f"{BP_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads', '{rep_id}' +'_' + '{bp_output_ext}')\
        ,bp_output_ext=bp_output_file_exts, regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,n_repetitions)])
        # read_region_sims = expand(os.path.join(f"{SIM_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads', '{rep_id}' +'_' + '{sim_output_ext}')\
        # ,sim_output_ext=sim_output_file_exts, regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,n_repetitions)]),
        
    run:
        print("rule all")

rule run_sim:
    params:
        sim_bin = config["simulate"]["bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        n_repetitions = n_repetitions,
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
        segmented_regions = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_regions.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt',
        
        sp_vec = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_sp_vec.csv',
        log_posterior_vec = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_log_posterior_vec.csv',
        expected_k_vec = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_expected_k_vec.csv',
        all_bps_comparison = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_all_bps_comparison.csv',
        aic_vec = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_aic_vec.csv'
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
            os.rename(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{filename}", os.path.join(f"{BP_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + wildcards.regions +'regions_'+ wildcards.reads + 'reads', wildcards.rep_id +'_' + filename))






        # bp_outputs = expand(os.path.join(f"{BP_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads', '{rep_id}' +'_' + '{bp_output_files}')\
        # ,bp_output_files=bp_output_file_exts)
        # segmented_regions = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_regions.txt',
        # segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'