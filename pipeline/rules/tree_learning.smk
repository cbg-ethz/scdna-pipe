import json

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

rule learn_candidate_cluster_tree:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        n_reps = config["inference"]["cluster_trees"]["n_reps"],
        c_penalise = config["inference"]["c_penalise"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        posfix = "",
        scicone_path = scicone_path,
        output_temp_path = output_temp_path,
    input:
        segmented_counts = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_candidate_counts.csv",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_candidate_region_sizes.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_candidate_neutral_states.txt"
    output:
        clustering_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustering_score_candidate.txt",
        ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_candidate.txt",
        ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_candidate.csv",
        ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_candidate.tsv",
        ctree_robustness_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_robustness_candidate.txt",
        ctree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_candidate.json",
        ctree_cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_candidate.csv",
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        segmented_neutral_states = np.loadtxt(input.segmented_neutral_states, delimiter=',')
        alpha = 1./segmented_counts.shape[1]
        n_bins = np.sum(segmented_region_sizes)
        n_cells = segmented_counts.shape[0]

        sci = SCICoNE(params.scicone_path, params.output_temp_path, persistence=False, postfix=params.posfix)

        # Run cluster trees
        sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=params.n_reps, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters,
                                max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, copy_number_limit=params.copy_number_limit,
                                c_penalise=params.c_penalise, region_neutral_states=segmented_neutral_states)
        sci.best_cluster_tree.read_tree_str(sci.best_cluster_tree.tree_str, num_labels=True)

        # Store best cluster tree
        with open(output.ctree, "w") as file:
            for line in sci.best_cluster_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.ctree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')
        np.savetxt(output.ctree_cell_node_assignments, sci.best_cluster_tree.outputs['cell_node_ids'], delimiter='\t')

        with open(output.ctree_robustness_score, "w") as file:
            file.write(f"{sci.cluster_tree_robustness_score}\n")

        with open(output.clustering_score, "w") as file:
            file.write(str(sci.clustering_score))

        with open(output.ctree_json, 'w') as file:
            json.dump(sci.best_cluster_tree.node_dict, file, cls=NumpyEncoder)

        np.savetxt(output.ctree_cell_labels, np.array(sci.best_cluster_tree.cell_node_labels).astype(int), delimiter='\t')

rule learn_final_cluster_tree:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        n_reps = config["inference"]["cluster_trees"]["n_reps"],
        c_penalise = config["inference"]["c_penalise"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        diploid_maximum_frac = config["breakpoint_detection"]["diploid_maximum_frac"],
        posfix = "",
        scicone_path = scicone_path,
        output_temp_path = output_temp_path,
    input:
        segmented_counts = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_final_counts.csv",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_final_region_sizes.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_final_neutral_states.txt",
        is_diploid = os.path.join(analysis_path, "inferred_cnvs", analysis_prefix) + "_is_diploid_candidate.txt",
    output:
        clustering_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustering_score_final_diploid.txt",
        ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_diploid.txt",
        ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_final_diploid.csv",
        ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_final_diploid.tsv",
        ctree_robustness_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_robustness_final_diploid.txt",
        ctree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_diploid.json",
        ctree_cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_final_diploid.csv",
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        segmented_neutral_states = np.loadtxt(input.segmented_neutral_states, delimiter=',')
        alpha = 1./segmented_counts.shape[1]
        n_bins = np.sum(segmented_region_sizes)
        n_cells = segmented_counts.shape[0]

        is_diploid = np.loadtxt(input.is_diploid).astype(bool)
        run = False if np.count_nonzero(is_diploid) <= n_cells * params.diploid_maximum_frac else True # Don't re-run if didn't re-run breakpoint_detection

        if run:
            sci = SCICoNE(params.scicone_path, params.output_temp_path, persistence=False, postfix=params.posfix)

            # Run cluster trees
            sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=params.n_reps, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters,
                                    max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, copy_number_limit=params.copy_number_limit,
                                    c_penalise=params.c_penalise, region_neutral_states=segmented_neutral_states)
            sci.best_cluster_tree.read_tree_str(sci.best_cluster_tree.tree_str, num_labels=True)

            # Store best cluster tree
            with open(output.ctree, "w") as file:
                for line in sci.best_cluster_tree.tree_str.splitlines():
                    file.write(f"{line}\n")
            np.savetxt(output.ctree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')
            np.savetxt(output.ctree_cell_node_assignments, sci.best_cluster_tree.outputs['cell_node_ids'], delimiter='\t')

            with open(output.ctree_robustness_score, "w") as file:
                file.write(f"{sci.cluster_tree_robustness_score}\n")

            with open(output.clustering_score, "w") as file:
                file.write(str(sci.clustering_score))

            with open(output.ctree_json, 'w') as file:
                json.dump(sci.best_cluster_tree.node_dict, file, cls=NumpyEncoder)

            np.savetxt(output.ctree_cell_labels, np.array(sci.best_cluster_tree.cell_node_labels).astype(int), delimiter='\t')
        else:
            # Copy previous
            from shutil import copyfile
            clustering_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustering_score_candidate.txt"
            ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_candidate.txt"
            ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_candidate.csv"
            ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_candidate.tsv"
            ctree_robustness_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_robustness_candidate.txt"
            ctree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_candidate.json"
            ctree_cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_candidate.csv"
            copyfile(clustering_score, output.clustering_score)
            copyfile(ctree, output.ctree)
            copyfile(ctree_inferred_cnvs, output.ctree_inferred_cnvs)
            copyfile(ctree_cell_node_assignments, output.ctree_cell_node_assignments)
            copyfile(ctree_robustness_score, output.ctree_robustness_score)
            copyfile(ctree_json, output.ctree_json)
            copyfile(ctree_cell_labels, output.ctree_cell_labels)


rule learn_tetraploid_cluster_tree:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        n_reps = config["inference"]["cluster_trees"]["n_reps"],
        c_penalise = config["inference"]["c_penalise"],
        copy_number_limit = config["inference"]["copy_number_limit"] + 2,
        posfix = "",
        scicone_path = scicone_path,
        output_temp_path = output_temp_path
    input:
        ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_diploid.txt", # to avoid conflicts
        segmented_counts = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_final_counts.csv",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_final_region_sizes.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_final_neutral_states.txt",
    output:
        clustering_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustering_score_final_tetraploid.txt",
        ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_tetraploid.txt",
        ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_final_tetraploid.csv",
        ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_final_tetraploid.tsv",
        ctree_robustness_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_robustness_final_tetraploid.txt",
        ctree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_tetraploid.json",
        ctree_cell_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cell_labels_final_tetraploid.csv",
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        segmented_neutral_states = np.loadtxt(input.segmented_neutral_states, delimiter=',')

        segmented_neutral_states = segmented_neutral_states * 2

        alpha = 1./segmented_counts.shape[1]
        n_bins = np.sum(segmented_region_sizes)
        n_cells = segmented_counts.shape[0]
        sci = SCICoNE(params.scicone_path, params.output_temp_path, persistence=False, postfix=params.posfix + 'TETRA')

        # Run cluster trees
        sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=params.n_reps, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters,
                                max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, copy_number_limit=params.copy_number_limit,
                                c_penalise=params.c_penalise, region_neutral_states=segmented_neutral_states)
        sci.best_cluster_tree.read_tree_str(sci.best_cluster_tree.tree_str, num_labels=True)

        # Adjust to diploid root
        sci.best_cluster_tree.adjust_to_wgd(threshold=0.98)
        # sci.best_cluster_tree.node_dict['0']['label'] = '0'
        sci.best_cluster_tree.update_tree_str()

        # Store best cluster tree
        with open(output.ctree, "w") as file:
            for line in sci.best_cluster_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.ctree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')
        np.savetxt(output.ctree_cell_node_assignments, sci.best_cluster_tree.outputs['cell_node_ids'], delimiter='\t')

        with open(output.ctree_robustness_score, "w") as file:
            file.write(f"{sci.cluster_tree_robustness_score}\n")

        with open(output.clustering_score, "w") as file:
            file.write(str(sci.clustering_score))

        with open(output.ctree_json, 'w') as file:
            json.dump(sci.best_cluster_tree.node_dict, file, cls=NumpyEncoder)

        try:
            np.savetxt(output.ctree_cell_labels, np.array(sci.best_cluster_tree.cell_node_labels).astype(int), delimiter='\t')
        except ValueError:
            labs = [ord(char.lower()) - 97 for char in sci.best_cluster_tree.cell_node_labels]
            np.savetxt(output.ctree_cell_labels, np.array(labs).astype(int), delimiter='\t')

rule get_unique_cnvs:
    input:
        cluster_tree_json = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_final_{root}.json",
    output:
        unique_cnv_profiles = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__unique_cluster_tree_cnvs_final_{root}.csv",
        node_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_node_labels_final_{root}.csv",
        node_sizes = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_node_sizes_final_{root}.csv",
    run:
        with open(input.cluster_tree_json) as json_file:
            cluster_tree = json.load(json_file)

        unique_cnvs = []
        node_labels = []
        node_sizes = []
        for node in cluster_tree:
            if cluster_tree[node]['size'] > 0:
                node_labels.append(cluster_tree[node]['label'])
                unique_cnvs.append(np.array(cluster_tree[node]['cnv']))
                node_sizes.append(cluster_tree[node]['size'])
        order = np.argsort(np.array(node_labels))
        unique_cnvs = np.array(unique_cnvs)[order]
        node_labels = np.array(node_labels)[order]
        node_sizes = np.array(node_sizes)[order]

        print("saving the unique cnv profiles...")
        np.savetxt(
            output.unique_cnv_profiles,
            unique_cnvs.astype('int'),
            delimiter=",",
            fmt='%d'
        )

        print("saving the node labels...")
        np.savetxt(
            output.node_labels,
            node_labels.astype('int'),
            delimiter=" ",
            fmt='%d'
        )

        print("saving the node sizes...")
        np.savetxt(
            output.node_sizes,
            node_sizes.astype('int'),
            delimiter=",",
            fmt='%d'
        )
