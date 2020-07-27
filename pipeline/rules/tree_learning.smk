rule learn_cluster_tree:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        n_reps = config["inference"]["cluster_trees"]["n_reps"],
        c_penalise = config["inference"]["c_penalise"],
        copy_number_limit = config["inference"]["copy_number_limit"],
        posfix = "",
        remove_cluster_outliers = lambda w: True if w.stage == 'candidate' else False,
        scicone_path = scicone_path,
        output_temp_path = output_temp_path
    input:
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts_{stage}.csv",
        segmented_region_sizes = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes_{stage}.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_neutral_states_{stage}.txt"
    output:
        clustering_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustering_score_{stage}.txt",
        ctree = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_{stage}.txt",
        ctree_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_{stage}.csv",
        ctree_cell_node_assignments = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_cell_node_ids_{stage}.tsv",
        ctree_robustness_score = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_robustness_{stage}.txt",
        is_outlier = os.path.join(analysis_path, "filtering", analysis_prefix) + "_is_outlier_{stage}.txt"
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
        sci.learn_tree(segmented_counts, segmented_region_sizes, remove_cluster_outliers=params.remove_cluster_outliers, n_reps=params.n_reps, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters,
                                max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, copy_number_limit=params.copy_number_limit,
                                c_penalise=params.c_penalise, region_neutral_states=segmented_neutral_states)

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

        is_outlier = np.isnan(sci.best_cluster_tree.outputs['cell_node_ids'].ravel()).reshape(n_cells,1).astype(bool)
        np.savetxt(outputs.is_outlier, is_outlier)

rule cell_assignment:
    input:
        inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__cluster_tree_inferred_cnvs_{stage}.csv",
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions_{stage}.csv",
        normalised_bins = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_bins_{stage}.csv"
    output:
        unique_cnv_profiles = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__unique_cluster_tree_cnvs_{stage}.csv",
        tree_node_sizes =  os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__tree_node_sizes_cluster_tree_{stage}.txt",
        clustered_inferred_cnvs = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustered_cluster_tree_inferred_cnvs_{stage}.csv",
        clustered_normalised_regions = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustered_normalised_regions_{stage}.csv",
        clustered_normalised_bins = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustered_normalised_bins_{stage}.csv",
        clustered_labels = os.path.join(analysis_path, "tree_learning", analysis_prefix) + "__clustered_labels_{stage}.csv"
    benchmark:
        "benchmarks/cell_assignments_{stage}.tsv"
    run:
        inferred_cnvs = np.loadtxt(input.inferred_cnvs, delimiter=',') # cell-wise profiles
        normalised_bins = np.loadtxt(input.normalised_bins, delimiter=',')
        normalised_regions = np.loadtxt(input.normalised_regions, delimiter=',')

        unique_cnvs, tree_node_sizes = np.unique(inferred_cnvs, axis=0, return_counts=True) # clone-wise profiles
        if len(unique_cnvs.shape) == 1: # if only one cluster
            unique_cnvs = unique_cnvs.reshape(1, -1)

        # Sort clones by distance to diploid profile
        dist_to_diploid = []
        diploid_profile = np.ones([unique_cnvs.shape[1]]) * 2
        for c_id in range(unique_cnvs.shape[0]):
            dist_to_diploid.append(np.linalg.norm(unique_cnvs[c_id]-diploid_profile))
        order = np.argsort(dist_to_diploid)
        unique_cnvs = unique_cnvs[order]
        tree_node_sizes = tree_node_sizes[order]

        labels = np.empty(inferred_cnvs.shape[0])
        for c_id in range(unique_cnvs.shape[0]):
            cells = np.where(np.all(inferred_cnvs==unique_cnvs[c_id], axis=1))[0]
            labels[cells] = c_id

        # Sort cells by sorted clone index
        cell_order = np.argsort(labels)
        clustered_labels = labels[cell_order]
        inferred_cnvs = inferred_cnvs[cell_order]
        normalised_bins = normalised_bins[cell_order]
        normalised_regions = normalised_regions[cell_order]

        print("saving the unique cnv profiles...")
        np.savetxt(
            output.unique_cnv_profiles,
            unique_cnvs,
            delimiter=",",
            fmt='%d'
        )
        print("saving the tree cluster sizes...")
        np.savetxt(
            output.tree_node_sizes,
            tree_node_sizes,
            delimiter=",",
            fmt='%d'
        )
        print("saving the sorted cnv profiles...")
        np.savetxt(
            output.clustered_inferred_cnvs,
            inferred_cnvs,
            delimiter=",",
            fmt='%d'
        )
        print("saving the sorted normalised bins...")
        np.savetxt(
            output.clustered_normalised_bins,
            normalised_bins,
            delimiter=",",
            fmt='%d'
        )
        print("saving the sorted normalised regions...")
        np.savetxt(
            output.clustered_normalised_regions,
            normalised_regions,
            delimiter=",",
            fmt='%d'
        )
        print("saving the sorted cell labels...")
        np.savetxt(
            output.clustered_labels,
            labels,
            delimiter=",",
            fmt='%d'
        )
