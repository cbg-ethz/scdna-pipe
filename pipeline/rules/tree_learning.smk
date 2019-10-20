rule clustering:
    input:
        normalised_regions = os.path.join(analysis_path, "normalisation", analysis_prefix) + "__normalised_regions.csv"
    output:
        clustering_score = os.path.join(analysis_path, "clustering", analysis_prefix) + "__clustering_score.txt",
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
        posfix = "cluster_trees" + "_{cluster_tree_rep}",
        alpha = config["inference"]["cluster_trees"]["alpha"]
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
        "benchmark/learn_cluster_trees_{cluster_tree_rep}.tsv"
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
                f"--move_probs={move_probs_str}", f"--seed={wildcards.cluster_tree_rep}", f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}",
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
            os.rename(f"{params.posfix}_cell_node_ids.tsv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_cell_node_ids.tsv")
            os.rename(f"{params.posfix}_cell_region_cnvs.csv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_cell_region_cnvs.csv")
            os.rename(f"{params.posfix}_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_markov_chain.csv")
            os.rename(f"{params.posfix}_rel_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_rel_markov_chain.csv")
            os.rename(f"{params.posfix}_acceptance_ratio.csv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_acceptance_ratio.csv")
            os.rename(f"{params.posfix}_gamma_values.csv", f"{debug_info_with_ap}__{wildcards.cluster_tree_rep}_gamma_values.csv")

rule cluster_tree_robustness:
    params:
        robustness_thr = config["inference"]["robustness_thr"],
        n_iters = config["inference"]["cluster_trees"]["n_iters"],
        n_reps = config["inference"]["cluster_trees"]["n_reps"],
        verbosity = config["inference"]["verbosity"]
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
            warnings.warn("The trees found are not robust, you may want to change the configurations")

        if params.verbosity > 0:
            debug_info_out_path = os.path.join(analysis_path, "tree_learning", "debug_info", "cluster_tree")
            debug_info_with_ap = os.path.join(debug_info_out_path, analysis_prefix)

            n_iters = params.n_iters
            all_chains = np.empty((n_iters,0))
            print("MCMC convergence plots...")
            for i in tqdm(range(params.n_reps)):
                cluster_chain = np.loadtxt(f"{debug_info_with_ap}__{i}_markov_chain.csv",delimiter=",")
                cluster_chain = cluster_chain.reshape(-1,1)
                all_chains = np.append(all_chains, cluster_chain, axis = 1)

            plt.figure(figsize=(20, 6))
            plt.plot(all_chains)
            plt.savefig(
                f"{debug_info_with_ap}__cluster_tree_mcmcs.png"
            )
            plt.close()

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
        "benchmark/learn_full_trees_{full_tree_rep}.tsv"
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
