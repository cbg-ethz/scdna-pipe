rule detect_breakpoints:
    params:
        binary = config["breakpoint_detection"]["bin"],
        window_size = config["breakpoint_detection"]["window_size"],
        verbosity = config["breakpoint_detection"]["verbosity"],
        threshold = config["breakpoint_detection"]["threshold"],
        bp_limit = config["breakpoint_detection"]["bp_limit"],
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
                f"--verbosity={params.verbosity}", f"--threshold={params.threshold}", f"--bp_limit={params.bp_limit}"])
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
