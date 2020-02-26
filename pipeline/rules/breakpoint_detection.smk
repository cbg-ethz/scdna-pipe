def filter_lr(lr_matrix, H=None):
    freq = np.fft.fftfreq(lr_matrix.shape[-1], 1) # 1 Hz sampling rate

    filtered_lr = np.empty(lr_matrix.shape)
    for c in range(lr_matrix.shape[0]):
        X = np.fft.fft(lr_matrix[c])
        Y = X * H
        y = np.fft.ifft(Y)
        filtered_lr[c] = y

    return filtered_lr

def set_region_neutral_states(all_region_stops, known_region_stops, known_region_neutral_states):
    """
        Creates a numpy array of size n_regions with the neutral state of each region.
        all_region_stops and known_region_stops must include the final bin
    """
    all_region_stops = np.array(all_region_stops)
    known_region_stops = np.array(known_region_stops)
    known_region_neutral_states = np.array(known_region_neutral_states)
    all_region_neutral_states = np.zeros(all_region_stops.shape)
    start_idx = 0
    for i, known_stop in enumerate(known_region_stops):
        end_idx = np.where(all_region_stops==known_stop)[0][0]
        all_region_neutral_states[start_idx:end_idx+1] = known_region_neutral_states[i]
        start_idx = end_idx

    return all_region_neutral_states

rule detect_breakpoints:
    params:
        window_size = config["breakpoint_detection"]["window_size"],
        verbosity = config["breakpoint_detection"]["verbosity"],
        threshold = config["breakpoint_detection"]["threshold"],
        bp_limit = config["breakpoint_detection"]["bp_limit"],
        bp_detection_path = os.path.join(analysis_path, "breakpoint_detection"),
        postfix = analysis_prefix,
        scicone_path = scicone_path,
        output_temp_path = output_temp_path
    input:
        d_matrix_file = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        matrix_shape = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts_shape.txt",
        chr_stops_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv",
        excluded_bins_path = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt"
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

        chr_stops = pd.read_csv(input.chr_stops_path, sep="\t", index_col=0)

        # Account for filtered out bins
        excluded_bins = np.loadtxt(input.excluded_bins_path)
        bin_chr_indicator = np.genfromtxt(input.bin_chr_indicator_path, dtype=None).astype(str)
        bin_chr_indicator_filtered = bin_chr_indicator[np.where(excluded_bins==0)[0]]
        chr_list = np.array(chr_stops['chr']).astype(str) # in order
        chr_stop_bins_filtered = []
        for chr in chr_list:
            chr_stop_bins_filtered.append(np.where(bin_chr_indicator_filtered==chr)[0][-1])
        chr_stop_bins_filtered = np.array(chr_stop_bins_filtered)[:-1]

        freq = np.fft.fftfreq(n_bins, 1)
        df = 0.015
        gpl = np.exp(- ((freq-1/(2*params.window_size))/(2*df))**2)  # pos. frequencies
        gmn = np.exp(- ((freq+1/(2*params.window_size))/(2*df))**2)
        g = gpl + gmn

        sci = SCICoNE(params.scicone_path, params.output_temp_path, persistence=True, postfix=params.postfix)

        data = np.loadtxt(input.d_matrix_file, delimiter=',')
        bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=params.threshold, bp_limit=params.bp_limit, compute_sp=False, evaluate_peaks=False)

        filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
        bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=params.threshold, bp_limit=params.bp_limit, lr=filtered_lr, input_breakpoints=chr_stop_bins_filtered)

        os.rename(params.postfix+"bp_segmented_regions.txt", output.segmented_regions)
        os.rename(params.postfix+"bp_segmented_region_sizes.txt", output.segmented_region_sizes)

rule segment_regions:
    input:
        filtered_counts = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
        segmented_regions = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_regions.txt",
        segmented_region_sizes = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_region_sizes.txt",
        chr_stops_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__chr_stops.tsv",
        excluded_bins_path = os.path.join(analysis_path, "filtering", analysis_prefix) + "__excluded_bins.csv",
        bin_chr_indicator_path = os.path.join(analysis_path, "genomic_coordinates", analysis_prefix) + "__bin_chr_indicator.txt"
    output:
        segmented_counts = os.path.join(analysis_path,\
                "breakpoint_detection", analysis_prefix) + "_segmented_counts.csv",
        segmented_counts_shape = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_counts_shape.txt",
        segmented_neutral_states = os.path.join(analysis_path, "breakpoint_detection", analysis_prefix) + "__segmented_neutral_states.txt"
    benchmark:
        "benchmark/segment_regions.tsv"
    run:
        print("loading the filtered counts...")
        filtered_counts = np.loadtxt(input.filtered_counts, delimiter=',')
        n_cells = filtered_counts.shape[0]
        n_bins = filtered_counts.shape[1]
        region_sizes = np.loadtxt(input.segmented_region_sizes)
        regions = np.loadtxt(input.segmented_regions)
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

        print("loading the chromosomal information...")
        chr_stops = pd.read_csv(input.chr_stops_path, sep="\t", index_col=0)
        chr_stop_bins = chr_stops.index
        chr_neutral_states = chr_stops['neutral_state']

        # Account for filtered out bins
        excluded_bins = np.loadtxt(input.excluded_bins_path)
        bin_chr_indicator = np.genfromtxt(input.bin_chr_indicator_path, dtype=None).astype(str)
        bin_chr_indicator_filtered = bin_chr_indicator[np.where(excluded_bins==0)[0]]
        chr_list = np.array(chr_stops['chr']).astype(str) # in order
        chr_stop_bins_filtered = []
        for chr in chr_list:
            chr_stop_bins_filtered.append(np.where(bin_chr_indicator_filtered==chr)[0][-1])
        chr_stop_bins_filtered = np.array(chr_stop_bins_filtered)
        # append final bin to regions
        regions = np.append(regions, chr_stop_bins_filtered[-1])
        segmented_neutral_states = set_region_neutral_states(regions, chr_stop_bins_filtered, chr_neutral_states)

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

        print("saving the region neutral states...")
        np.savetxt(
            output.segmented_neutral_states,
            segmented_neutral_states,
            delimiter=",",
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
