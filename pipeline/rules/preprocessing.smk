import phenograph

rule identify_outliers:
    input:
        d_matrix_file = os.path.join(analysis_path, "filtering", analysis_prefix) + "__filtered_counts.csv",
    output:
        is_outlier_initial = os.path.join(analysis_path, "filtering", analysis_prefix) + "_is_outlier_initial.txt",
    run:
        data = np.loadtxt(input.d_matrix_file, delimiter=',')

        n_cells = data.shape[0]
        n_neighbours = max(int(n_cells / 10), 2) # avoid errors

        normalized_data = data / np.sum(data, axis=1).reshape(-1,1)
        communities, graph, Q = phenograph.cluster(data=normalized_data, k=n_neighbours, n_jobs=1, jaccard=True)

        is_outlier_initial = np.array([[True if c==-1 else False for c in communities]])
        is_outlier_initial = is_outlier_initial.reshape(n_cells, 1)

        np.savetxt(output.is_outlier_initial, is_outlier_initial)
