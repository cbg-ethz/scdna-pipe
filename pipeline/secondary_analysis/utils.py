import numpy as np
import h5py
import re


def merge_chromosomes(h5, key="normalized_counts"):

    n_cells = h5["cell_barcodes"][:].shape[0]
    all_chromosomes = list(h5[key].keys())
    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(
            h5[key][chr][:][0:n_cells, :]
        )  # select only the cells, not cell groups

    cell_all_chrs = np.concatenate(cnv_matrices, axis=1)
    return cell_all_chrs


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


def get_gene_cn_df(genes, cnvs_arr, bin_size, chr_stops):
    """
        Creates and returns the dataframe of copy numbers, genes by cluster ids
        :param genes: The input list of genes to be specified
        :param cnvs_arr: Copy number matrix per cell per bin
        :param bin_size: Integer constant speciying size of each bin
        :param chr_stops: Pandas DataFrame
            Contains the final bins of each chromosome (the df index)
        :return: CN dataframe of genes by clusters
    """
    print(f"genes shape: {genes.shape}")
    print(f"cnvs_arr shape: {cnvs_arr.shape}")
    if len(cnvs_arr.shape) == 1:
        cnvs_arr = cnvs_arr.reshape(1, -1)
    cluster_ids = range(cnvs_arr.shape[0])
    gene_cn_df = pd.DataFrame(index=cluster_ids)

    # for each gene
    for index, row in tqdm(genes.iterrows(), total=genes.shape[0]):
        start_bin = int(row["Gene start (bp)"] / bin_size)
        stop_bin = int(row["Gene end (bp)"] / bin_size)
        chromosome = str(row["Chromosome/scaffold name"])
        gene_name = row["Gene name"]

        if chromosome != "1":  # coordinates are given by chromosome
            chr_start = (
                chr_stops.iloc[
                    np.where(chr_stops.index == chromosome)[0][0] - 1
                ].values[0]
                + 1
            )
            start_bin = start_bin + chr_start
            stop_bin = stop_bin + chr_start

        gene_cn_per_cluster = []
        for c_id in cluster_ids:
            cn_states = cnvs_arr[c_id, start_bin : stop_bin + 1]
            median_cn_cell_bin = np.nanmedian(
                cn_states
            )  # all the bins within the gene, min due to biology
            gene_cn_per_cluster.append(median_cn_cell_bin)
        gene_cn_df[gene_name] = gene_cn_per_cluster

    print("Transposing the dataframe...")
    gene_cn_df = gene_cn_df.T
    print("Sorting the genes...")
    gene_cn_df.sort_index(inplace=True)

    print(gene_cn_df.head())
    return gene_cn_df


def get_bin_gene_region_df(
    bin_size, genes, chr_stops, region_stops, bin_is_excluded, important_genes=None
):
    """
        Creates a DataFrame with the gene and region corresponding to each bin
        :param bin_size: integer
        :param chr_stops: df indicating the final (unfiltered) bin of each chromosome
        :param region_stops: df indicating the final (filtered) bin of each chromosome
        :param bin_is_excluded: df indicating unmappable bins
        (optional) :param important_genes: list of important genes to be flagged
        :return: DataFrame of (gene, region, original_bin)
    """

    bin_gene_region_df = pd.DataFrame(index=range(chr_stops.iloc[-1].values[0]))
    bin_gene_region_df["gene"] = None
    bin_gene_region_df["chr"] = None
    bin_gene_region_df["region"] = None
    if important_genes is None:
        bin_gene_region_df["is_important"] = True
    else:
        bin_gene_region_df["is_important"] = False

    # for each gene
    for index, row in tqdm(genes.iterrows(), total=genes.shape[0]):
        start_bin = int(row["Gene start (bp)"] / bin_size)
        stop_bin = int(row["Gene end (bp)"] / bin_size)
        chromosome = str(row["Chromosome/scaffold name"])

        if chromosome != "1":  # coordinates are given by chromosome
            chr_start = (
                chr_stops.iloc[
                    np.where(chr_stops.index == chromosome)[0][0] - 1
                ].values[0]
                + 1
            )
            start_bin = start_bin + chr_start
            stop_bin = stop_bin + chr_start

        # Check if bins already contain an important gene: if they do, skip
        if (important_genes is None) or (
            (important_genes is not None)
            and np.all(
                bin_gene_region_df.loc[start_bin : stop_bin + 1, "is_important"]
                == False
            )
        ):
            gene_name = row["Gene name"]
            bin_gene_region_df.loc[start_bin : stop_bin + 1, "gene"] = gene_name
            bin_gene_region_df.loc[start_bin : stop_bin + 1, "chr"] = chromosome
            if important_genes is not None:
                if gene_name in important_genes:
                    bin_gene_region_df.loc[
                        start_bin : stop_bin + 1, "is_important"
                    ] = True
                else:
                    # Need this to deal with overlapping gene coordinates
                    bin_gene_region_df.loc[
                        start_bin : stop_bin + 1, "is_important"
                    ] = False

    # Now we have a (bin, gene) dataframe. Let's remove the unmappable bins
    bin_gene_region_df = bin_gene_region_df.iloc[np.where(bin_is_excluded == 0)[0]]

    # Reset the index to indicate new bin indices after filtering
    bin_gene_region_df = bin_gene_region_df.reset_index(drop=False)
    bin_gene_region_df = bin_gene_region_df.rename(columns={"index": "original_bin"})

    # Get the regions
    start_bin = 0
    for index, row in tqdm(region_stops.iterrows(), total=region_stops.shape[0]):
        stop_bin = row.values[0]
        bin_gene_region_df.loc[
            start_bin : stop_bin + 1, "region"
        ] = index  # regions are 0 indexed
        start_bin = row.values[0]

    return bin_gene_region_df


def get_genes_in_region(region, bin_gene_region_df, important_only=False):
    """
        Returns a list of genes present in the given region index
        :param region: integer
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param important_only: boolean
            indicating if only important genes should be returned
        :return: list of gene names in region
    """
    gene_list = bin_gene_region_df["gene"][
        np.where(bin_gene_region_df["region"] == region)[0]
    ].unique()

    gene_list = gene_list[gene_list != None].tolist()

    if important_only:
        # Subset only the important genes
        important_genes = bin_gene_region_df["gene"][
            bin_gene_region_df["is_important"]
        ].unique()
        gene_list = [str(gene) for gene in gene_list if gene in important_genes]

    return gene_list


def get_region_with_gene(gene, bin_gene_region_df):
    """
        Returns the region index containing a gene
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: region index (integer)
    """
    gene = bin_gene_region_df["region"][
        np.where(bin_gene_region_df["gene"] == gene)[0]
    ].values[0]

    return gene


def convert_event_region_to_gene(
    region_event_str,
    bin_gene_region_df,
    important_only=False,
    genes_to_highlight=None,
    highlight_color="red",
):
    """
        Returns a string indicating gene-wise events in affected region
        Examples:
                "+2R174"     -> ["+2BRAF", "+2MALAT1"]
                "-1R656:658" -> ["-1JAK2", "-1MLNA", "-1CDK4"]
        :param region_event_str: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param important_only: boolean
            indicating if only important genes should be returned
        :param genes_to_highlight: list
            genes that should be displayed in a different color
        :param highlight_color: str
            color to use in genes to highlight
        :return: list of str
    """
    # Get event (-2, -1, +1, +2, etc)
    event_str = region_event_str[:2]
    region_str = region_event_str[3:]
    if ":" in region_str:  # multiple regions: "-1R656:658"
        aux = [int(region) for region in region_str.split(":")]
        region_list = np.arange(aux[0], aux[1] + 1)
    else:
        region_list = [int(region_str)]

    gene_list = []
    for region in region_list:
        genes_in_region = get_genes_in_region(
            region, bin_gene_region_df, important_only=important_only
        )
        gene_list.append(genes_in_region)

    gene_list = [item for sublist in gene_list for item in sublist]

    # Highlight some genes
    if genes_to_highlight is not None:
        for index, gene in enumerate(gene_list):
            if gene in genes_to_highlight:
                gene_list[index] = (
                    "<font color="
                    + "'"
                    + highlight_color
                    + "'"
                    + ">"
                    + gene
                    + "</font>"
                )

    gene_string = "[" + ",".join(gene_list) + "]"
    if len(gene_list) == 0:
        gene_event_str = ""
    else:
        gene_event_str = event_str + gene_string

    return gene_event_str


def convert_node_regions_to_genes(
    node_str,
    bin_gene_region_df,
    important_only=False,
    genes_to_highlight=None,
    **kwargs,
):
    """
        Returns a string indicating gene events and total number of
        amplifications and deletions in node
        Examples:
                "+2R174 -1R656:658" -> "+2[BRAF,MALAT1] -1[JAK2,MLNA,CDK4]\n(1+, 1-)"
        :param node_str: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param important_only: boolean
            indicating if only important genes should be returned
        :param genes_to_highlight: list
            genes that should be displayed in a different color
        :return: str
    """
    region_event_strs = node_str.split(" ")
    num_events = len(region_event_strs)

    num_amplifications = 0
    gene_event_str = []

    for region_event_str in region_event_strs:
        s = convert_event_region_to_gene(
            region_event_str,
            bin_gene_region_df,
            important_only=important_only,
            genes_to_highlight=genes_to_highlight,
            highlight_color=kwargs["highlight_color"],
        )
        gene_event_str.append(s)

        if region_event_str[0] == "+":
            num_amplifications += 1

    gene_event_str = [x for x in gene_event_str if x]
    # Add newline after each x events
    gene_event_str = " ".join(
        f"{x}<br/>" if i % 2 == 0 and i > 0 else str(x)
        for i, x in enumerate(gene_event_str)
    )
    # gene_event_str = ' '.join(gene_event_str)

    num_deletions = num_events - num_amplifications
    num_events, num_amplifications, num_deletions
    num_events_str = "{} +, {} -".format(num_amplifications, num_deletions)
    num_events_str

    node_str = gene_event_str + " <br/>" + "(" + num_events_str + ")"
    if gene_event_str == "":
        node_str = num_events_str

    return node_str


def tree_to_graphviz(
    tree_path,
    gene_labels=False,
    bin_gene_region_df=None,
    genes_to_highlight=None,
    **kwargs,
):
    """
        reads the file containing trees converts it to graphviz format
        :param tree_path: path to the tree file.
        :param gene_labels: whether to label nodes with genes
        :param bin_gene_region_df: Pandas DataFrame with gene-region correspondence
        :param genes_to_highlight: List containing genes to highlight
        :return: string object containing the graphviz formatted tree
    """
    with open(tree_path) as f:
        list_tree_file = list(f)

    graphviz_header = [
        "digraph { \n",
        'node [style=filled,color="#D4C0D6"]' 'edge [arrowhead=none, color="#602A86"]',
    ]

    graphviz_labels = []
    graphviz_links = []

    graphviz_labels.append('0[label="Neutral"]')  # root

    for line in list_tree_file:
        if line.startswith("node 0:"):
            continue
        elif line.startswith("node"):
            comma_splits = line.split(",")

            comma_first = re.split(" |:", comma_splits[0])
            node_id = comma_first[1]
            p_id = comma_first[4]
            comma_rest = comma_splits[1:]
            comma_rest[0] = comma_rest[0].lstrip("[")
            comma_rest[-1] = comma_rest[-1].rstrip("]\n")
            merged_labels = []
            [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
            k_end = k_begin
            for term in comma_rest[1:]:  # events vector
                [k, v] = (int(x) for x in term.split(":"))
                if k == k_end + 1 and v == previous_v:
                    k_end = k  # update the end
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

            str_merged_labels = " ".join(
                f"{x}<br/>" if i % 10 == 0 and i > 0 else str(x)
                for i, x in enumerate(merged_labels)
            )
            if gene_labels and bin_gene_region_df is not None:
                node_str = " ".join(merged_labels)  # "+1R75 +1R218:219 +1R221:223"
                str_merged_labels = convert_node_regions_to_genes(
                    node_str,
                    bin_gene_region_df,
                    important_only=True,
                    genes_to_highlight=genes_to_highlight,
                    **kwargs,
                )

            graphviz_labels.append(
                f"{node_id}[label=<{str_merged_labels}>]"
            )  # use < > to allow HTML
            graphviz_links.append(f"{p_id} -> {node_id}")

    return graphviz_header + graphviz_labels + graphviz_links + ["}"]


def rename_fastq(s_name):
    """
        renames the merged fastqs according to the bcl2fastq naming convention
        Sample input name: MERGED_BSSE_QGF_123456_ZXVN2SHG5_1_QWEERTY_T_scD_250c_r1v1_0_SI-GA-H5_S1_L003_I1_001.fastq.gz
    """
    split_name = s_name.split("_")
    new_name = "_".join(split_name[6:7] + split_name[-4:])
    return new_name
