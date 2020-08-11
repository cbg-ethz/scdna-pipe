import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist
from tqdm import tqdm


def cluster_clones(inferred_cnvs, normalised_bins, normalised_regions):
    # label cells according to their clone
    labels = np.empty(inferred_cnvs.shape[0])
    unique_cnvs, cluster_sizes, idx = np.unique(
        inferred_cnvs, return_counts=True, return_index=True, axis=0
    )

    # Sort clones by distance to diploid profile
    dist_to_diploid = []
    diploid_profile = np.ones([unique_cnvs.shape[1]]) * 2
    for c_id in range(unique_cnvs.shape[0]):
        dist_to_diploid.append(np.linalg.norm(unique_cnvs[c_id] - diploid_profile))
    order = np.argsort(dist_to_diploid)
    unique_cnvs = unique_cnvs[order]
    cluster_sizes = cluster_sizes[order]

    for c_id in range(unique_cnvs.shape[0]):
        cells = np.where(np.all(inferred_cnvs == unique_cnvs[c_id], axis=1))[0]
        if np.any(np.isnan(unique_cnvs[c_id])):
            labels[cells] = "-"
        labels[cells] = c_id

    # Sort cells by sorted clone index
    cell_order = np.argsort(labels)
    labels = labels[cell_order]
    inferred_cnvs = inferred_cnvs[cell_order]
    normalised_bins = normalised_bins[cell_order]
    normalised_regions = normalised_regions[cell_order]

    # Within each clone, sort the cells's normalised regions via hierarchical clustering
    for c_id in range(unique_cnvs.shape[0]):
        if len(np.where(labels == c_id)[0]) > 1:
            Z = ward(pdist(normalised_regions[np.where(labels == c_id)[0]]))
            hclust_index = leaves_list(Z)
            normalised_bins[np.where(labels == c_id)[0]] = normalised_bins[
                np.where(labels == c_id)[0]
            ][hclust_index]

    return inferred_cnvs, normalised_bins, labels


def get_bin_gene_region_df(
    bin_size,
    gene_coordinates,
    chr_stops,
    region_stops,
    bin_is_excluded,
    cnvs=None,
    priority_genes=None,
):
    """
        Creates a DataFrame with the gene and region corresponding to each bin
        :param bin_size: integer
        :param gene_coordinates: coordinates of all genes
        :param chr_stops: df indicating the final (unfiltered) bin of each chromosome
        :param region_stops: df indicating the final (filtered) bin of each chromosome
        :param bin_is_excluded: df indicating unmappable bins
        (optional) :param cnvs: np.array of cnv of clone per bin
        (optional) :param priority_genes: list of priority genes to be flagged
        :return: DataFrame of (gene, chr, region, is_priority, filtered_bin, cnv_{})
    """

    bin_gene_region_df = pd.DataFrame(index=range(chr_stops.index[-1] + 1))

    bin_gene_region_df["region"] = None
    bin_gene_region_df["gene"] = [list() for _ in range(bin_gene_region_df.shape[0])]
    bin_gene_region_df["gene_id"] = [list() for _ in range(bin_gene_region_df.shape[0])]
    bin_gene_region_df["chr"] = [list() for _ in range(bin_gene_region_df.shape[0])]
    bin_gene_region_df["is_priority"] = [
        list() for _ in range(bin_gene_region_df.shape[0])
    ]
    bin_gene_region_df["biotype"] = [list() for _ in range(bin_gene_region_df.shape[0])]

    # for each gene
    for index, row in tqdm(
        gene_coordinates.iterrows(), total=gene_coordinates.shape[0]
    ):
        start_bin = int(row["start"] / bin_size)
        stop_bin = int(row["end"] / bin_size)
        chromosome = str(row["contig"])
        biotype = str(row["biotype"])

        if chromosome != "1":  # coordinates are given by chromosome
            chr_start = (
                chr_stops.iloc[np.where(chr_stops["chr"] == chromosome)[0][0] - 1].name
                + 1
            )
            start_bin = start_bin + chr_start
            stop_bin = stop_bin + chr_start

        gene = row["gene_name"]
        gene_id = row["gene_id"]
        for bin in range(start_bin, stop_bin + 1):
            bin_gene_region_df.loc[bin, "gene"].append(gene)
            bin_gene_region_df.loc[bin, "gene_id"].append(gene_id)
            bin_gene_region_df.loc[bin, "chr"].append(chromosome)
            bin_gene_region_df.loc[bin, "biotype"].append(biotype)

        if priority_genes is not None:
            if gene in priority_genes:
                for bin in range(start_bin, stop_bin + 1):
                    bin_gene_region_df.loc[bin, "is_priority"].append(True)
                print(bin_gene_region_df.loc[start_bin : stop_bin + 1, "gene"])
            else:
                # Need this to deal with overlapping gene coordinates
                for bin in range(start_bin, stop_bin + 1):
                    bin_gene_region_df.loc[bin, "is_priority"].append(False)
        else:
            for bin in range(start_bin, stop_bin + 1):
                bin_gene_region_df.loc[bin, "is_priority"].append(True)

    # Turn columns of lists into columns of strings with comma-separated values
    bin_gene_region_df["gene"] = [
        ",".join(map(str, l)) for l in bin_gene_region_df["gene"]
    ]
    bin_gene_region_df["gene_id"] = [
        ",".join(map(str, l)) for l in bin_gene_region_df["gene_id"]
    ]
    bin_gene_region_df["chr"] = [
        ",".join(map(str, l)) for l in bin_gene_region_df["chr"]
    ]
    bin_gene_region_df["is_priority"] = [
        ",".join(map(str, l)) for l in bin_gene_region_df["is_priority"]
    ]
    bin_gene_region_df["biotype"] = [
        ",".join(map(str, l)) for l in bin_gene_region_df["biotype"]
    ]

    # Indicate original_bin-filtered_bin correspondence
    bin_gene_region_df["filtered_bin"] = None
    bin_gene_region_df["filtered_bin"].iloc[
        np.where(np.array(bin_is_excluded) == False)[0]
    ] = np.arange(np.count_nonzero(np.array(bin_is_excluded) == False))

    # Get the regions
    start_bin = 0
    for index, row in tqdm(region_stops.iterrows(), total=region_stops.shape[0]):
        stop_bin = row.values[0]
        original_start_bin = np.where(bin_gene_region_df["filtered_bin"] == start_bin)[
            0
        ][0]
        original_stop_bin = np.where(bin_gene_region_df["filtered_bin"] == stop_bin)[0][
            0
        ]
        bin_gene_region_df.loc[
            original_start_bin : original_stop_bin + 1, "region"
        ] = index  # regions are 0 indexed
        start_bin = row.values[0]

    if cnvs is not None:
        # Add the CN values
        if len(cnvs.shape) == 1:
            cnvs = cnvs.reshape(1, -1)

        for c_id in range(cnvs.shape[0]):
            bin_gene_region_df["cnv_{}".format(c_id)] = cnvs[c_id]

    # Make sure excluded bins have no info
    bin_gene_region_df.loc[np.where(bin_is_excluded)[0], "region"] = None

    return bin_gene_region_df


def get_genes_in_region(region, bin_gene_region_df, priority_only=False):
    """
        Returns a list of genes present in the given region index
        :param region: integer
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param priority_only: boolean
            indicating if only priority genes should be returned
        :return: list of gene names in region
    """
    gene_lists = bin_gene_region_df["gene"][
        np.where(bin_gene_region_df["region"] == region)[0]
    ].values.tolist()
    gene_lists = ["" if x is np.nan else x for x in gene_lists]
    gene_lists = [sublist.split(",") for sublist in gene_lists]
    gene_list = [item for sublist in gene_lists for item in sublist]
    gene_list = [x for x in gene_list if x]

    if priority_only:
        is_priority_lists = bin_gene_region_df["is_priority"][
            np.where(bin_gene_region_df["region"] == region)[0]
        ].values.tolist()
        is_priority_lists = ["" if x is np.nan else x for x in is_priority_lists]
        is_priority_lists = [sublist.split(",") for sublist in is_priority_lists]
        is_priority_list = [item for sublist in is_priority_lists for item in sublist]
        is_priority_list = [x for x in is_priority_list if x]

        # Subset only the priority genes
        priority_gene_list = []
        for i, gene in enumerate(gene_list):
            if is_priority_list[i] == "True":
                priority_gene_list.append(gene)
        gene_list = priority_gene_list

    # Remove duplicates
    gene_list = np.unique(gene_list).tolist()

    return gene_list


def get_region_with_gene(gene, bin_gene_region_df, all=False):
    """
        Returns the region index containing a gene
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: region index (integer)
    """
    bin_gene_region_df = bin_gene_region_df.copy(deep=True)
    bin_gene_region_df["gene"] = (
        bin_gene_region_df["gene"]
        .astype(str)
        .apply(lambda x: x.split(","))
        .apply(lambda x: set(x))
    )
    bins_with_gene = bin_gene_region_df[
        bin_gene_region_df["gene"].apply(lambda x: gene in x)
    ].index.values

    region = np.array(bin_gene_region_df.loc[bins_with_gene, "region"].values)
    region = np.array([x for x in region if x is not None])

    if len(region[~np.isnan(region)]) > 0:
        # Get first non-None value
        region = region[~np.isnan(region)][0] if not all else region
        region = region.astype(int)
    else:
        region = np.nan

    return region


def get_surrounding_regions(gene, bin_gene_region_df):
    """
        Finds the regions to the left and to the right of a gene's bins
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: Tuple containing (left_region, right_region)
    """
    bin_gene_region_df = bin_gene_region_df.copy(deep=True)
    bin_gene_region_df["gene"] = (
        bin_gene_region_df["gene"]
        .astype(str)
        .apply(lambda x: x.split(","))
        .apply(lambda x: set(x))
    )
    bins = bin_gene_region_df[
        bin_gene_region_df["gene"].apply(lambda x: gene in x)
    ].index.values

    left_bin = bins[0]
    right_bin = bins[-1]

    left_regions = bin_gene_region_df["region"].values[:left_bin]
    left_regions = np.array([x for x in left_regions if x is not None])
    left_region = left_regions[~np.isnan(left_regions)][-1]

    right_regions = bin_gene_region_df["region"].values[right_bin + 1 :]
    right_regions = np.array([x for x in right_regions if x is not None])
    right_region = right_regions[~np.isnan(right_regions)][0]

    return (int(left_region), int(right_region))


def get_gene_cn_df(gene_list, bin_gene_region_df, impute=False):
    """
        Creates and returns the dataframe of copy numbers, genes by cluster ids
        :param gene_list: the input list of genes to be specified
        :param bin_gene_region_df: pd.DataFrame with cnv per bin, region and gene
        :param impute: replace NaN values with median of left and rightmost regions
        :return: CN dataframe of genes by subclone
    """
    n_subclones = np.count_nonzero(
        ["cnv" in column for column in bin_gene_region_df.columns]
    )
    cluster_ids = range(n_subclones)
    gene_cn_df = pd.DataFrame(index=cluster_ids)

    df = bin_gene_region_df.copy(deep=True)
    df["gene"] = (
        df["gene"].astype(str).apply(lambda x: x.split(",")).apply(lambda x: set(x))
    )

    gene_list = list(dict.fromkeys(gene_list))  # remove duplicates

    is_imputed = np.empty(len(gene_list))

    # for each gene
    i = -1
    for gene in tqdm(gene_list):
        i += 1
        gene_cn_per_cluster = []
        is_imputed[i] = False
        gene_exists = True

        for c_id in cluster_ids:
            bins = df[df["gene"].apply(lambda x: gene in x)].index.values
            median_cn = np.nanmedian(
                bin_gene_region_df["cnv_{}".format(c_id)][bins].values
            )
            # If NaN, impute with median value of regions surrounding it
            if np.isnan(median_cn) and impute:
                if len(bins) > 0:
                    # Get regions surrounding gene
                    left_region, right_region = get_surrounding_regions(
                        gene, bin_gene_region_df
                    )

                    # get CNV values of region surronding gene
                    left_cn = bin_gene_region_df["cnv_{}".format(c_id)][
                        bin_gene_region_df["region"] == left_region
                    ].iloc[0]
                    right_cn = bin_gene_region_df["cnv_{}".format(c_id)][
                        bin_gene_region_df["region"] == right_region
                    ].iloc[0]

                    median_cn = np.nanmedian([left_cn, right_cn])

                    is_imputed[i] = True
                else:
                    gene_exists = False
                    is_imputed[i] = np.nan
                    print(f"Gene {gene} not found.")

            if not np.isnan(median_cn):
                if median_cn > 2:
                    median_cn = int(np.floor(median_cn))
                elif median_cn < 2:
                    median_cn = int(np.ceil(median_cn))

            gene_cn_per_cluster.append(median_cn)

        if gene_exists:
            gene_cn_df[gene] = gene_cn_per_cluster

    print("Transposing the dataframe...")
    gene_cn_df = gene_cn_df.T
    if impute:
        gene_cn_df["is_imputed"] = is_imputed[~np.isnan(is_imputed)].tolist()
    # gene_cn_df = gene_cn_df.rename(columns = {'two':'new_name'})
    print("Sorting the genes...")
    gene_cn_df.sort_index(inplace=True)

    return gene_cn_df
