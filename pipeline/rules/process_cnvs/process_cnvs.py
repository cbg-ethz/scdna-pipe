import numpy as np
import pandas as pd
from tqdm import tqdm

def get_bin_gene_region_df(bin_size, gene_coordinates, chr_stops, region_stops, bin_is_excluded, cnvs=None, priority_genes=None):
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

    bin_gene_region_df = pd.DataFrame(index=range(chr_stops.iloc[-1].values[0]+1))
    bin_gene_region_df["gene"] = None
    bin_gene_region_df["chr"] = None
    bin_gene_region_df["region"] = None
    bin_gene_region_df["is_priority"] = None
    if priority_genes is not None:
        bin_gene_region_df["gene"] = [list() for _ in range(bin_gene_region_df.shape[0])]
        bin_gene_region_df["chr"] = [list() for _ in range(bin_gene_region_df.shape[0])]
        bin_gene_region_df["is_priority"] = [list() for _ in range(bin_gene_region_df.shape[0])]

    # for each gene
    for index, row in tqdm(gene_coordinates.iterrows(), total=gene_coordinates.shape[0]):
        start_bin = int(row["Gene start (bp)"] / bin_size)
        stop_bin = int(row["Gene end (bp)"] / bin_size)
        chromosome = str(row["Chromosome/scaffold name"])

        if chromosome != '1': # coordinates are given by chromosome
            chr_start = chr_stops.iloc[np.where(chr_stops.index==chromosome)[0][0]-1].values[0] + 1
            start_bin = start_bin + chr_start
            stop_bin = stop_bin + chr_start

        gene_name = row["Gene name"]
        for bin in range(start_bin, stop_bin+1):
            bin_gene_region_df.loc[bin, "gene"].append(gene_name)
            bin_gene_region_df.loc[bin, "chr"].append(chromosome)

        if priority_genes is not None:
            if gene_name in priority_genes:
                for bin in range(start_bin, stop_bin+1):
                    bin_gene_region_df.loc[bin, "is_priority"].append(True)
                print(bin_gene_region_df.loc[start_bin:stop_bin+1, "gene"])
            else:
                # Need this to deal with overlapping gene coordinates
                for bin in range(start_bin, stop_bin+1):
                    bin_gene_region_df.loc[bin, "is_priority"].append(False)
        else:
            for bin in range(start_bin, stop_bin+1):
                bin_gene_region_df.loc[bin, "is_priority"].append(True)

    # Turn columns of lists into columns of strings with comma-separated values
    bin_gene_region_df['gene'] = [','.join(map(str, l)) for l in bin_gene_region_df['gene']]
    bin_gene_region_df['chr'] = [','.join(map(str, l)) for l in bin_gene_region_df['chr']]
    bin_gene_region_df['is_priority'] = [','.join(map(str, l)) for l in bin_gene_region_df['is_priority']]

    # Indicate original_bin-filtered_bin correspondence
    bin_gene_region_df['filtered_bin'] = None
    bin_gene_region_df['filtered_bin'].iloc[np.where(np.array(bin_is_excluded)==False)[0]] = np.arange(np.count_nonzero(np.array(bin_is_excluded)==False))

    # Get the regions
    start_bin = 0
    for index, row in tqdm(region_stops.iterrows(), total=region_stops.shape[0]):
        stop_bin = row.values[0]
        original_start_bin = np.where(bin_gene_region_df['filtered_bin'] == start_bin)[0][0]
        original_stop_bin = np.where(bin_gene_region_df['filtered_bin'] == stop_bin)[0][0]
        bin_gene_region_df.loc[original_start_bin:original_stop_bin+1, 'region'] = index # regions are 0 indexed
        start_bin = row.values[0]

    if cnvs is not None:
        # Add the CN values
        if len(cnvs.shape) == 1:
            cnvs = cnvs.reshape(1, -1)

        for c_id in range(cnvs.shape[0]):
            bin_gene_region_df["cnv_{}".format(c_id)] = cnvs[c_id]

    # Make sure excluded bins have no info
    bin_gene_region_df.loc[np.where(bin_is_excluded)[0], 'region'] = None

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
    gene_lists = df['gene'][np.where(df['region']==region)[0]].values.tolist()
    gene_lists = [sublist.split(',') for sublist in gene_lists]
    gene_list = [item for sublist in gene_lists for item in sublist]

    if priority_only:
        is_priority_lists = df['is_priority'][np.where(df['region']==region)[0]].values.tolist()
        is_priority_lists = [sublist.split(',') for sublist in is_priority_lists]
        is_priority_list = [item for sublist in is_priority_lists for item in sublist]

        # Subset only the priority genes
        priority_gene_list = []
        for i, gene in enumerate(gene_list):
            if is_priority_list[i] == 'True':
                priority_gene_list.append(gene)
        gene_list = priority_gene_list

    # Remove duplicates
    gene_list = np.unique(gene_list).tolist()

    return gene_list

def get_region_with_gene(gene, bin_gene_region_df):
    """
        Returns the region index containing a gene
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: region index (integer)
    """
    bins_with_gene = np.where([gene in l for l in bin_gene_region_df['gene']])[0]
    region = np.array(bin_gene_region_df.loc[bins_with_gene,'region'].values)

    if len(region[region!=None]) > 0:
        # Get first non-None value
        region = region[region!=None][0]
    else:
        region = None

    return region

def convert_event_region_to_gene(region_event_str, bin_gene_region_df, priority_only=False, genes_to_highlight=None, highlight_color='red', genes_only=False):
    """
        Returns a string indicating gene-wise events in affected region
        Examples:
                "+2R174"     -> ["+2BRAF", "+2MALAT1"]
                "-1R656:658" -> ["-1JAK2", "-1MLNA", "-1CDK4"]
        :param region_event_str: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param priority_only: boolean
            indicating if only priority genes should be returned
        :param genes_to_highlight: list
            genes that should be displayed in a different color
        :param highlight_color: str
            color to use in genes to highlight
        :return: list of str
    """
    # Get event (-2, -1, +1, +2, etc)
    event_str = region_event_str[:2]
    region_str = region_event_str[3:]
    if ":" in region_str: # multiple regions: "-1R656:658"
        aux = [int(region) for region in region_str.split(":")]
        region_list = np.arange(aux[0], aux[1]+1)
    else:
        region_list = [int(region_str)]

    gene_list = []
    for region in region_list:
        genes_in_region = get_genes_in_region(region, bin_gene_region_df, priority_only=priority_only)
        gene_list.append(genes_in_region)

    gene_list = [item for sublist in gene_list for item in sublist]

    # Highlight some genes
    if genes_to_highlight is not None:
        for index, gene in enumerate(gene_list):
            if gene in genes_to_highlight:
                gene_list[index] = "<font color=" + "\'" + highlight_color + "\'" + ">" + gene + "</font>"

    gene_string = '[' + ','.join(gene_list) + ']'
    if len(gene_list) == 0:
        gene_event_str = ""
    else:
        if not genes_only:
            gene_event_str = event_str + gene_string
        else:
            gene_event_str = gene_string

    return gene_event_str

def convert_node_regions_to_genes(node_str, bin_gene_region_df, priority_only=False, genes_to_highlight=None, highlight_color='red', max_genes_per_line=10):
        """
            Returns a string indicating gene events and total number of
            amplifications and deletions in node
            Examples:
                    "+2R174 +2R291:293 -1R656:658" -> "+2[BRAF,MALAT1] -1[JAK2,MLNA,CDK4]\n(3+, 1-)"
            :param node_str: str
            :param bin_gene_region_df: DataFrame
                with (gene_name, region, original_bin) fields
            :param priority_only: boolean
                indicating if only priority genes should be returned
            :param genes_to_highlight: list
                genes that should be displayed in a different color
            :return: str
        """
        region_event_strs = node_str.split(' ')
        num_events = len(region_event_strs)

        num_amplifications = 0

        str_dict = dict()
        possible_events = ['+4', '+3', '+2', '+1', '-1', '-2', '-3', '+4']
        for event in possible_events:
            str_dict[event] = []

        for region_event_str in region_event_strs:
            # Get event (-2, -1, +1, +2, etc)
            event_str = region_event_str[:2]
            region_str = region_event_str[3:]
            if ":" in region_str: # multiple regions: "-1R656:658"
                aux = [int(region) for region in region_str.split(":")]
                region_list = np.arange(aux[0], aux[1]+1)
            else:
                region_list = [int(region_str)]

            gene_list = []
            for region in region_list:
                genes_in_region = get_genes_in_region(region, bin_gene_region_df, priority_only=priority_only)
                gene_list.append(genes_in_region)

            gene_list = [item for sublist in gene_list for item in sublist]

            str_dict[event_str].append(gene_list)

            if region_event_str[0]=='+':
                num_amplifications += 1

        for key in str_dict:
            str_dict[key] = [item for sublist in str_dict[key] for item in sublist]

        gene_event_str = ''
        newline = '<br/>'
        for key in str_dict:
            if len(str_dict[key]) != 0:
                gene_event_str = gene_event_str + '<B>' + key + '</B>'
                gene_event_str = gene_event_str + '[' + ','.join(f"{x}"+newline if (i+1)%max_genes_per_line == 0 else str(x) for i, x in enumerate(str_dict[key])) + ']'
                gene_event_str = gene_event_str + " " + newline + " " + newline

        num_deletions = num_events - num_amplifications
        num_events, num_amplifications, num_deletions
        num_events_str = "{} +, {} -".format(num_amplifications, num_deletions)

        node_str = gene_event_str + "(" + num_events_str + ")"

        # If newline followed by ']', replace with newline after ']'
        node_str = node_str.replace(newline+"]", "]"+newline)

        # If newline followed by ',', replace with newline after ','
        node_str = node_str.replace(newline+",", ","+newline)

        # If newline followed by newline, remove one
        for m in re.finditer(newline, node_str):
            index = m.start()
            if node_str[index+len(newline):index+2*len(newline)] == newline:
                node_str = "".join((node_str[:index+len(newline)], "", node_str[index+2*len(newline):]))

        # highlight genes
        if genes_to_highlight is not None:
            for gene in genes_to_highlight:
                node_str = node_str.replace(gene, "<font color=" + "\'" + highlight_color + "\'" + ">" + gene + "</font>")

        if gene_event_str == "":
            node_str = num_events_str

        return node_str

def get_region_with_gene(gene, bin_gene_region_df):
    """
        Returns the region index containing a gene
        :param gene: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :return: region index (integer)
    """
    bins_with_gene = np.where([gene in l for l in bin_gene_region_df['gene']])[0]
    region = np.array(bin_gene_region_df.loc[bins_with_gene,'region'].values)

    if len(region[~np.isnan(region)]) > 0:
        # Get first non-None value
        region = int(region[~np.isnan(region)][0])
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
    bin_gene_region_df['gene'] = bin_gene_region_df['gene'].astype(str).apply(lambda x: x.split(',')).apply(lambda x: set(x))
    bins = bin_gene_region_df[bin_gene_region_df['gene'].apply(lambda x: gene in x)].index.values

    left_bin = bins[0]
    right_bin = bins[-1]

    left_regions = bin_gene_region_df['region'].values[:left_bin]
    left_region = left_regions[~np.isnan(left_regions)][-1]

    right_regions = bin_gene_region_df['region'].values[right_bin+1:]
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
    n_subclones = np.count_nonzero(['cnv' in column for column in bin_gene_region_df.columns])
    cluster_ids = range(n_subclones)
    gene_cn_df = pd.DataFrame(index=cluster_ids)

    # for each gene
    for gene in gene_list:
        gene_cn_per_cluster = []
        for c_id in cluster_ids:
            median_cn = np.nanmedian(
                bin_gene_region_df['cnv_{}'.format(c_id)][bin_gene_region_df['gene']==gene].values
            )

            # If NaN, impute with median value of regions surrounding it
            if np.isnan(median_cn) and impute:
                # Get regions surrounding gene
                left_region, right_region = get_surrounding_regions(gene, bin_gene_region_df)

                # get CNV values of region surronding gene
                left_cn = bin_gene_region_df['cnv_{}'.format(c_id)][bin_gene_region_df['region']==left_region].iloc[0]
                right_cn = bin_gene_region_df['cnv_{}'.format(c_id)][bin_gene_region_df['region']==right_region].iloc[0]

                median_cn = np.nanmedian([left_cn, right_cn])

            if median_cn > 2:
                median_cn = int(np.floor(median_cn))
            elif median_cn < 2:
                median_cn = int(np.ceil(median_cn))

            gene_cn_per_cluster.append(median_cn)

        gene_cn_df[gene] = gene_cn_per_cluster

    print("Transposing the dataframe...")
    gene_cn_df = gene_cn_df.T
    print("Sorting the genes...")
    gene_cn_df.sort_index(inplace=True)

    return gene_cn_df
