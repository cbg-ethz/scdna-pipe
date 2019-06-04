# DNA-Pipe 
[![CircleCI](https://circleci.com/gh/anilbey/dna-pipeline.svg?style=svg&circle-token=7d59442470c38d05f7d1661a97da237d482684ef)](https://circleci.com/gh/anilbey/dna-pipeline)
[![License](http://img.shields.io/:license-Apache%202-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0.txt)




Single cell dna data analysis workflow

# Analysis: Filtering bins
- Description: Removal of low quality bins, using umbilical cord data as a reference.
- Script/Software: 
- Input: list of low quality bins and cnv_data.h5
## filtered_counts.tsv
- Format: Tab-separated file.
- Content description: Filtered normalized count matrix (counts per bin per cell).
## filtered_cnvs.tsv
- Format: Tab-separated file.
- Content description: Filtered copy number states matrix.
## bins_genome.tsv
- Format: Tab-separated table.
- Content description: Rows correspond to the bins in the genome. Column "bin_ids" denotes the bin ids. Column "start" is the starting genomic position of that bin. Column "end" is the end genomic position of that bin. 
## chr_stops.tsv
- Format: Tab-separated file. Column "chr" is denoting the name of the chromosome and the index is the stop position for the chromosomes.
- Content description: Chromosomes and their stop positions.
# Analysis: Clustering
- Description: Cluster cells based on the similarity of their filtered normalized counts. Number of neighbours parameter is derived from the number of cells divided by 10 and rounded to integer.
- Script/Software: PhenoGraph (https://github.com/anilbey/PhenoGraph) 
- Input: __filtered_counts.tsv
## clusters_phenograph_assignment.tsv
- Format: Tab-separated table presenting the phenograph-based clustering assignment. 2 columns indicate the cell barcode ("cell_barcode") and the cluster id ("cluster"). Rows correspond to cells.
## clusters_phenograph_count_profiles.tsv
- Format: Tab-separated table presenting the phenograph-based clustering profiles. Rows for each cluster. Columns correspond to bins and the cells (of the matrix) contain the average counts for the bin, across all of the cells belonging to the corresponding cluster. The last column "cluster_ids" denotes the id of the cluster. -1 value for cluster id refers to outliers.
## clusters_phenograph_cn_profiles.tsv
- Format: Tab-separated table presenting the phenograph-based clustering profiles. Rows for each cluster. Columns correspond to bins and the cells (of the matrix) contain the average copy number state for the bin, across all of the cells belonging to the corresponding cluster. The last column "cluster_ids" denotes the id of the cluster. -1 value for cluster id refers to outliers.
## clustering_score.txt
- Content description: The Louvain Modularity score of the graph-based clustering algorithm. 
- Format: Text file with a single floating point number.
## cluster_sizes.txt
- Content description: Represents the sizes of each cluster.
- Format: Single line dictionary format with key-value pairs. 
## cluster_profile_files.txt
- Type: file list
- Format: Portable Network Graphics (PNG)
- Content description: The average copy number state plot for each cluster across all bins. The x axis of the plot is labeled by the chromosome stop positions.
## cluster_profile_overlapping.png
- Format: Portable Network Graphics (PNG)
- Content description: The overlapping plot of average counts across bins for all clusters. The x axis of the plot is labeled by the chromosome stop positions.
## tsne_output.png
- Format: Portable Network Graphics (PNG)
- Content description: t-SNE plot of the cells based on Phenograph distance matrix and labeling.
# Analysis: Mapping to genes
- Description: Maps the inferred CN states to the genes of interest for each cluster. 
- Input: clusters_phenograph_cn_profiles.tsv and list of interesting genes for the disease at hand.
## cn_gene_cluster.tsv
- Format: Tab-separated file. Rows for each cluster id. Columns correspond to genes and the cells (of the matrix) contain the average copy number state for the gene, across all of the cells belonging to the corresponding cluster. The last column "cluster_ids" denotes the id of the cluster. -1 value for cluster id refers to outliers.
## cn_genes_clusters_heatmap.png
- Format: Portable Network Graphics (PNG)
- Content description: Heatmap showing the average copy number states per cluster for clinically relevant genes. 
