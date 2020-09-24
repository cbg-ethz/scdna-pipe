# Title: Standard operating procedure single cell genomics data analysis
# Version: 1.14
# Date: 2020-08-11
# Author: Pedro Ferreira and Monica Dragan
# Reviewer: Jack Kuipers
# Approver: Niko Beerenwinkel

# Introduction:
This SOP describes the computational analysis performed on the scDNA data from the 10X Genomics platform.
All analyses are based on the reference genome hg19 with annotation GENCODE V25 (Ensembl V87).

# Analysis: CellRanger DNA
- Description:Cell Ranger DNA allows for the analysis of a single 10x-barcoded sequencing library prepared from a sample that is sequenced across one or more flowcells.
- Script/Software: CellRanger DNA (<https://support.10xgenomics.com/single-cell-dna)>
- Input: raw .fastq.gz files and the reference genome
## cnv_data.h5
- Format: HDF5
```code
FILE_CONTENTS {
group      /
dataset    /cell_barcodes
group      /cnvs
dataset    /cnvs/chr1
dataset    /cnvs/chr10
dataset    /cnvs/chr11
dataset    /cnvs/chr12
dataset    /cnvs/chr13
dataset    /cnvs/chr14
dataset    /cnvs/chr15
dataset    /cnvs/chr16
dataset    /cnvs/chr17
dataset    /cnvs/chr18
dataset    /cnvs/chr19
dataset    /cnvs/chr2
dataset    /cnvs/chr20
dataset    /cnvs/chr21
dataset    /cnvs/chr22
dataset    /cnvs/chr3
dataset    /cnvs/chr4
dataset    /cnvs/chr5
dataset    /cnvs/chr6
dataset    /cnvs/chr7
dataset    /cnvs/chr8
dataset    /cnvs/chr9
dataset    /cnvs/chrX
dataset    /cnvs/chrY
group      /constants
dataset    /constants/bin_size
dataset    /constants/chroms
dataset    /constants/num_bins_per_chrom
dataset    /constants/num_cells
dataset    /constants/num_chroms
dataset    /constants/num_nodes
group      /genome_tracks
group      /genome_tracks/gc_fraction
dataset    /genome_tracks/gc_fraction/chr1
dataset    /genome_tracks/gc_fraction/chr10
dataset    /genome_tracks/gc_fraction/chr11
dataset    /genome_tracks/gc_fraction/chr12
dataset    /genome_tracks/gc_fraction/chr13
dataset    /genome_tracks/gc_fraction/chr14
dataset    /genome_tracks/gc_fraction/chr15
dataset    /genome_tracks/gc_fraction/chr16
dataset    /genome_tracks/gc_fraction/chr17
dataset    /genome_tracks/gc_fraction/chr18
dataset    /genome_tracks/gc_fraction/chr19
dataset    /genome_tracks/gc_fraction/chr2
dataset    /genome_tracks/gc_fraction/chr20
dataset    /genome_tracks/gc_fraction/chr21
dataset    /genome_tracks/gc_fraction/chr22
dataset    /genome_tracks/gc_fraction/chr3
dataset    /genome_tracks/gc_fraction/chr4
dataset    /genome_tracks/gc_fraction/chr5
dataset    /genome_tracks/gc_fraction/chr6
dataset    /genome_tracks/gc_fraction/chr7
dataset    /genome_tracks/gc_fraction/chr8
dataset    /genome_tracks/gc_fraction/chr9
dataset    /genome_tracks/gc_fraction/chrX
dataset    /genome_tracks/gc_fraction/chrY
group      /genome_tracks/is_mappable
dataset    /genome_tracks/is_mappable/chr1
dataset    /genome_tracks/is_mappable/chr10
dataset    /genome_tracks/is_mappable/chr11
dataset    /genome_tracks/is_mappable/chr12
dataset    /genome_tracks/is_mappable/chr13
dataset    /genome_tracks/is_mappable/chr14
dataset    /genome_tracks/is_mappable/chr15
dataset    /genome_tracks/is_mappable/chr16
dataset    /genome_tracks/is_mappable/chr17
dataset    /genome_tracks/is_mappable/chr18
dataset    /genome_tracks/is_mappable/chr19
dataset    /genome_tracks/is_mappable/chr2
dataset    /genome_tracks/is_mappable/chr20
dataset    /genome_tracks/is_mappable/chr21
dataset    /genome_tracks/is_mappable/chr22
dataset    /genome_tracks/is_mappable/chr3
dataset    /genome_tracks/is_mappable/chr4
dataset    /genome_tracks/is_mappable/chr5
dataset    /genome_tracks/is_mappable/chr6
dataset    /genome_tracks/is_mappable/chr7
dataset    /genome_tracks/is_mappable/chr8
dataset    /genome_tracks/is_mappable/chr9
dataset    /genome_tracks/is_mappable/chrX
dataset    /genome_tracks/is_mappable/chrY
group      /genome_tracks/mappability
dataset    /genome_tracks/mappability/chr1
dataset    /genome_tracks/mappability/chr10
dataset    /genome_tracks/mappability/chr11
dataset    /genome_tracks/mappability/chr12
dataset    /genome_tracks/mappability/chr13
dataset    /genome_tracks/mappability/chr14
dataset    /genome_tracks/mappability/chr15
dataset    /genome_tracks/mappability/chr16
dataset    /genome_tracks/mappability/chr17
dataset    /genome_tracks/mappability/chr18
dataset    /genome_tracks/mappability/chr19
dataset    /genome_tracks/mappability/chr2
dataset    /genome_tracks/mappability/chr20
dataset    /genome_tracks/mappability/chr21
dataset    /genome_tracks/mappability/chr22
dataset    /genome_tracks/mappability/chr3
dataset    /genome_tracks/mappability/chr4
dataset    /genome_tracks/mappability/chr5
dataset    /genome_tracks/mappability/chr6
dataset    /genome_tracks/mappability/chr7
dataset    /genome_tracks/mappability/chr8
dataset    /genome_tracks/mappability/chr9
dataset    /genome_tracks/mappability/chrX
dataset    /genome_tracks/mappability/chrY
group      /genome_tracks/n_fraction
dataset    /genome_tracks/n_fraction/chr1
dataset    /genome_tracks/n_fraction/chr10
dataset    /genome_tracks/n_fraction/chr11
dataset    /genome_tracks/n_fraction/chr12
dataset    /genome_tracks/n_fraction/chr13
dataset    /genome_tracks/n_fraction/chr14
dataset    /genome_tracks/n_fraction/chr15
dataset    /genome_tracks/n_fraction/chr16
dataset    /genome_tracks/n_fraction/chr17
dataset    /genome_tracks/n_fraction/chr18
dataset    /genome_tracks/n_fraction/chr19
dataset    /genome_tracks/n_fraction/chr2
dataset    /genome_tracks/n_fraction/chr20
dataset    /genome_tracks/n_fraction/chr21
dataset    /genome_tracks/n_fraction/chr22
dataset    /genome_tracks/n_fraction/chr3
dataset    /genome_tracks/n_fraction/chr4
dataset    /genome_tracks/n_fraction/chr5
dataset    /genome_tracks/n_fraction/chr6
dataset    /genome_tracks/n_fraction/chr7
dataset    /genome_tracks/n_fraction/chr8
dataset    /genome_tracks/n_fraction/chr9
dataset    /genome_tracks/n_fraction/chrX
dataset    /genome_tracks/n_fraction/chrY
group      /normalized_counts
dataset    /normalized_counts/chr1
dataset    /normalized_counts/chr10
dataset    /normalized_counts/chr11
dataset    /normalized_counts/chr12
dataset    /normalized_counts/chr13
dataset    /normalized_counts/chr14
dataset    /normalized_counts/chr15
dataset    /normalized_counts/chr16
dataset    /normalized_counts/chr17
dataset    /normalized_counts/chr18
dataset    /normalized_counts/chr19
dataset    /normalized_counts/chr2
dataset    /normalized_counts/chr20
dataset    /normalized_counts/chr21
dataset    /normalized_counts/chr22
dataset    /normalized_counts/chr3
dataset    /normalized_counts/chr4
dataset    /normalized_counts/chr5
dataset    /normalized_counts/chr6
dataset    /normalized_counts/chr7
dataset    /normalized_counts/chr8
dataset    /normalized_counts/chr9
dataset    /normalized_counts/chrX
dataset    /normalized_counts/chrY
group      /raw_counts
dataset    /raw_counts/chr1
dataset    /raw_counts/chr10
dataset    /raw_counts/chr11
dataset    /raw_counts/chr12
dataset    /raw_counts/chr13
dataset    /raw_counts/chr14
dataset    /raw_counts/chr15
dataset    /raw_counts/chr16
dataset    /raw_counts/chr17
dataset    /raw_counts/chr18
dataset    /raw_counts/chr19
dataset    /raw_counts/chr2
dataset    /raw_counts/chr20
dataset    /raw_counts/chr21
dataset    /raw_counts/chr22
dataset    /raw_counts/chr3
dataset    /raw_counts/chr4
dataset    /raw_counts/chr5
dataset    /raw_counts/chr6
dataset    /raw_counts/chr7
dataset    /raw_counts/chr8
dataset    /raw_counts/chr9
dataset    /raw_counts/chrX
dataset    /raw_counts/chrY
group      /tree
dataset    /tree/Z
group      /tree/heterogeneity
dataset    /tree/heterogeneity/chr1
dataset    /tree/heterogeneity/chr10
dataset    /tree/heterogeneity/chr11
dataset    /tree/heterogeneity/chr12
dataset    /tree/heterogeneity/chr13
dataset    /tree/heterogeneity/chr14
dataset    /tree/heterogeneity/chr15
dataset    /tree/heterogeneity/chr16
dataset    /tree/heterogeneity/chr17
dataset    /tree/heterogeneity/chr18
dataset    /tree/heterogeneity/chr19
dataset    /tree/heterogeneity/chr2
dataset    /tree/heterogeneity/chr20
dataset    /tree/heterogeneity/chr21
dataset    /tree/heterogeneity/chr22
dataset    /tree/heterogeneity/chr3
dataset    /tree/heterogeneity/chr4
dataset    /tree/heterogeneity/chr5
dataset    /tree/heterogeneity/chr6
dataset    /tree/heterogeneity/chr7
dataset    /tree/heterogeneity/chr8
dataset    /tree/heterogeneity/chr9
dataset    /tree/heterogeneity/chrX
dataset    /tree/heterogeneity/chrY
dataset    /tree/is_cell_in_group
    }
}
```
## web_summary.html
- Format: Webpage
- Content description: Summary statistics of the cellranger run (basic cell and read counting), including figures.
## summary.csv
- Format: Comma-separated table including header.
- Content description: Summary statistics of the cellranger run (basic cell and read counting).
## alarms_summary.txt
- Format: Unstructured human readable data
- Content description: Text file containing the sequencing alerts.
## possorted_bam.bam
- Format: BAM file
- Content description: Barcode-corrected reads aligned to the reference, sorted by reference position.

# Analysis: Extract genomic coordinates
- Description: Extracts the genomic coordinates from cellranger output
- Input: List of genes of interest and cnv_data.h5
## chr_stops.tsv
- Format: Tab-separated value
- Content description: contains chromosome names and their corresponding bin positions
## bins_genome.tsv
- Format: Tab-separated value
- Content description: Rows correspond to the bins in the genome. Column "bin_ids" denotes the bin ids. Column "start" is the starting genomic position of that bin. Column "end" is the end genomic position of that bin.

# Analysis: Filtering bins and cells
- Description: Removal of low quality bins and unmappable bins (using umbilical cord data as a reference) and cells with outlier DIMAPD statistic.
- Input: list of low quality bins and cnv_data.h5
## filtered_counts.csv
- Format: Comma-separated value
- Content description: Filtered normalized count matrix (counts per bin per cell).
## excluded_bins.csv
- Format: Comma-separated value
- Content description: Boolean mask for the bins to exclude

# Analysis: Breakpoint Detection
- Description: Detects breakpoints on the whole genome in order to create informative genomic regions
- Input: filtered_counts.csv
## segmented_regions.txt
- Format: One value per line
- Content description: Ids of the breakpoints
## segmented_region_sizes.txt
- Format: One value per line
- Content description: Sizes information per region
## segmented_counts.csv
- Format: Comma-separated value
- Content description: cells per region matrix

# Analysis: Normalisation
- Description: Normalises the counts matrices per rows (cells on rows)
- Input: filtered_counts.csv, segmented_counts.csv
## normalised_bins.csv
- Format: Comma-separated value
- Content description: normalised counts per cell per bin
## normalised_regions.csv
- Format: Comma-separated value
- Content description: normalised counts per cell per region with chromosome annotations

# Analysis: Cluster tree
- Description: Cluster cells based on the similarity of their filtered normalized counts. Number of neighbours parameter is derived from the number of cells divided by 10 and rounded to integer. This is used to learn the Phylogenetic tree on the cluster averages.
- Script/Software: SCICoNE (<https://github.com/cbg-ethz/SCICoNE>)
- Input: filtered_counts.tsv, segmented_region_sizes.txt

## clustering_score.txt
- Format: Text file with a single floating point number.
- Content description: The Louvain Modularity score of the graph-based clustering algorithm.
## cluster_tree.txt
- Format: Custom human-readable text format.
- Content description: represents the phylogenetic tree structure and its mutations per node.
## cluster_tree.json
- Format: JSON
- Content description: represents the phylogenetic tree structure and its mutations per node.
## unique_cluster_tree_cnvs.csv
- Format: Comma-separated value
- Content description: the unique inferred copy number values per cell per bin
## inferred_cnvs.csv
- Format: Comma-separated value
- Content description: the copy number values inferred per cell per bin, including the filtered out bins.
This step is required to map the genomic information back to the inferred results.
## cluster_profile_files.txt
- Type: file list
- Format: Portable Network Graphics (PNG)
- Content description: The inferred copy number state plot for each cluster across all bins. The x axis of the plot contains chromosome annotations.
## cluster_profile_overlapping.png
- Format: Portable Network Graphics (PNG)
- Content description: The overlapping plot of inferred copy number states across bins for all clusters. The x axis of the plot contains chromosome annotations.
## cn_gene_df.csv
- Format: Comma separated value
- Content description: contains the copy number values per disease-specific genes of interest.
## heatmap_cnvs.png
- Format: Portable Network Graphics (PNG)
- Content description: Heatmap showing the average copy number states per cluster per chromosome for disease-specific genes of interest.
## cn_gene_df_roche_gene_list.csv
- Format: Comma separated value
- Content description: contains the copy number values per genes of interest provided by Roche (available at <https://github.com/cbg-ethz/scdna-pipe/blob/master/required_files/genes_of_interest/roche_gene_list.txt>).
## cluster_tree_genes.png
- Format: Portable Network Graphics (PNG)
- Content description: Phylogenetic tree figure displaying all the drug info package genes affected by copy number alterations with disease-specific genes highlighted.
## cluster_tree_sorted_normalised_counts_bins.png
- Format: Portable Network Graphics (PNG)
- Content description: Plots the normalised counts matrix per cell per bin with chromosome annotations in the columns and clone annotations in the rows.
## cluster_tree_sorted_cnvs_bins.png
- Format: Portable Network Graphics (PNG)
- Content description: Plots the inferred CNV matrix per cell per bin with chromosome annotations in the columns and clone annotations in the rows.
## clone_lib_sizes.png
- Format: Portable Network Graphics (PNG)
- Content description: Plots the library sizes of each clone.
