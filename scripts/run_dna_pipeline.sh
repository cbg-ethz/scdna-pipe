#!/bin/bash

# Load the necessary modules and launch the DNA pipeline.
#
# Example: ./run_dna_pipeline.sh MADEGOD melanoma
#

echo "Usage: ./run_dna_pipeline.sh base_path sample_id cancer_type"

if [ $# -ne 2 ]
then
  echo "Please provide two arguments."
  exit 1
fi

module load python/3.6.0 graphviz hdf5/1.8.12 gcc/6.2.0 eth_proxy
module load /cluster/work/bewi/modules/cellranger-dna/1.1.0/cellranger-dna-1.1.0/cellranger-dna.modulefile

cd /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_${2}/${1}-T/singlecell_dna/snake_analysis_files/

bsub -J ${1} -n 48 -W 23:57 -R fullnode -R "rusage[mem=5000]" snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline-novaseq/pipeline/Snakefile --configfile ./config.json -j 48 -p -k 
