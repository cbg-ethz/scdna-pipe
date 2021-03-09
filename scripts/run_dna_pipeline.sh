#!/bin/bash

# Load the necessary modules and launch the DNA pipeline.
#
# Example: ./run_dna_pipeline.sh MADEGOD melanoma
#

echo "Usage: ./run_dna_pipeline.sh sample_id"

if [ $# -ne 1 ]
then
  echo "Please provide sample_id."
  exit 1
fi

module load gcc/6.3.0
module load /cluster/work/tumorp/analysis/beerenwinkellab/scDNA/code/module_files/cellranger-dna/1.1.0

cd /cluster/work/tumorp/analysis/beerenwinkellab/scDNA/analysis/${1}-T/snake_analysis_files/

bsub -J ${1} -n 48 -W 23:57 -R "rusage[mem=5000]" snakemake -s /cluster/work/tumorp/analysis/beerenwinkellab/scDNA/code/scdna-pipe/pipeline/Snakefile --configfile ./config.json -j 48 -p -k
