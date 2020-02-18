"""
Script to generate the config.json file used for running the TuPro DNA pipeline.
Usage:
python create_dnapipeline_config.py -f BSSE_QGF_131049_HNY5VBGXC_1_MADEGOD_T_scD_250c_r1v1_0_SI-GA-E9_S2_L001_I1_001 -o config_MADEGOD.json -t melanoma
"""

import argparse
import json
import os
import pprint
import re
import sys

parser = argparse.ArgumentParser(
    description='Generate the json config for the setup of the dna pipeline')
parser.add_argument(
    "-f",
    "--openbis_fastq_filename",
    type=str,
    required=True,
    help="Filename of one of the openBIS fastq files, with or without extension, e.g., BSSE_QGF_131049_HNY5VBGXC_1_MADEGOD_T_scD_250c_r1v1_0_SI-GA-E9_S2_L001_I1_001")
parser.add_argument(
    "-t",
    "--cancer_type",
    type=str,
    required=True,
    help="cancer type, e.g. melanoma")
parser.add_argument(
    "-n",
    "--is_novaseq",
    action='store_true',
    help='True when the data was sequenced with NovaSeq (default). When False, we expect data sequenced with NextSeq.')
parser.add_argument(
    "-v",
    "--pipeline_version",
    type=str,
    default='v1.12',
    help='Version of the DNA pipeline')
parser.add_argument("-a", "--analysis_path", type=str,
                    default='/cluster/work/bewi/ngs/projects/tumorProfiler/analysis', help='Path to the analisys directory, where the results and intermediate files will be stored.')
parser.add_argument("-p", "--project_path", type=str,
                    default='/cluster/work/bewi/ngs/projects/tumorProfiler/', help='Project path, where all the code and tools are set up.')
parser.add_argument("-o", "--out", type=str,
                    default="config.json", help='Path of the directory where to store the output config. By default it will be saved in the current directory.')
args = parser.parse_args()

# Argument validation
if not re.match(r'v\d+\.\d+', args.pipeline_version):
    sys.exit(
        'Argument \"--pipeline_version ' +
        args.pipeline_version +
        '\" does not have the expected pattern.')

# Parse the input
pattern_1 = re.search('_([^_]+)_T_', args.openbis_fastq_filename)
pattern_2 = re.search('^BSSE_QGF_[0-9]+_(.+)_', args.openbis_fastq_filename)
pattern_3 = re.search('_S([0-9]+)_L00[0-9]_.._001$', args.openbis_fastq_filename)
if not (pattern_1 and pattern_2 and pattern_3):
    sys.exit(
        'Argument \"--openbis_fastq_filename ' +
        args.openbis_fastq_filename +
        '\" does not have the expected pattern.')
sample_name = pattern_1.group(1)
sample_annotation = pattern_2.group(1)
sample_number = pattern_3.group(1)

singlecell_dna_path = os.path.join(
    args.analysis_path, "trial_" + args.cancer_type, sample_name + "-T", "singlecell_dna/")
analysis_path = os.path.join(singlecell_dna_path, "analysis")
dna_pipeline_code_path = os.path.join(args.project_path, "code/dna-pipeline-novaseq")
sc_dna_code_path = os.path.join(dna_pipeline_code_path, "sc-dna/bin")

# Build the json config
config = {}

config['sample_name'] = sample_name + "_S" + sample_number
config['sequencing_prefix'] = sample_name + \
    "-T_scD_250c-r1v1.0_r1v1.0-A" + sample_annotation
config['analysis_prefix'] = sample_name + \
    "-T" + "_scD_Ar1" + args.pipeline_version
config['disease'] = args.cancer_type
config['ref_genome_version'] = "GRCh37"
config['ref_genome_path'] = os.path.join(
    args.project_path, "data/refdata-GRCh37-1.0.0_dna")
config['fastqs_path'] = os.path.join(singlecell_dna_path, "openbis")
config['analysis_path'] = os.path.join(analysis_path)
config['scripts_dir'] = os.path.join(dna_pipeline_code_path, "scripts/")
config['10x_artifacts'] = os.path.join(dna_pipeline_code_path, "required_files/10x_artifacts")
config['bin_size'] = 20000

# By default we expect data sequenced with NovaSeq.
config['n_lanes'] = 2
insert_length = 91
if not args.is_novaseq:
  config['n_lanes'] = 4
  insert_length = 58
config['tricking_fastqs'] = {
    "insert_length": insert_length,
    "mem": 1000,
    "time": 1438}

config['cellranger_dna'] = {
    "local_cores": 48,
    "local_mem": 192,
    "mem_per_core": 4,
    "mem": 4096,
    "time": 1438}
config['secondary_analysis'] = {"h5_path": os.path.join(analysis_path, "cellranger", sample_name, "outs/cnv_data.h5"),
                                "bins_to_remove": os.path.join(dna_pipeline_code_path, "required_files/10x_artifacts/GRCh37_10x_artifacts_mask_whole_genome.tsv"),
                                "genes_path": os.path.join(dna_pipeline_code_path, "required_files/genes_of_interest/"),
                                "general_main_gene_list": "dna_long_gene_list.txt"}

config['plotting'] = {"profiles": {
    "offset_sizes": 0.1,
    "s": 5
},
    "trees": {
    "highlight_color": "yellow",
    "max_genes_per_line": 6
}
}
config['breakpoint_detection'] = {"window_size": 50,
                                  "verbosity": 1,
                                  "threshold": 2,
                                  "bp_limit": 200,
                                  "bin": os.path.join(sc_dna_code_path, "breakpoint_detection")
                                  }
config['inference'] = {"ploidy": 2,
                       "verbosity": 1,
                       "copy_number_limit": 2,
                       "bin": os.path.join(sc_dna_code_path, "inference"),
                       "robustness_thr": 0.3,
                       "learn_nu":
                       {
                           "n_iters": 2000,
                           "n_nodes": 0,
                           "move_probs": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
                       },
                       "cluster_trees":
                       {
                           "n_iters": 500000,
                           "n_nodes": 0,
                           "move_probs": [0, 1, 0, 1, 0, 10, 0, 1, 0, 1, 0.4, 1, 0.01],
                           "n_reps": 10,
                           "alpha": 0.
                       },
                       "learn_nu_cluster_trees":
                       {
                           "n_iters": 20000,
                           "move_probs": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
                       },
                       "full_trees":
                       {
                           "n_iters": 1,
                           "n_nodes": 0,
                           "move_probs": [0, 1, 0, 1, 0, 10, 0, 1, 0, 1, 0.4, 0, 0.01],
                           "n_reps": 10,
                           "cluster_fraction": 1
                       },
                       "seed": 41
                       }

with open(args.out, 'w') as outfile:
    outfile.write(json.dumps(config, indent=2, sort_keys=False))
