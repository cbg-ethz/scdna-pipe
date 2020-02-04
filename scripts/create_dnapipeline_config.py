"""
Script to generate the config.json file used for running the TuPro DNA pipeline.

Usage:
python create_dnapipeline_config.py -s MADEGOD-T_scD_250c-r1v1.0_r1v1.0-AHNY5VBGXC -o config_MADEGOD.json -t melanoma
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
    "-s",
    "--labkey_sample_name",
    type=str,
    required=True,
    help="Name of one scD sample as reported in LabKey, e.g., MADEGOD-T_scD_250c-r1v1.0_r1v1.0-AHNY5VBGXC")
parser.add_argument(
    "-t",
    "--cancer_type",
    type=str,
    required=True,
    help="cancer type, e.g. melanoma")
parser.add_argument(
    "-v",
    "--pipeline_version",
    type=str,
    default='v1.12',
    help='Version of the DNA pipeline')
parser.add_argument("-a", "--analysis_path", type=str,
                    default='/cluster/work/bewi/ngs/projects/tumorProfiler/analysis', help='Path to the analisys directory.')
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
pattern_1 = re.search('^(.+)-T_', args.labkey_sample_name)
pattern_2 = re.search('-([^-]+)$', args.labkey_sample_name)
pattern_3 = re.search('scD_250c-r1v1.0', args.labkey_sample_name)
if not (pattern_1 and pattern_2 and pattern_3):
    sys.exit(
        'Argument \"--labkey_sample_name ' +
        args.labkey_sample_name +
        '\" does not have the expected pattern.')
sample_name = pattern_1.group(1)
sample_annotation = pattern_2.group(1)

# Prepare the folder structure: analysis
singlecell_dna_path = os.path.join(
    args.analysis_path, "trial_" + args.cancer_type, sample_name + "-T", "singlecell_dna/")
analysis_path = os.path.join(singlecell_dna_path, "analysis")
code_path = os.path.join(args.project_path, "code/dna-pipeline/sc-dna/bin")

# Build the json config
config = {}

config['sample_name'] = sample_name + "_S2"
config['sequencing_prefix'] = sample_name + \
    "-T_scD_250c-r1v1.0_r1v1.0-" + sample_annotation
config['analysis_prefix'] = sample_name + \
    "-T" + "_scD_Ar1" + args.pipeline_version
config['disease'] = args.cancer_type
config['ref_genome_version'] = "GRCh37"
config['ref_genome_path'] = os.path.join(
    args.project_path, "data/refdata-GRCh37-1.0.0_dna")
config['fastqs_path'] = os.path.join(singlecell_dna_path, "openbis")
config['analysis_path'] = os.path.join(analysis_path)
config['scripts_dir'] = "/cluster/work/bewi/members/pedrof/dna-pipeline/scripts/"
config['10x_artifacts'] = "/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/10x_artifacts"
config['bin_size'] = 20000
config['tricking_fastqs'] = {"mem": 1000, "time": 1438}
config['cellranger_dna'] = {
    "local_cores": 48,
    "local_mem": 192,
    "mem_per_core": 4,
    "mem": 4096,
    "time": 1438}
config['secondary_analysis'] = {"h5_path": os.path.join(analysis_path, "cellranger", sample_name, "outs/cnv_data.h5"),
                                "bins_to_remove": "/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/10x_artifacts/GRCh37_10x_artifacts_mask_whole_genome.tsv",
                                "genes_path": "/cluster/work/bewi/members/pedrof/dna-pipeline/required_files/genes_of_interest/",
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
                                  "bin": os.path.join(code_path, "breakpoint_detection")
                                  }
config['inference'] = {"ploidy": 2,
                       "verbosity": 1,
                       "copy_number_limit": 2,
                       "bin": os.path.join(code_path, "inference"),
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
