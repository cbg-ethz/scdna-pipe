from secondary_analysis import SecondaryAnalysis
import argparse
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--configfile", required=True)
args = parser.parse_args()

with open(args.configfile) as json_config_file:
    config = json.load(json_config_file)

print(config)


try:
    sample_name = config['sample_name']
    output_path = config['secondary_analysis']['output_path']
    h5_path = config['secondary_analysis']['h5_path']
    genes_path = config['secondary_analysis']['genes_path']
    all_genes_path = config['secondary_analysis']['all_genes_path']
    bins = config['secondary_analysis']['bins_to_remove']
except Exception as e:
    print("Error while parsing the config")
    print(e)
    sys.exit(1)

sa = SecondaryAnalysis(sample_name=sample_name, output_path=output_path, h5_path=h5_path, genes_path=genes_path, all_genes_path=all_genes_path)
sa.remove_tenx_genomics_artifacts(bins=bins)
sa.apply_phenograph()
sa.plot_clusters()
sa.plot_heatmap()
sa.create_cn_cluster_h5()

