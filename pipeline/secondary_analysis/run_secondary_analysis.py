from secondary_analysis import SecondaryAnalysis
import argparse
import json


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--configfile", required=True)
args = parser.parse_args()

with open(args.configfile) as json_config_file:
    config = json.load(json_config_file)

print(config)


try:
    sample_name = config['sample_name']
    output_path = config['output_path']
    h5_path = config['h5_path']
    genes_path = config['genes_path']
    bins = config['bins_to_remove']
except Exception as e:
    print("Error while parsing the config")
    print(e)

sa = SecondaryAnalysis(sample_name=sample_name, output_path=output_path, h5_path=h5_path, genes_path=genes_path)
sa.remove_tenx_genomics_artifacts(bins=bins)
sa.apply_phenograph()

