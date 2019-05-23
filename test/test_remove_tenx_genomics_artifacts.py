import pytest
import sys
import os
import logging
import json
logging.basicConfig(level=logging.DEBUG)

dirname = os.path.dirname(__file__)
logging.debug('__file__ variable: ' +  __file__)
logging.debug('dirname: ' + dirname)
scripts_path = os.path.join(dirname, '../pipeline/secondary_analysis')
sys.path.append(scripts_path)
from secondary_analysis import SecondaryAnalysis


def test_remove_tenx_genomics_artifacts():
    """
    Tests the artifact removal
    :return:
    """

    config_file = './test_config.json'
    with open(config_file) as json_config_file:
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
    assert(sa.filtered_normalized_counts.shape == sa.filtered_cnvs.shape)

