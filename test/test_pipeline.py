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


def test_pipeline():
    """
    Tests the artifact removal
    :return:
    """

    config_file = './test_config.json'
    with open(config_file) as json_config_file:
        config = json.load(json_config_file)

    print(config)

    try:
        sample_name = config['analysis_prefix']
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
    assert(sa.filtered_normalized_counts.shape == sa.filtered_cnvs.shape)
    sa.apply_phenograph()
    sa.plot_clusters()
    sa.plot_heatmap()
    # sa.create_cn_cluster_h5()


