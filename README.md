# Single-cell DNA Analysis <img src="figures/logo.png" align="right" width="160">

[![CircleCI](https://circleci.com/gh/cbg-ethz/scdna-pipe.svg?style=svg&circle-token=60921152e3353ae7cd6d5e13d14158bcbde57973)](https://circleci.com/gh/cbg-ethz/scdna-pipe)
[![License](http://img.shields.io/:license-Apache%202-green.svg)](http://www.apache.org/licenses/LICENSE-2.0.txt)
[![Python Version](https://img.shields.io/badge/python-3-blue.svg)](https://img.shields.io/badge/python-3-blue.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## About
Reproducible Python pipeline for genomic data analysis. Performs single-cell copy number variation calling by learning the underlying tumour evolution history by state-of-the-art phylogenetic tree reconstruction method: SCICoNE. The pipeline is built using Python, Conda environment management system and the Snakemake workflow management system. The pipeline starts from the raw sequencing files and a settings file for the parameter configurations. After the analysis, it produces a report and multiple figures to inform the treatment decision of the cancer patient.

The pipeline makes use of `scgenpy`, a package that exposes functions for preprocessing, postprocessing and plotting data, allowing you to interact with data outside the pipeline context.

## Installing
1. Clone the repository
2. Install and update using `pip`:
  ```bash
  pip install -e .
  ```
3. Install SCICoNE.

## Running
1. Prepare the configuration file according to your analysis
2. Run `snakemake` with:
  ```bash
  snakemake --configfile your_config_file
  ```
3. (Optional) Refer to https://snakemake.readthedocs.io to customise your `snakemake` for your environment

## Contributing

You are very welcome to contribute! You can start with the existing issues or create new issues.
Make sure to follow the CI checks. Use the [pre-commit hook](https://github.com/cbg-ethz/scdna-pipe/blob/master/.pre-commit-config.yaml) defined in the project to meet the code style. If you are adding new functionality, add the corresponding test as well in order to keep the code coverage high.

## License

This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details
