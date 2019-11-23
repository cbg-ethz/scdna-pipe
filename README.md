# Single-cell DNA Analysis <img src="figures/logo.png" align="right" width="160">

[![CircleCI](https://circleci.com/gh/anilbey/dna-pipeline.svg?style=svg&circle-token=7d59442470c38d05f7d1661a97da237d482684ef)](https://circleci.com/gh/anilbey/dna-pipeline)
[![License](http://img.shields.io/:license-Apache%202-green.svg)](http://www.apache.org/licenses/LICENSE-2.0.txt)
[![Python Version](https://img.shields.io/badge/python-3-blue.svg)](https://img.shields.io/badge/python-3-blue.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## About
Reproducible Python data analysis pipeline. Performs single-cell copy number variation calling by learning the underlying tumour evolution history by state-of-the art phylogenetic tree reconstruction method: SCICoNE.The pipeline is built using Python, Conda environment management system and the Snakemake workflow management system. The pipeline starts from the raw sequencing files and a settings file for parameter configurations. After running the data analysis, pipeline produces report and figures to inform the treatment decision of the cancer patient.


## Installing
1. Clone the repository
2. Install and update using `pip`:
  ```bash
  pip install -e .
  ```

## Running
1. Prepare the configuration file according to your analysis
2. Run `snakemake` with:
  ```bash
  snakemake --configfile your_config_file 
  ```
3. (Optional) Refer to https://snakemake.readthedocs.io to customise your `snakemake` for your environment

## Contributing

You are very welcome to contribute! You can start with the existing issues or create new issues.
Make sure to follow the CI checks. Use the [pre-commit hook](https://github.com/anilbey/dna-pipeline/blob/master/.pre-commit-config.yaml) defined in the project to meet the code style. If you are adding new functionality, add the corresponding test as well in order to keep the code coverage high.

## Authors

* **Mustafa Anıl Tuncel**  [:octocat:](https://github.com/anilbey)
* **Pedro Ferreira** [:octocat:](https://github.com/pedrofale)

## License

This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details
