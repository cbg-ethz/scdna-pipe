# Single-cell DNA Analysis <img src="figures/logo.png" align="right" width="160">

[![CircleCI](https://circleci.com/gh/anilbey/dna-pipeline.svg?style=svg&circle-token=7d59442470c38d05f7d1661a97da237d482684ef)](https://circleci.com/gh/anilbey/dna-pipeline)
[![License](http://img.shields.io/:license-Apache%202-green.svg)](http://www.apache.org/licenses/LICENSE-2.0.txt)
[![Python Version](https://img.shields.io/badge/python-3-blue.svg)](https://img.shields.io/badge/python-3-blue.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## About
Reproducible Python data analysis pipeline. Performs single-cell copy number variation calling by learning the underlying tumour evolution history by state-of-the art phylogenetic tree reconstruction method: SCICoNE.The pipeline is built using Python, Conda environment management system and the Snakemake workflow management system. The pipeline starts from the raw sequencing files and a settings file for parameter configurations. After running the data analysis, pipeline produces report and figures to inform the treatment decision of the cancer patient.


## Installing
----------
1. Clone the repository
2. Install and update using `pip`:
  ```bash
  pip install -e .
  ```

## Authors

* **Mustafa Anıl Tuncel**  [:octocat:](https://github.com/anilbey)
* **Pedro Ferreira** [:octocat:](https://github.com/pedrofale)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details
