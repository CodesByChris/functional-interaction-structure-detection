[![GitHub License](https://img.shields.io/github/license/CodesByChris/functional-interaction-structure-detection?label=License)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10601461.svg)](https://doi.org/10.5281/zenodo.10601461)


# Detecting the Functional Interaction Structure of Software Development Teams

Welcome to the repository containing the reproducibility scripts for the paper
*"Detecting the Functional Interaction Structure of Software Development Teams"*.
The repository provides (i) a reproducibility script for the paper and (ii) tools to detect a team's functional interaction structure from interaction networks.


## Usage

To run the reproducibility script, follow these steps.
They are written for a system with Ubuntu Linux but work accordingly on other OSes.

1. Install R and packages:
    ```bash
    sudo apt update
    sudo apt install r-base r-cran-devtools r-cran-irkernel jupyter-notebook
    Rscript -e "install.packages(readLines('packages_cran.txt'))"
    Rscript -e "devtools::install_github(readLines('packages_github.txt'))"
    ```
2. Open `Reproducibility.ipynb`:
    ```bash
    jupyter notebook Reproducibility.ipynb
    ```
3. Run all cells in the notebook.


## Contributors

The code in this repository has been developed by

- Christian Zingg
- Christoph Gote

at the Chair of Systems Design, ETH Zurich.


## Copyright

This repository is released under the GNU Affero General Public License v3.0

Copyright 2020-2024, ETH Zurich.
