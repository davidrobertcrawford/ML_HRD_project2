# Reproducing experiments from Davies, et al. (_Nature Medince_, 2017)

TODO: Fill in details of planned experiments to be reproduced.

Additional details of the experiments, data, and conclusions can be found in the [`report/`](report/).

## Setup

The source code is written in Python 3. We use `snakemake` to manage the workflow. We suggest using Conda to install dependencies, which you can do directly from the provided [`environment.yml`](environment.yml) file.

    conda env create -f environment.yml
    source activate reproducing-davies2017-env

We have included two git modules that you will need to download before you can run `snakemake all`

    git submodule init
    git submodule update

## Usage

To download the input data files and reproduce the experiments and..., simply run:

    snakemake all

TODO: insert some expectation of runtime.

### Configuration

Additional configuration options are detailed at the beginning of the [`Snakefile`](Snakefile).
There are also various command-line options for the provided scripts in [`src/`](src).

## References
1. Davies, et al. (2017) "HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures." _Nature Medicine_ **23**, pages 517-525. [doi: 10.1038/nm.4292](https://doi.org/10.1038/nm.4292)
