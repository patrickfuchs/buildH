# Installation

## Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

## Installation within a conda environment

**buildH** is also available through the [Bioconda](https://bioconda.github.io/) channel.

In case you want to install **buildH** within a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), first create a new conda env:

```
conda create -n buildH "python>=3.6"
```

Then activate your environment:

```
conda activate buildH
```

Last, install **buildH** within that environment:

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda install buildh
```

For installing a developement version, see [here](https://github.com/patrickfuchs/buildH/tree/master/devtools/install_dev.md).
