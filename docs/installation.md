# Installation

## Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

## Installation within a conda environment

**buildH** is also available through the [Bioconda](https://anaconda.org/bioconda/buildh) channel.

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

## For developpers

For installing a developement version, see [here](https://github.com/patrickfuchs/buildH/tree/master/devtools/install_dev.md).

## Testing

The tests rely on the `pytest` package. You need first to install a development version (see above). Once it's done, you can run the tests with just:

```
cd buildh # if it's not already the case
pytest
```

All tests should pass. If anything fails, please [open an issue](https://github.com/patrickfuchs/buildH/issues).
