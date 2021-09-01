# Installation

## Requirements and compatibility

buildH requires at least Python 3.6 and needs the following modules:
  - numpy
  - pandas
  - numba
  - MDAnalysis (with support of 2.0)

All the instructions below have been tested on Unix like platforms (e.g. Ubuntu), which we recommend for running **buildH**. We do not provide support for other platforms, but since **buildH** has been written in pure Python, it should work there provided its dependencies are supported. 

## Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

Note that this way of proceeding will install **buildH** and its dependencies within the python of your Unix system, which may lead to conflicts of version if you have other scientific packages installed. To avoid this you may want to create a specific conda or virtual environment for **buildH** (see below).

## Installation within a conda environment

**buildH** is also available through the [Bioconda](https://anaconda.org/bioconda/buildh) channel.

We recommend to install **buildH** within a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). First create a new conda env:

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

## Building from source

We recommend to use a specific environment, either by using [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands). All packages will be installed within that environment which avoid conflicts of version.

If you still do not want to create a specific environment for **buildH**, you can skip the first section `Create an environment` below.

In any case, the python version should be >= 3.6.

### Create an environment

First, create a new environment (we call it `env4buildH`):

- If you chose conda: `conda create -n env4buildH python=3.8`
- If you chose venv: `python3 -m venv /path/to/env4buildH`

Activate your environment:

- If you chose conda: `conda activate env4buildH`
- If you chose venv: `source /path/to/env4buildH/bin/activate`

### Install **buildH** from source

Clone the **buildH** repository:

```
git clone https://github.com/patrickfuchs/buildH.git
cd buildH
```

Install with `pip` the packages required by **buildH**, namely numpy, pandas, MDAnalysis and Numba, which are all specified in the file [requirements.txt](https://github.com/patrickfuchs/buildH/blob/master/requirements.txt):

```
pip install -r requirements.txt
```

Install **buildH** from source with `pip`:

```
pip install -e .
```

## For developers

For installing a development version from source with the full environment (allowing building the doc, launching tests, etc.), see [here](https://github.com/patrickfuchs/buildH/tree/master/devtools/install_dev.md).

## Testing

The tests rely on the `pytest` package. You need first to install a development version (see above). Once done, you can run the tests with just:

```
cd buildh # if it's not already the case
pytest
```

All tests should pass. If anything fails, please [open an issue](https://github.com/patrickfuchs/buildH/issues).
