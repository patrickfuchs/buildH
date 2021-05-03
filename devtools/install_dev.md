In this file is explained how to install buildH for developers.

## Installation (development)

1. Install conda (either with Miniconda or Anaconda, we recommend Miniconda)

2. Clone this GitHub repository:
```
$ git clone https://github.com/patrickfuchs/buildH.git
$ cd buildH
```

3. Create conda environment:
```
$ conda env create -f binder/environment.yml
$ conda activate buildh
```

If needed, update your conda env with
```
$ conda env update -f binder/environment.yml
```

4. Install the dev version of buildH:
```
$ pip install -e .
```
