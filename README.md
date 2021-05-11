# buildH
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4676218.svg)](https://doi.org/10.5281/zenodo.4676218)
[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/patrickfuchs/buildH/master?urlpath=lab)
[![Code CI Status](https://github.com/patrickfuchs/buildH/workflows/GitHub%20CI%20code/badge.svg)](https://github.com/patrickfuchs/buildH/actions?query=workflow%3A%22GitHub+CI+code%22)
[![Doc CI Status](https://github.com/patrickfuchs/buildH/workflows/GitHub%20CI%20doc/badge.svg)](https://github.com/patrickfuchs/buildH/actions?query=workflow%3A%22GitHub+CI+doc%22)
[![Documentation Status](https://readthedocs.org/projects/buildh/badge/?version=latest)](https://buildh.readthedocs.io/en/latest/?badge=latest)
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

> Build hydrogen atoms from a united-atom MD of lipids and calculate the order parameter.

## Features

BuildH can :
  - reconstruct hydrogens from a **united-atom** structure file (PDB, GRO) or a trajectory.
  - calculate the order parameter based on the reconstructed hydrogens
  - write a new structure/trajectory file with the reconstructed hydrogens


BuildH works in two modes :
  1.  A slow mode when an output trajectory (e.g. in xtc format) is requested by
     the user. In this case, the whole trajectory including newly built
     hydrogens are written to this trajectory.
  2. A fast mode without any output trajectory.


## Requirements

Python >= 3.6 is mandatory for running buildH.

buildH is written in Python 3 and needs the modules numpy, pandas, MDAnalysis and Numba.

## Installation

### Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

### Installation within a conda environment

In case you want to install buildH within a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), first create a new conda env:

```
conda create -n env_buildH python pip
```

Then activate your environment:

```
conda activate env_buildH
```

Last, install buildH within that environment using `pip`:

```
pip install buildh
```

*Note*: we recall that once the conda env is activated, when you use `pip` it is the version of `pip` within the conda env, not the one of your Unix system. It allows to embed buildH and all its dependencies within the env without interacting with the Python of the Unix system.

For installing a developement version, see [here](devtools/install_dev.md).

## Launching buildH

Once installed with `pip` as shown above, a simple invocation of `buildH` will launch the program (`$` represents the Unix prompt):

```
$ buildH
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP
              [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END] [-pi PICKLE]
buildH: error: the following arguments are required: -c/--coord, -l/--lipid, -d/--defop
```

The minimal command for launching **buildH** can ressemble this:

```
$ buildH -c start_128popc.pdb -t popc0-25ns_dt1000.xtc -l Berger_POPC -d order_parameter_definitions_MODEL_Berger_POPC.def
```

The different arguments mean the following: `-c start_128popc.pdb` is a pdb file with 128 POPC, `-t popc0-25ns_dt1000.xtc` is a trajectory with 25 frames, `-l Berger_POPC` indicates the united-atom force field and the type of lipid to be analyzed, `-d order_parameter_definitions_MODEL_Berger_POPC.def` indicates what C-H are considered for H building and order parameter calculation. This latter file can be found [here](https://github.com/patrickfuchs/buildH/blob/master/def_files/order_parameter_definitions_MODEL_Berger_POPC.def). The final order parameters averaged over the trajectory will be written to the default output name `OP_buildH.out`

Some other commented examples can be found on the [documentation part](https://github.com/patrickfuchs/buildH#documentation). There is also a Notebook showing a full analysis on a trajectory of 2500 frames (**TODO**: link to be added).

**Important**: sometimes, when performing MD, some molecules are split over periodic boundary conditions (PBC). **buildH** takes as input whole structure (pdb, gro, xtc, etc.). If broken molecules are supplied, it will most likely generate nonsense results. So it is up to the user to take care of making molecules whole before running **buildH** (e.g. by using a tool like *trjconv* in GROMACS).

Invoking **buildH** with the `-h` flag will display some help to the screen and tell the user which lipids are supported by **buildH**.

```
$ buildH -h
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP
              [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END] [-pi PICKLE]
[...]
The list of supported lipids (-l option) are: Berger_POPC, Berger_PLA, Berger_POP, CHARMM_POPC.
```

## Documentation

The full documentation is available on [readthedocs](https://buildh.readthedocs.io/en/latest/index.html).

## Contributors

  - Patrick Fuchs
  - Amélie Bacle
  - Hubert Santuz
  - Pierre Poulain

## Licence

buildH is licensed under the [BSD License](LICENSE).
