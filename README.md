# buildH

[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/patrickfuchs/buildH/master?urlpath=lab)

> Build hydrogens from a united-atom MD of lipids and calculate the order parameter. 

## Motivation

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), so far there is a lack of suitable program for reconstructing hydrogens. In the past, we used to use g_protonate in GROMACS 3.* versions. But now, this program has been removed in recent versions. The idea is to build our own using python and a package such as MDAnalysis for reading a trajectory, as well as numpy and possibly others such as pandas.


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

buildH is written in Python 3 and need the following modules :
  - numpy
  - pandas
  - MDAnalysis.

This is automatically taken into account if you follow the procedure below.

## Installation (development)

1. Install conda (either with Miniconda or Anaconda, we recommend Miniconda)

2. Clone this GitHub repository:
```
$ git clone https://github.com/patrickfuchs/buildH.git
$ cd buildH
```

3. Create environment:
```
$ conda env create -f binder/environment.yml
$ conda activate buildh
```

4. Install the dev version of buildH:
```
$ pip install -e .
```
## Usage

```
$ buildH
usage: buildH [-h] [-x XTC] -l LIPID -d DEFOP [-opx OPDBXTC] [-o OUT] [-b BEGIN]
              [-e END] [-pi PICKLE]
              topfile

This program builds hydrogens and calculate the order parameters (OP) from a
united-atom trajectory. If -opx is requested, pdb and xtc output files with
hydrogens are created but OP calculation will be slow. If no trajectory output
is requested (no use of flag -opx), it uses a fast procedure to build
hydrogens and calculate the OP.

positional arguments:
  topfile               Topology file (pdb or gro).

optional arguments:
  -h, --help            show this help message and exit
  -x XTC, --xtc XTC     Input trajectory file in xtc format.
  -l LIPID, --lipid LIPID
                        Residue name of lipid to calculate the OP on (e.g.
                        POPC).
  -d DEFOP, --defop DEFOP
                        Order parameter definition file. Can be found on
                        NMRlipids MATCH repository:https://github.com/NMRLipid
                        s/MATCH/tree/master/scripts/orderParm_defs
  -opx OPDBXTC, --opdbxtc OPDBXTC
                        Base name for trajectory output with hydrogens. File
                        extension will be automatically added. For example
                        -opx trajH will generate trajH.pdb and trajH.xtc. So
                        far only xtc is supported.
  -o OUT, --out OUT     Output base name for storing order parameters.
                        Extention ".out" will be automatically added. Default
                        name is OP_buildH.out.
```

The program needs two mandatory files (present in this repo):
- `dic_lipids.py` (option `-l`, present in this repo) ;
- `order_parameter_definitions_MODEL_Berger_POPC.def` (option `-d`, present in this repo).

## Further documentation

Some more documentation can be found in the directory `docs` :

- Explanation of the different file formats.
- Examples of how to launch buildH.
- Validation of buildH.
- The geometric algorithm on how H are rebuilt.


## Contributors

  - Patrick Fuchs
  - Amélie Bacle
  - Hubert Santuz
  - Pierre Poulain


## Licence

buildH is licensed under the [BSD License](LICENSE).
