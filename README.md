# buildH

This repository contains a project for reconstructing hydrogens from a united-atom trajectory and calculate the order parameter.

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), so far there is a lack of suitable program for reconstructing hydrogens. In the past, we used to use g_protonate in GROMACS 3.* versions. But now, this program has been removed in recent versions. The idea is to build our own using python and a package such as MDAnalysis for reading a trajectory, as well as numpy and possibly others such as pandas.

It is writte in Python 3.

Prerequisites : 
- internal: argparse, pickle
- external: numpy, pandas, MDAnalysis.

Easiest way to install all of these is via conda :

```
conda config --add channels conda-forge
conda install numpy pandas mdanalysis
```

Launching the program without arguments shows its usage:

```
python ./buildH_calcOP.py
usage: buildH_calcOP.py [-h] [-x XTC] [-opx OPDBXTC] [-o OUT] topfile
buildH_calcOP.py: error: the following arguments are required: topfile
```

The program needs two files (present in this repo):
- `dic_lipids.py`
- `order_parameter_definitions_MODEL_Berger_POPC.def`

Examples of ways of launching the program:

```
python ./buildH_calcOP.py popc.pdb

python ./buildH_calcOP.py popc.pdb -x traj.xtc

python ./buildH_calcOP.py popc.pdb -x traj.xtc -opx popc_with_H
```
