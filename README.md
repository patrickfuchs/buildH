# buildH

This repository contains a project for reconstructing hydrogens from a united-atom trajectory and calculate the order parameter.

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), so far there is a lack of suitable program for reconstructing hydrogens. In the past, we used to use g_protonate in GROMACS 3.* versions. But now, this program has been removed in recent versions. The idea is to build our own using python and a package such as MDAnalysis for reading a trajectory, as well as numpy and possibly others such as pandas.

It is writte in Python 3.

Prerequisites : argparse, io, numpy, pandas, MDAnalysis.

Easiest way to install all of these is via conda :

```
conda config --add channels conda-forge
conda install numpy pandas mdanalysis
```

Launching the program without arguments shows its usage:

```
$ python ./add_hydrogens.py
usage: add_hydrogens.py [-h] [--xtc XTC] [--pdbout PDBOUT] [--xtcout XTCOUT]
                        topfile
add_hydrogens.py: error: the following arguments are required: topfile
```

Examples of ways of lauching the program:

```
python ./add_hydrogens.py popc.pdb

python ./add_hydrogens.py popc.pdb --xtc traj.xtc

python ./add_hydrogens.py popc.pdb --xtc traj.xtc --pdbout popc_with_H.pdb

python ./add_hydrogens.py popc.pdb --xtc traj.xtc --pdbout popc_with_H.pdb --xtcout traj_with_H.xtc
```
