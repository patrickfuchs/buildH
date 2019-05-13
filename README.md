# buildH

This repository contains a project for reconstructing hydrogens from a united-atom trajectory.
The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), so far there is a lack of suitable program for reconstructing hydrogens. In the past, we used to use g_protonate in GROMACS 3.* versions. But now, this program has been removed in recent versions. The idea is to build our own using python and a package such as MDAnalysis for reading a trajectory, as well as numpy and possibly others such as pandas.
