# buildH

**buildH** is a software that reads a united-atom (UA) trajectory of lipids, build the hydrogens on it and calculate the order parameter on each C-H bond. **buildH** also allows to output the trajectory with the new reconstructed hydrogens.

In this page you will find general 

## Motivation

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), there was a lack of suitable program for reconstructing hydrogens. In the past, we used to use `g_protonate` in GROMACS 3.*. But this program has been removed in recent versions. Our idea was to build our own implementation in Python upon libraries such MDAnalysis, Numpy and Pandas.
**buildH** is used actively in the recent projects of NMRlipids such as [NMRlipidsIVPEandPG](https://github.com/NMRLipids/NMRlipidsIVPEandPG) or [Databank](https://github.com/NMRLipids/Databank). **buildH** can also be used by anyone willing to analyze the order parameter from a UA trajectory, or if one needs to have explicit hydrogens for some further analyzes.

## Features

**buildH** can:

  - reconstruct hydrogens from a **united-atom** structure file (PDB, GRO) or a trajectory.
  - calculate the order parameter based on the reconstructed hydrogens
  - write a new structure/trajectory file with the reconstructed hydrogens

**buildH** works in two modes:

  1.  A slow mode when an output trajectory is requested by the user. In this case, the whole trajectory including newly built hydrogens is written to this trajectory file. So far, only the xtc format is supported.
  2. A fast mode without any output trajectory.

In both modes, the order parameter is calculated.

It is possible to select only a part of the lipid on which **buildH** will do his job (e.g. the polar head, the sn-1 aliphatic chain, etc.).

## Validation of buildH

**buildH** has been thoroughly validated using a CHARMM36 all-atom trajetory. Everything is detailed in the folder `docs/CHARMM36_POPC_validation`. You can get started with the [README file](CHARMM36_POPC_validation/README.md). There is also a [report](CHARMM36_POPC_validation/report_buildH.pdf) and an [animated gif](CHARMM36_POPC_validation/CHARMM_vs_buildH.gif).

## Simple examples

**TODO**: update this!!!

A couple of test cases on Berger POPC are available in the `docs/Berger_POPC_test_case` folder from the [GitHub repository of the project](https://github.com/patrickfuchs/buildH).

Here are some examples on how to run `buildH`:

- Basic run on a single structure (default name for output OPs will be used):
  ```bash
  buildH -c start_128popc.pdb -l Berger_POPC \
  -d order_parameter_definitions_MODEL_Berger_POPC.def
  ```
- Same but an output file for OPs name is given:
  ```bash
  buildH -c start_128popc.pdb -l Berger_POPC \
  -d order_parameter_definitions_MODEL_Berger_POPC.def \
  -o OP_buildH.out
  ```
- Run `buildH` on a trajectory `traj.xtc`:
  ```bash
  buildH -c start_128popc.pdb -l Berger_POPC \
  -d order_parameter_definitions_MODEL_Berger_POPC.def \
  -t traj.xtc
  ```
- Run `buildH` on a trajectory `traj.xtc` with the trajecory outputs with hydrogens (`traj_with_H.xtc` and `traj_with_H.pdb`), and a default file name for the OP:
  ```bash
  buildH -c start_128popc.pdb -l Berger_POPC \
  -d order_parameter_definitions_MODEL_Berger_POPC.def \
  -t traj.xtc -opx traj_with_H
  ```
  Note that in this last case, the `.def` file **must** contain all possible pairs of C-H to reconstruct. (since the whole trajectory with Hs will be reconstructed).

## Further documentations

Some more detailed are available about:

- [How **buildH** builds hydrogens](algorithms_Hbuilding.md)
- [Command line options](command_line_options.md)
- [The def file format](def_format.md)
- [The lipid json format](json_format.md)
