# buildH

**buildH** is a software that reads a united-atom (UA) trajectory of lipids, build the hydrogens on it and calculate the order parameter on each C-H bond. **buildH** also allows to output the trajectory with the new reconstructed hydrogens.

**buildH** works in two modes:

  1.  A slow mode when an output trajectory is requested by the user. In this case, the whole trajectory including newly built hydrogens is written to this trajectory file. So far, only the xtc format is supported.
  2. A fast mode without any output trajectory.

In both modes, the order parameter is calculated.

It is possible to select only a part of the lipid on which **buildH** will do his job (e.g. the polar head, the sn-1 aliphatic chain, etc.) thanks to the [def file](def_format.md).

In this page you will find general information about **buildH**.

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

## Motivation

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), there was a lack of suitable program for reconstructing hydrogens. In the past, we used to use `g_protonate` in GROMACS 3.*. But this program has been removed in recent versions. Our idea was to build our own implementation in Python using libraries such MDAnalysis, Numpy and Pandas.
**buildH** is used actively in the recent projects of NMRlipids such as [NMRlipidsIVPEandPG](https://github.com/NMRLipids/NMRlipidsIVPEandPG) or [Databank](https://github.com/NMRLipids/Databank). **buildH** can also be used by anyone willing to analyze the order parameter from a UA trajectory, or if one needs to have explicit hydrogens for some further analyzes.

## Validation of buildH

**buildH** has been thoroughly validated using a CHARMM36 all-atom trajetory. Everything is detailed in the folder `docs/CHARMM36_POPC_validation`. You can get started with the [README file](CHARMM36_POPC_validation/README.md). There is also a [report](CHARMM36_POPC_validation/report_buildH.pdf) and an [animated gif](CHARMM36_POPC_validation/CHARMM_vs_buildH.gif).

## Simple examples

No more BLABLA, please show me how to run **buildH**! OK, the examples below are based on a simple test case using Berger POPC. The files can be found on github in the directory [docs/Berger_POPC_test_case](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case). You will need 3 files :

- `start_128popc.pdb`: contains 128 POPC.
- `popc0-25ns_dt1000.xtc`: contains a small trajectory of 25 frames.
- `order_parameter_definitions_MODEL_Berger_POPC.def`: contains a list of C-H which tells **buildH** what hydrogens to reconstruct and what C-H to calculate the order parameter on.


Here are some examples on how to run `buildH`:

### Basic run on a single structure

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d order_parameter_definitions_MODEL_Berger_POPC.def
```

**buildH** can be used on a single structure (OK not very common for reasearch, but useful for debugging ;-)). The pdb structure is passed with option `-c` (it also works with gro files), the def file with `-d`. The flag `-l` is mandatory, it tells buildH the force field and the lipid: here it is `Berger_POPC`. The order parameters will be written in `OP_buildH.out` which is the default name.

### Same but with a chosen output name

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d order_parameter_definitions_MODEL_Berger_POPC.def \
-o my_OP_buildH.out
```

Here we add a `-o` flag which tells **buildH** to output the results in a file with name `my_OP_buildH.out`.

### Run on a trajectory

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d order_parameter_definitions_MODEL_Berger_POPC.def \
-t popc0-25ns_dt1000.xtc
```

Here the flag `-t` indicates a trajectory. The final order parameters will be averaged over all lipids and all frames for each C-H present in the def file. More can be found on how the averaging is done [below](buildH.md#statistics). The default name `OP_buildH.out` will be used.

### Same with an output trajectory with reconstructed hydrogens

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d order_parameter_definitions_MODEL_Berger_POPC.def \
-t popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_with_H
```

Here we had the flag `-opx` to request a pdb and an xtc file of the system with all the reconstructed hydrogens. Note that the flag takes a base name without extension since it will create a pdb and an xtc, here `popc0-25ns_dt1000_with_H.pdb` and `popc0-25ns_dt1000_with_H.xtc`. The use of this flag `-opx` requires the `.def` file to contain **all** possible pairs of C-H to reconstruct (since the trajectory with all Hs will be reconstructed).

## Additional features

### Supported lipids

The list of supported lipids can be requested with `buildH -h`, it is indicated at the last line. If you want to analyze a lipid that is not present in **buildH**, you will have to create your own def file as well as a json file which explains to **buildH** how the hydrogens will be reconstructed. This user json file is passed with option `-lt`. More documentation on how to [create your own def file](def_format.md) and how to [create your own json file](json_format.md).

### Mixtures of lipids

If you have a mixture of lipids, you will have to run **buildH** for each lipid separately.

### Statistics

The order parameter output of buildH (default name `OP_buildH.out`) looks like this:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
gamma1_1             POPC    C1    H11    0.01304  0.12090  0.01069
gamma1_2             POPC    C1    H12    0.00666  0.09279  0.00820
gamma1_3             POPC    C1    H13   -0.01531  0.09141  0.00808
[...]
```

Each line corresponds to a given CH. The 4 first columns contain the generic name, residue name, carbon and hydrogen names respectively. The other column contains different statistics:

- `OP_mean` is the order parameter averaged over all lipids and all frames of the trajectory.
- `OP_stddev` is the standard deviation of the order parameter; first we average each C-H over the whole trajectory, then we calculate the standard deviation over all residues: 
$$ OP\_stddev(CH_j) = \frac{1}{nres} \sum_{i=1}^{i=nres} 
\left[ \frac{1}{nframes} \sum_{t=0}^{t=nframes} OP(CH_j)(i)(t) \right]$$
where $CH_i$ is the $j^{th}$ C-H, $nframes$ is the total number of frames, $nres$ is the total number of residues (i.e. lipids).
- `OP_stem` is the standard error of the mean, averaged in the same spirit: 
$$OP\_stem(CH_j) = \frac{OP\_stddev(CH_j)}{\sqrt{nres}}$$

## Further documentations

Some more detailed are available about:

- [More on command line options](command_line_options.md)
- [How **buildH** builds hydrogens](algorithms_Hbuilding.md)
- [The def file format](def_format.md)
- [The lipid json format](json_format.md)
