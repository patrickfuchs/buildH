# buildH

**buildH** is a software that reads a united-atom (UA) trajectory of lipids, builds the hydrogens on it and calculates the order parameter on each C-H bond. **buildH** also allows to output the trajectory with the new reconstructed hydrogens.

**buildH** works in two modes:

  1.  A slow mode when an output trajectory is requested by the user. In this case, the whole trajectory including newly built hydrogens is written to this trajectory file. If other molecules are present (e.g. water, ions, etc.), they will just be copied to the output trajectory with the same coordinates. So far, only the xtc format is supported.
  2. A fast mode without any output trajectory.

In both modes, the order parameters are calculated.

It is possible to select only a part of the lipid on which **buildH** will do his job (e.g. the polar head, the sn-1 aliphatic chain, etc.) thanks to the [def file](def_format.md).

In this page you will find general information about **buildH**.

**buildH** is hosted on [github](https://github.com/patrickfuchs/buildH).

## Installation

### Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

### Installation within a conda environment

**buildH** is also available through the [Bioconda](https://bioconda.github.io/) channel.

In case you want to install **buildH** within a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), first create a new conda env:

```
conda create -n env_buildH "python>=3.6"
```

Then activate your environment:

```
conda activate env_buildH
```

Last, install **buildH** within that environment:

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda install buildh
```

For installing a developement version, see [here](https://github.com/patrickfuchs/buildH/tree/master/devtools/install_dev.md).

## Motivation

The initial motivation comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project. As stated in this [post](https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html), there was a lack of suitable program for reconstructing hydrogens. In the past, we used to use `g_protonate` in GROMACS 3.*. But this program has been removed in recent versions. Our idea was to build our own implementation in Python using libraries such MDAnalysis, Numpy and Pandas.
**buildH** is used actively in the recent projects of NMRlipids such as [NMRlipidsIVPEandPG](https://github.com/NMRLipids/NMRlipidsIVPEandPG) or [Databank](https://github.com/NMRLipids/Databank). **buildH** can also be used by anyone willing to analyze the order parameter from a UA trajectory, or if one needs to have explicit hydrogens for some further analyzes.

## Simple examples

No more BLABLA, please show me how to run **buildH**! OK, the examples below are based on a simple test case using Berger POPC. The files can be found on github in the directory [docs/Berger_POPC_test_case](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case). You will need 3 files :

- [`start_128popc.pdb`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/start_128popc.pdb): contains 128 POPC.
- [`popc0-25ns_dt1000.xtc`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/popc0-25ns_dt1000.xtc): contains a small trajectory of 25 frames.
- [`Berger_POPC.def`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/Berger_POPC.def): contains a list of C-H which tells **buildH** what hydrogens to reconstruct and what C-H to calculate the order parameter on.


Here are some examples on how to run **buildH** with these 3 files:

### Basic run on a single structure

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d Berger_POPC.def
```

**buildH** can be used on a single structure (OK not very common for research, but useful for debugging ;-)). The pdb structure is passed with option `-c` (it also works with gro files), the def file with `-d`. The flag `-l` is mandatory, it tells **buildH** what force field and lipid to use: here it is `Berger_POPC`. The order parameters will be written to `OP_buildH.out` which is the default name.

### Same but with a chosen output name

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d Berger_POPC.def \
-o my_OP_buildH.out
```

Here we add a `-o` flag which tells **buildH** to output the results in a file named `my_OP_buildH.out`.

### Run on a trajectory

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc
```

Here the flag `-t` indicates a trajectory. The final order parameters will be averaged over all lipids and all frames for each C-H present in the def file. More can be found on how the averaging is done [below](https://github.com/patrickfuchs/buildH/blob/master/docs/buildh.md#statistics). The default name `OP_buildH.out` will be used.

### Same with an output trajectory with reconstructed hydrogens

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_with_H
```

Here we added the flag `-opx` to request a pdb and an xtc file of the system with all the reconstructed hydrogens. Note that the flag takes a base name without extension since it will create a pdb and an xtc, here `popc0-25ns_dt1000_with_H.pdb` and `popc0-25ns_dt1000_with_H.xtc`. The use of this flag `-opx` requires the `.def` file to contain **all possible pairs of C-H** to reconstruct (since the trajectory with all Hs will be reconstructed). The order parameters will be written in `OP_buildH.out` (default name).

### Get a single pdb file with reconstructed hydrogens

If you do not provide a trajectory with the `-t` flag and you use the `opx` flag, **buildH** will only output a pdb file with hydrogens (no xtc will be produced):

```bash
buildH -c start_128popc.pdb -l Berger_POPC \
-d Berger_POPC.def \
-opx start_128popc_wH
```

In this case, the file `start_128popc_wH.pdb` with reconstructed hydrogens will be created as well as `OP_buildH.out` with the order parameters.

## Additional features

### Why do I need a def file?

A def file looks like this:

```
gamma1_1 POPC C1  H11
gamma1_2 POPC C1  H12
gamma1_3 POPC C1  H13
[...]
```

Each line corresponds to a given C-H. The 4 columns correspond to the generic name, residue name, carbon name and hydrogen name, respectively, for that C-H.

In **buildH**, the def file has two main purposes:

- Tell what are the C-H we want to consider for H reconstruction and order parameter calculation.
- Give a generic name to each C-H (which will appear in the output) and make the correspondance with the PDB names (e.g. `gamma1_1` stands for the C-H which have `C1` and `H11` atom names in the pdb file.


For example, if you want to calculate the order parameters only on the polar head (excluding the CH3s of choline) of a Berger POPC, you can use:

```
beta1 POPC C5  H51
beta2 POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1 POPC C12 H121
g3_2 POPC C12 H122
g2_1 POPC C13 H131
g1_1 POPC C32 H321
g1_2 POPC C32 H322
```

Using the [Berger POPC trajectory](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case) of 25 frames, the output `OP_buildH.out` will contain the order parameters of the C-H specified in the def file:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
beta1                POPC    C5    H51    0.04934  0.11999  0.01061
beta2                POPC    C5    H52    0.07162  0.12108  0.01070
alpha1               POPC    C6    H61    0.11839  0.15261  0.01349
alpha2               POPC    C6    H62    0.13903  0.19003  0.01680
g3_1                 POPC    C12   H121  -0.28674  0.09135  0.00807
g3_2                 POPC    C12   H122  -0.16195  0.14832  0.01311
g2_1                 POPC    C13   H131  -0.15159  0.14511  0.01283
g1_1                 POPC    C32   H321   0.21133  0.22491  0.01988
g1_2                 POPC    C32   H322   0.09638  0.16189  0.01431
```

The def files of the lipids supported by **buildH** can be found [here](https://github.com/patrickfuchs/buildH/tree/master/def_files).

More on def files and creating your own ones can be found [here](def_format.md).

### Supported lipids

The list of supported lipids by **buildH** can be requested with `buildH -h`. This command will throw a detailed help to the screen, the list will be indicated at the last line. If you want to analyze a lipid that is not present in **buildH**, you will have to create your own def file as well as a json file which explains to **buildH** how the hydrogens will be reconstructed. This user json file is passed with option `-lt`. Here is more documentation on how to [create your own def file](def_format.md) and how to [create your own json file](json_format.md).

### Mixtures of lipids

If you have a mixture of lipids, you will have to run **buildH** for each lipid separately. If you request an output trajectory, this will have to be done iteratively as well. An example is shown in **TODO**.

### Order parameters and statistics

The order parameter of bond $CH_j$ is calculated using the standard formula:

$$S_{CH_j} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the $CH_j$ bond and the normal to the membrane (usually the *z* axis), <...> means averaging over molecules and frames. $S_{CH}$ can be measured by NMR which is useful to validate simulation results, as largely described in the [NMRlipids project](http://nmrlipids.blogspot.com).

The order parameter output of buildH (default name `OP_buildH.out`) looks like this:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
gamma1_1             POPC    C1    H11    0.01304  0.12090  0.01069
gamma1_2             POPC    C1    H12    0.00666  0.09279  0.00820
gamma1_3             POPC    C1    H13   -0.01531  0.09141  0.00808
[...]
```

Each line corresponds to a given CH. The 4 first columns contain the generic name, residue name, carbon and hydrogen names respectively. The other columns contains different statistics on order parameters (OP):

- `OP_mean` is the OP of bond $CH_j$ averaged over all lipids and all frames of the trajectory, we shall write it $\overline{S_{CH_j}}$.
- `OP_stddev` is the standard deviation of the OP, we shall write it $\sigma(S_{CH_j})$; first we average each OP of bond $CH_j$ (e.g. the CH of beta1) of residue $i$ (i.e. lipid $i$) over the whole trajectory:

$$ \overline{S_{CH_j}(i)} = \frac{1}{nframes} \sum_{t=0}^{t=nframes} S_{CH_j}(i)(t) $$

where $nframes$ is the total number of frames, then we calculate the standard deviation of those means over all residues:

$$ \sigma(S_{CH_j}) =
\sqrt{
\frac{1}{nres} \sum_{i=1}^{i=nres} (\overline{S_{CH_j}(i)} - \overline{S_{CH_j}})^2
}$$

where $nres$ is the total number of residues (i.e. lipids).
- `OP_stem` is the standard error of the mean averaged in the same spirit, let's call it $err(S_{CH_j})$:

$$err(S_{CH_j}) = \frac{\sigma(S_{CH_j})}{\sqrt{nres}}$$

### Periodic boundary conditions

Sometimes, when performing MD, some molecules are split over periodic boundary conditions (PBC). **buildH** takes as input whole structures (pdb, gro, xtc, etc.). If broken molecules are supplied, it will most likely generate nonsense results. So it is up to the user to take care of making molecules whole before running **buildH** (e.g. by using a tool like [trjconv](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html) in GROMACS with flag `-pbc mol`).
