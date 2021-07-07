# Usage


## Simple examples

The examples below are based on a simple test case using Berger POPC. The files can be found on github in the directories [docs/Berger_POPC_test_case](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case) and [def_files](https://github.com/patrickfuchs/buildH/tree/master/def_files). You will need 3 files :

- [`start_128popc.pdb`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/start_128popc.pdb): contains 128 POPC.
- [`popc0-25ns_dt1000.xtc`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/popc0-25ns_dt1000.xtc): contains a small trajectory of 25 frames.
- [`Berger_POPC.def`](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def): contains a list of C-H which tells **buildH** what hydrogens to reconstruct, what C-H to calculate the order parameters on.


Here are some examples on how to run **buildH** with these 3 files:

### Basic run on a single structure

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def
```

**buildH** can be used on a single structure (OK not very common for research, but useful for debugging ;-)). The pdb structure is passed with option `-c` (it also works with gro files), the def file with `-d`. The flag `-l` is mandatory, it tells **buildH** what force field and lipid to use: here it is `Berger_POPC`. The order parameters will be written to `OP_buildH.out` which is the default name.

### Same but with a chosen output name

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-o my_OP_buildH.out
```

Here we add a `-o` flag which tells **buildH** to output the results in a file named `my_OP_buildH.out`.

### Run on a trajectory

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc
```

Here the flag `-t` indicates a trajectory. The final order parameters will be averaged over all lipids and all frames for each C-H present in the def file. More can be found on how the averaging is done [here](order_parameter.md). The default name `OP_buildH.out` will be used.

### Same with an output trajectory with reconstructed hydrogens

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_with_H
```

Here we added the flag `-opx` to request a pdb and an xtc file of the system with all the reconstructed hydrogens. Note that the flag takes a base name without extension since it will create a pdb and an xtc, here `popc0-25ns_dt1000_with_H.pdb` and `popc0-25ns_dt1000_with_H.xtc`. The use of this flag `-opx` requires the `.def` file to contain **all possible pairs of C-H** to reconstruct (since the trajectory with all Hs will be reconstructed). Importantly, the newly built hydrogens in the output pdb will be named according to the names written in the def file. See more about this [here](def_format.md). The order parameters will be written in `OP_buildH.out` (default name).

### Get a single pdb file with reconstructed hydrogens

If you do not provide a trajectory with the `-t` flag and you use the `opx` flag, **buildH** will only output a pdb file with hydrogens (no xtc will be produced):

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-opx start_128popc_wH
```

In this case, the file `start_128popc_wH.pdb` with reconstructed hydrogens will be created as well as `OP_buildH.out` with the order parameters.

## Additional details

### Why do I need a def file?

A def file looks like this:

```
gamma1_1 POPC C1  H11
gamma1_2 POPC C1  H12
gamma1_3 POPC C1  H13
[...]
```

Each line corresponds to a given C-H. The 4 columns correspond to the generic name, residue name, carbon name and hydrogen name, respectively, for that C-H.

In **buildH**, the def file has three main purposes:

- Tell what are the C-H we want to consider for H reconstruction and order parameter calculation.
- Give a generic name to each C-H (which will appear in the output) and make the correspondance with the PDB names (e.g. `gamma1_1` stands for the C-H which have `C1` and `H11` atom names in the pdb file.
- If an output file with the newly built hydrogens is requested, their names will follow the 4th column of the def file. For example, the 3 hydrogens reconstructed on atom `C1` will be named `H11`, `H12` and `H13`.


In the following example dealing with a Berger POPC, the order parameters will be calculated on the polar head only (excluding the CH3s of choline):

```
beta1  POPC C5  H51
beta2  POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1   POPC C12 H121
g3_2   POPC C12 H122
g2_1   POPC C13 H131
g1_1   POPC C32 H321
g1_2   POPC C32 H322
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

The list of supported lipids can be requested with `buildH -h`. This command will throw a detailed help to the screen, the list will be indicated at the last line. If you want to analyze a lipid that is not present in **buildH**, you will have to create your own def file as well as a json file which explains to **buildH** how the hydrogens will be reconstructed. This user json file is passed with option `-lt`. Here is more documentation on how to [create your own def file](def_format.md) and how to [create your own json file](json_format.md).

### What about polar hydrogens?

When a lipid contains polar hydrogens, such as the 3 Hs of ethanolamine in PE or the H of the hydroxyl group in cholesterol, these Hs are handled explicitely by the force field. Thus they already exist in the input pdb (and possibly xtc) given as input to **buildH**. In this case, those Hs will be ignored by **buildH** and no order parameter will be calculated on these ones. Usually, these Hs are exchangeable and we do not have experimental order parameters for them. If an output trajectory is requested, **buildH** will just copy the coordinates of these Hs as it does for the heavy atoms.

There is an exception for the force field CHARMM36UA. In this force field, only the apolar Hs of the sn-1 and sn-2 aliphatic tails are in a united-atom representation (starting from the 3rd carbon up to the end of the chain). The other apolar Hs (choline, glycerol, second carbon of sn-1 and sn-2) are explicit. Thus for these latter, **buildH** will ignore them as it does for polar Hs as explained above. Again, if an output pdb (or xtc) is requested, those Hs will be copied to the output pdb and xtc files.

### Mixtures of lipids

If you have a mixture of lipids, you will have to run **buildH** for each lipid separately. If you request an output trajectory, this will have to be done iteratively as well. A guided example on a POPC/POPE mixture can be found in [Notebook03](notebooks/Notebook_03_mixture_POPC_POPE.ipynb). Another one on a POPC/cholesterol mixture can be found in [Notebook05](notebooks/Notebook_05_mixture_POPC_cholesterol.ipynb).

### Periodic boundary conditions

Sometimes, when performing MD, some molecules are split over periodic boundary conditions (PBC). **buildH** takes as input whole structures (pdb, gro, xtc, etc.). If broken molecules are supplied, it will most likely generate nonsense results. So it is up to the user to take care of making molecules whole before running **buildH** (e.g. by using a tool like [trjconv](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html) in GROMACS with flag `-pbc mol`).


### BuildH as a module

buildH is intended to be used mainly in the Unix command line. It is also possible to use it as a module but to a lesser extent.
The features available are minimal: you can just call the main function (`buildh.launch()`) and result files are still written.

It's not a proper API but more a way to call buildH inside larger analysis python scripts.

A guided example can be found on [Notebook04](notebooks/Notebook_04_library.ipynb).

