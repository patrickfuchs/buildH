# Command line options

We explain in this document the details of the different options on the command line. 

## General usage

When **buildH** is invoked with the flag `-h`, it displays some quite detailed help to the screen:

```
usage: buildH [-h] [-v] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END]
              [-igch3]

This program builds hydrogens and calculates the order parameters (OP) from a united-atom trajectory of lipids. If -opx is requested, pdb and xtc
output files with hydrogens are created but OP calculation will be slow. If no trajectory output is requested (no use of flag -opx), it uses a fast
procedure to build hydrogens and calculate the OP.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -c COORD, --coord COORD
                        Coordinate file (pdb or gro format).
  -t TRAJ, --traj TRAJ  Input trajectory file. Could be in XTC, TRR or DCD format.
  -l LIPID, --lipid LIPID
                        Combinaison of ForceField name and residue name for the lipid to calculate the OP on (e.g. Berger_POPC).It must match with
                        the internal topology files or the one(s) supplied.A list of supported terms is printed when calling the help.
  -lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...], --lipid_topology LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]
                        User topology lipid json file(s).
  -d DEFOP, --defop DEFOP
                        Order parameter definition file. Can be found on https://github.com/patrickfuchs/buildH/tree/master/def_files.
  -opx OPDBXTC, --opdbxtc OPDBXTC
                        Base name for trajectory output with hydrogens. File extension will be automatically added. For example -opx trajH will
                        generate trajH.pdb and trajH.xtc. So far only xtc is supported.
  -o OUT, --out OUT     Output file name for storing order parameters. Default name is OP_buildH.out.
  -b BEGIN, --begin BEGIN
                        The first frame (ps) to read from the trajectory.
  -e END, --end END     The last frame (ps) to read from the trajectory.
  -igch3, --ignore-CH3s
                        Ignore CH3s groups for the construction of hydrogens and the calculation of the OP.

The list of supported lipids (-l option) are: Berger_CHOL, Berger_DOPC, Berger_DPPC, Berger_POPC, Berger_PLA, Berger_POP, Berger_POPE, Berger_POPS,
CHARMM36UA_DPPC, CHARMM36UA_DPUC, CHARMM36_POPC, GROMOS53A6L_DPPC, GROMOSCKP_POPC, GROMOSCKP_POPS. More documentation can be found at
https://buildh.readthedocs.io.
```

Importantly, the `-h` option also displays the list of supported lipids at the end. Note that they are always written using the naming convention `ForceField_Lipid`.

## Description of options

Each option is explained below. Some of them are mandatory.

### Help

`-h` or `--help`: display the help message shown above and exit.

(**optional flag**)

### Coordinates

`-c COORD` or `--coord COORD`: `COORD` is the main coordinate file in pdb or gro format describing the system. This is strictly required.

(**mandatory flag**)

### Trajectory

`-t TRAJ` or `--traj TRAJ`: `TRAJ` is an input trajectory file in XTC, TRR or DCD format. If not provided, the H reconstruction and order parameter calculation will be done solely on the `COORD` file. If a trajecotry is provided with `-t` flag, the resulting order parameters will be averaged over that trajectory. More on how this averaging is performed can be found [here](https://buildh.readthedocs.io/en/latest/buildh.html#order-parameters-and-statistics).

(**optional flag**)

### Lipid requested

`-l LIPID` or `--lipid LIPID`: `LIPID` is the name of the lipid to calculate the OP on. It must follow the naming convention `ForceField_Lipid`, for example `Berger_POPC`. The list of supported lipids can be queried with the `-h` flag.

(**mandatory flag**)

### User lipid (json file)

`-lt LIPID_TOPOLOGY` or `--lipid_topology LIPID_TOPOLOGY`: `LIPID_TOPOLOGY` is a user supplied topology lipid json file(s). When you want to analyze a lipid not present in **buildH** you can [build your own json file](json_format.md) and supply it with this option. Again, it has to follow the naming convention `ForceField_Lipid.json`. For example, if you build your own json file for butane with the Berger force field, you can use `-lt Berger_BUTA.json`; in this case, you will have to use also `-l Berger_BUTA` flag. 

(**optional flag**)

### Definition file

`-d DEFOP` or `--defop DEFOP`: `DEFOP` is the order parameter definition file. It tells **buildH** what C-H will be considered for H reconstruction and order parameter calculation. You can find some [def files](https://github.com/patrickfuchs/buildH/tree/master/def_files) on the **buildH** repository for the supported lipids. You can also create your [own def file](def_format.md).

(**mandatory flag**)

### Output trajectory

`-opx OPDBXTC` or `--opdbxtc OPDBXTC`: if you want a trajectory with all hydrogens reconstructed, `OPDBXTC` is a base name. If you supply `-opx traj_wH`, two files will be created: `traj_wH.pdb` and `traj_wH.xtc` both containing all hydrogens. File extension will be automatically added. If no trajectory is supplied with the option `-t`, **buildH** will only create a pdb but not an xtc. So far only pdb and xtc are supported. Note that this option is slow.

(**optional flag**)

### Output name for the order parameter

`-o OUT` or `--out OUT`: `OUT` will be the output name for storing the calculated order parameters. If this option is not supplied, the default output name is `OP_buildH.out`.

(**optional flag**)

### Precising beginning and end of trajectory

`-b BEGIN` or `--begin BEGIN`: The first frame (in ps) to read from the trajectory.

`-e END` or `--end END`: The last frame (in ps) to read from the trajectory.

**buildH** checks whether the `BEGIN` and `END` make sense with the supplied trajectory.

(**optional flag**)

### Ignoring CH3 in the output

`-igch3` or `--ignore-CH3s`: when this flag is on, the OP of each C-H belonging to any CH3 will not be computed and not be written in the ouput file even if they are present in the def file. Individual reconstruted C-H of methyl groups are not very interesting for OP calculation since we cannot precisely know where they are because of the methyl rotation. With this option, one can avoid their evaluation without having to remove these C-H in the def file (recall, the [def files from the **buildH** website](https://github.com/patrickfuchs/buildH/tree/master/def_files) contain all possible C-Hs). However, note that this `-igch3` option is not usable with the `-opx` option (which outputs the structure pdb and xtc files with hydrogens) since this latter needs all hydrogen atoms to reconstruct.
