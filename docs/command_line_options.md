# Command line options

We explain in this document the details of the different options on the command line. 

## General usage

When **buildH** is invoked with the flag `-h`, it displays some quite detailed help to the screen:

```
$ buildH
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]]
               -d DEFOP [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END] [-pi PICKLE]

This program builds hydrogens and calculate the order parameters (OP) from a
united-atom trajectory. If -opx is requested, pdb and xtc output files with
hydrogens are created but OP calculation will be slow. If no trajectory output
is requested (no use of flag -opx), it uses a fast procedure to build
hydrogens and calculate the OP.

optional arguments:
  -h, --help            show this help message and exit
  -c COORD, --coord COORD
                        Coordinate file (pdb or gro).
  -t TRAJ, --traj TRAJ  Input trajectory file. Could be in XTC, TRR or DCD format.
  -l LIPID, --lipid LIPID
                        Residue name of lipid to calculate the OP on (e.g.
                        POPC).
  -lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...], --lipid_topology LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]
                        User topology lipid json file(s). Mandatory to build hydrogens.
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
  -b BEGIN, --begin BEGIN
                        The first frame (ps) to read from the trajectory.
  -e END, --end END     The last frame (ps) to read from the trajectory.

The list of supported lipids (-l option) are: CHARMM_POPC, Berger_POPC, Berger_PLA, Berger_POP.
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

`-t TRAJ` or `--traj TRAJ`: `TRAJ` is an input trajectory file in XTC, TRR or DCD format. If not provided, the H reconstruction and order parameter calculation will be done solely on the `COORD` file. In this case, the resulting order parameters will be averaged over the trajectory. More on how this averaging is performed can be found [here](buildh.md#statistics).

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

`-opx OPDBXTC` or `--opdbxtc OPDBXTC`: if you want a trajectory with all hydrogens reconstructed, `OPDBXTC` is a base name. If you supply `-opx traj_wH`, two files will be created: `traj_wH.pdb` and `traj_wH.xtc` both containing all hydrogens. File extension will be automatically added. So far only pdb and xtc are supported. Note that this option is slow.

(**optional flag**)

### Output name for the order parameter

`-o OUT` or `--out OUT`: `OUT` will be the output name for storing the calculated order parameters. If this option is not supplied, the default output name is `OP_buildH.out`.

(**optional flag**)

### Precising beginning and end of trajectory

`-b BEGIN` or `--begin BEGIN`: The first frame (in ps) to read from the trajectory.
`-e END` or `--end END`: The last frame (in ps) to read from the trajectory.
**buildH** checks whether the `BEGIN` and `END` make sense with the supplied trajectory.

(**optional flag**)

### Pickle

*Sorry, this option is not implemented yet*.

`-pi PICKLE` or `--pickle PICKLE`: Output pickle filename. The structure pickled is a dictonnary containing for each Order parameter, the value of each lipid and each frame as a matric

(**optional flag**)
