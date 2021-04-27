# Command line options

We explain in this document the details of the different options on the command line. 

## Usage

First, when **buildH** is invoked with the flag `-h`, it gives quite some details:

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
  -pi PICKLE, --pickle PICKLE
                        Output pickle filename. The structure pickled is a dictonnary containing for each Order parameter,
                        the value of each lipid and each frame as a matric

The list of supported lipids (-l option) are: CHARMM_POPC, Berger_POPC, Berger_PLA, Berger_POP.
```

## Description of options

Each option is explained below. Some of them are mandatory.

### Help

-h, --help            show this help message and exit

### Coordinates

-c COORD, --coord COORD
                        Coordinate file (pdb or gro).
(**mandatory option**)

### Trajectory 

-t TRAJ, --traj TRAJ  Input trajectory file. Could be in XTC, TRR or DCD format.

(**optional**)

### Lipid requested

-l LIPID, --lipid LIPID
                        Residue name of lipid to calculate the OP on (e.g.
                        POPC).
(**mandatory option**)

### User lipid (json file)

  -lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...], --lipid_topology LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]
                        User topology lipid json file(s). Mandatory to build hydrogens.
(**optional**)

### Definition file

  -d DEFOP, --defop DEFOP
                        Order parameter definition file. Can be found on
                        NMRlipids MATCH repository:https://github.com/NMRLipid
                        s/MATCH/tree/master/scripts/orderParm_defs
(**mandatory option**)

### Output trajectory
  -opx OPDBXTC, --opdbxtc OPDBXTC
                        Base name for trajectory output with hydrogens. File
                        extension will be automatically added. For example
                        -opx trajH will generate trajH.pdb and trajH.xtc. So
                        far only xtc is supported.
(**optional**)

### Output name for the order parameter

  -o OUT, --out OUT     Output base name for storing order parameters.
                        Extention ".out" will be automatically added. Default
                        name is OP_buildH.out.
(**optional**)

### Precising beginning and end of trajectory

  -b BEGIN, --begin BEGIN
                        The first frame (ps) to read from the trajectory.
  -e END, --end END     The last frame (ps) to read from the trajectory.
(**optional**)

### Pickle

*This option is not implemented yet*.

  -pi PICKLE, --pickle PICKLE
                        Output pickle filename. The structure pickled is a dictonnary containing for each Order parameter,
                        the value of each lipid and each frame as a matric
(**optional**)
