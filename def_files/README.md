This directory contains def files which tell **buildH** which C-H to analyze and give them a generic name for a more readable output. Initially, the use of def files comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project, some def files for all-atom force fields can be found [here](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs).

The global format is the following:

```
name            residue   carbon atom    hydrogen atom 
 
gamma1_1         POPC        C13             H13A 
gamma1_2         POPC        C13             H13B 
gamma1_3         POPC        C13             H13C 
. 
. 
.
```

Each column has to be separated by any combination of whitespaces (at least one).

Last, in the [NMRlipids](https://nmrlipids.blogspot.com/) project, the def files are generally named like this: `order_parameter_definitions_MODEL_CHARMM36_POPC.def`. In **buildH** we use merely `CHARMM36_POPC.def`, or more generally `Forcefield_Lipid.def` where `Lipid` is a residue name as found in the pdb or gro file.


More on def files can be found in the [documentation](https://buildh.readthedocs.io/en/latest/def_format.html).
