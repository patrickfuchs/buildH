# Test cases

Here are some examples on how to use buildH with POPC Berger. The file `order_parameter_definitions_MODEL_Berger_POPC.def` comes from the [MATCH repository](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs).

All output files (`OUT*`) were obtained by lauching buildH in the following way:

```bash

# On a file with a single POPC (1POPC.pdb)
python ../buildH_calcOP.py 1POPC.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-o OUT.buildH.1POPC.pdb

# On a file with 128 POPC (start_128popc.pdb).
python ../buildH_calcOP.py start_128popc.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-o OUT.buildH.start_128popc.pdb

# On a small trajectory with 25 frames (popc0-25ns_dt1000.xtc).
python ../buildH_calcOP.py start_128popc.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-x popc0-25ns_dt1000.xtc -o OUT.buildH.popc0-25ns_dt1000.xtc

```
