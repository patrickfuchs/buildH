# Test case

Here are some examples on how to use buildH with a simple system made of POPC Berger lipids. The file `order_parameter_definitions_MODEL_Berger_POPC.def` comes from the [MATCH repository](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs).

All output files (`OUT*`) were obtained by lauching buildH in the following way:

```bash

# On a file with a single POPC (1POPC.pdb)
buildH -c 1POPC.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-o OUT.buildH.1POPC.pdb.out

# On a file with 128 POPC (start_128popc.pdb).
buildH -c start_128popc.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-o OUT.buildH.start_128popc.pdb.out

# On a small trajectory with 25 frames (popc0-25ns_dt1000.xtc).
buildH -c start_128popc.pdb -l Berger_POPC \
	-d order_parameter_definitions_MODEL_Berger_POPC.def \
	-t popc0-25ns_dt1000.xtc -o OUT.buildH.popc0-25ns_dt1000.xtc.out

```
