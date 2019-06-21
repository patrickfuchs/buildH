# Test cases

Here some examples on how to use buildH.

For each system, you can find the reference output file containing the order parameters values.

```bash

#Reference file : OUT.buildH.correct.popc_start.pdb
python ../buildH_calcOP.py -l Berger_POPC -d ../order_parameter_definitions_MODEL_Berger_POPC.def popc_start.pdb

#Reference file : OUT.buildH.correct.1POPC
python ../buildH_calcOP.py -l Berger_POPC -d ../order_parameter_definitions_MODEL_Berger_POPC.def  1POPC.pdb

#Reference file : OUT.buildH.correct.popc0-25ns_dt1000.xtc
python ../buildH_calcOP.py -l Berger_POPC -d ../order_parameter_definitions_MODEL_Berger_POPC.def popc_start.pdb -x popc0-25ns_dt1000.xtc

```
