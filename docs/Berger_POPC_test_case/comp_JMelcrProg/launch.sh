#!/bin/bash

# Launch buildH to reconstruct traj w Hs.
buildH -c ../start_128popc.pdb -l Berger_POPC \
       -d ../order_parameter_definitions_MODEL_Berger_POPC.def \
       -t ../popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_H

# Rename file.
mv popc0-25ns_dt1000_H.pdb start_128popc_H.pdb

# Remove buildH OP output.
rm OP_buildH.out

# Get prog from JMelcr on NMRlipids (last version in Python 3, downloaded on Jan 13 2020).
# The specific commit where this file has been pushed is here: https://raw.githubusercontent.com/NMRLipids/MATCH/40dc2d29437cc8b9a442fb6e46fb73731a1662fa/scripts/calcOrderParameters.py
wget https://raw.githubusercontent.com/NMRLipids/MATCH/master/scripts/calcOrderParameters.py

# Script calcOrderParameters.py modified by P. Fuchs to make it work with Python 3 and renamed to calcOrderParameters_python3.py.

# Launch JMelcr Prog.
python ./calcOrderParameters_python3.py -i ../order_parameter_definitions_MODEL_Berger_POPC.def \
       -t start_128popc_H.pdb -x popc0-25ns_dt1000_H.xtc -o OUT.progJMeclr.popc0-25ns_dt1000_H.out

# Remove file.
OUT.progJMeclr.popc0-25ns_dt1000_H.out.line

# To see numerical comparison.
meld OUT.progJMeclr.popc0-25ns_dt1000_H.out ../OUT.buildH.popc0-25ns_dt1000.xtc.out

# Comparison graph with R.
R --vanilla < comp_JMelcrProg.R
