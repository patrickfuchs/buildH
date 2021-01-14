#!/bin/bash

# Launch buildH to reconstruct traj w Hs
buildH ../start_128popc.pdb -l Berger_POPC \
       -d ../order_parameter_definitions_MODEL_Berger_POPC.def \
       -x ../popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_H

# Rename files
mv popc0-25ns_dt1000_H.pdb start_128popc_H.pdb

# Remove buildH OP output
rm OP_buildH.out

# Get prog from JMelcr on NMRlipids (last version in Python 3, downloaded on Jan 13 2020)
wget https://raw.githubusercontent.com/NMRLipids/MATCH/master/scratch/scriptsBYmelcr/calcOrderParametersPYTHON3.py

# Launch JMelcr Prog
python ./calcOrderParametersPYTHON3.py -i ../order_parameter_definitions_MODEL_Berger_POPC.def \
       -t start_128popc_H.pdb -x popc0-25ns_dt1000_H.xtc -o OUT.progJMeclr.popc0-25ns_dt1000_H.out

# Comp with R
R --vanilla < comp_JMelcrProg.R
