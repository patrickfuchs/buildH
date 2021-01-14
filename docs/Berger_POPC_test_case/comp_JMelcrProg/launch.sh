#!/bin/bash

# Launch buildH to reconstruct traj w Hs
buildH ../start_128popc.pdb -l Berger_POPC -d ../order_parameter_definitions_MODEL_Berger_POPC.def -x ../popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_H.xtc

# Rename files
mv popc0-25ns_dt1000_H.xtc.pdb start_128popc_H.pdb
mv popc0-25ns_dt1000_H.xtc.xtc popc0-25ns_dt1000_H.xtc

# Remove buildH OP output
rm OP_buildH.out

# Get prog from JMelcr on NMRlipids (last version in Python 3)
wget https://raw.githubusercontent.com/NMRLipids/MATCH/master/scratch/scriptsBYmelcr/calcOrderParametersPYTHON3.py

# Launch JMelcr Prog
python ./calcOrderParametersPYTHON3.py -i ../order_parameter_definitions_MODEL_Berger_POPC.def \
       -t start_128popc_H.pdb -x popc0-25ns_dt1000_H.xtc -o OUT.progJMeclr.popc0-25ns_dt1000_H.out

# Comp with R
R --vanilla < comp_JMelcrProg.R
