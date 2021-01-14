#!/usr/bin/env python
"""
 calculation of order parameters of lipid bilayers
 from a MD trajectory
 
 meant for use with NMRlipids projects

 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/02/02
------------------------------------------------------------
 input: Order parameter definitions
        gro and xtc file (or equivalents)
 output: order parameters (2 textfiles)

--------------------------------------------------------
"""

# coding: utf-8

import MDAnalysis as mda
import numpy as np
import math
import os, sys
from optparse import OptionParser


#k_b = 0.0083144621  #kJ/Mol*K
#f_conc=55430  # factor for calculating concentrations of salts from numbers of ions/waters; in mM/L

bond_len_max=1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq=bond_len_max**2

#%%
class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """
    def __init__(self, name, resname, atom_A_name, atom_B_name, *args):
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        for field in self.__dict__:
            if not isinstance(field, str):
                raise (UserWarning, "provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.".format(field))
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.".format(field))
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
        else:
            raise UserWarning("Number of optional positional arguments is {len}, not 2 or 0. Args: {args}\nWrong file format?".format(len=len(args), args=args))
        self.traj = []  # for storing OPs


    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
        if d2>bond_len_max_sq:
            raise UserWarning("Atomic distance for atoms \
 {at1} and {at2} in residue no. {resnr} is suspiciously \
 long: {d}!\nPBC removed???".format(at1=atoms[0].name, at2=atoms[1].name, resnr=atoms[0].resid, d=math.sqrt(d2)))
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        return S


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))


def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.

    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use in trajectory
    for op in ordPars.values():
        # selection = pairs of atoms, split-by residues
        selection = mol.select_atoms("resname {rnm} and name {atA} {atB}".format(
                                    rnm=op.resname, atA=op.atAname, atB=op.atBname)
                                    ).atoms.split("residue")
        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                print(res.resnames, res.resids)
                for atom in res.atoms:
                    print(atom.name, atom.id)
                raise UserWarning("Selection >> name {atA} {atB} << \
                contains {nat} atoms, but should contain exactly 2!".format(
                atA=op.atAname, atB=op.atBname, nat=res.n_atoms))
        op.selection = selection

    # go through trajectory frame-by-frame
    for frame in mol.trajectory:
        for op in ordPars.values():
            for residue in op.selection:
                S = op.calc_OP(residue)
                op.traj.append(S)
        #print "--", mol.atoms[0].position


def parse_op_input(fname):
    """
    parses input file with Order Parameter definitions
    file format is as follows:
    OP_name    resname    atom1    atom2
    (flexible cols)

    fname : string
        input file name

    returns : dictionary 
        with OrderParameters class instances
    """
    ordPars = {}
    try:
        with open(fname,"r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    items = line.split()
                    ordPars[items[0]] = OrderParameter(*items)
    except:
        raise RuntimeError("Couldn't read input file >> {inpf} <<".format(inpf=opts.inp_fname))
    return ordPars



#%%

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-i', '--inp',  dest='inp_fname',  help='input (OP definitions) file name', default="Headgroup_Glycerol_OPs.def")
    parser.add_option('-t', '--top',  dest='top_fname',  help='topology (gro, pdb) file name', default="last_frame_nonwat.gro")
    parser.add_option('-x', '--traj', dest='traj_fname', help='beginning of trajectory (xtc) files names (will use all that begin with this).', default="traj")
    parser.add_option('-o', '--out',  dest='out_fname',  help='output (OPs mean&std) file name', default="Headgroup_Glycerol_OPs.dat")
    opts, args = parser.parse_args()

#%%
    # dictionary for storing of OrderParameter class instances (name-wise, of course)
    print("\nReading OP definitions ...\n")
    ordPars = parse_op_input(opts.inp_fname)

#%%
    # get all parts of trajectories
    trajs = []
    for file_name in os.listdir(os.getcwd()):
             if file_name.startswith(opts.traj_fname):
                 trajs.append(file_name)

    # read trajectory and calculate all OPs
    print("Reading trajectories and calculating OPs ...\n")
    read_trajs_calc_OPs(ordPars, opts.top_fname, trajs)

#%%

    print("OP Name     mean     stddev ")
    print("----------------------------")
    for op in ordPars.values():
        (op.avg, op.std) = op.get_avg_std_OP
        print(op.name, op.avg, op.std)


    try:
        with open(opts.out_fname,"w") as f:
            f.write("# OP_name    resname    atom1    atom2    OP_mean   OP_stddev\n\
#-------------------------------------------------------------\n")
            for op in ordPars.values():
                f.write( "   ".join([op.name, op.resname, op.atAname, op.atBname, str(op.avg), str(op.std), "\n"]) )
        print("\nOrderParameters written to >> {fname} <<".format(fname=opts.out_fname))
    except:
        print("ERROR: Problems writing main output file.")


    #conc = f_conc*nion/nwat
    try:
        conc_formatted_line = "conc  {b1} 0  {b2} 0    {a1} 0  {a2} 0".format(
                              b1=ordPars['beta1'].avg, b2=ordPars['beta2'].avg,
                              a1=ordPars['alpha1'].avg, a2=ordPars['alpha2'].avg)
        print
        print("Single line format:\nconc  beta1 0  beta2 0  alpha1 0  alpha2 0")
        print(conc_formatted_line)
        with open(opts.out_fname+".line","w") as f:
            f.write(conc_formatted_line)
    except:
        print("ERROR: Problems writing the beta-alpha single line format file.")

