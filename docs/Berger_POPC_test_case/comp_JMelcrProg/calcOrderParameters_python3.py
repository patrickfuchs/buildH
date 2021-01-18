#!/usr/bin/env python3

"""
 calculation of order parameters
 from a MD trajectory
 useful for example for lipid bilayers

 meant for use with NMRlipids projects

------------------------------------------------------------
 Made by J. Melcr,  Last edit 2018/03/26

 with contributions from:
  H. Antila

 edited by P. Fuchs (2020/01/18) to convert it to Python 3
------------------------------------------------------------
 input: Order parameter definitions [order_parameter_file],
        gro/tpr [topology_file], and xtc/trr/pbb file [trajectory_file] (or equivalents)
 output: order parameters [output_file] (2 textfiles)
 usage : python calcOrderParameters.py -i [order_parameter_file] -t [topology_file] -x [trajectory_file] -o [output_file]
--------------------------------------------------------
 [order_parameter_file] should contain the definitions of order parameters in format:

 OP_name1 Residue Carbon_name Hydrogen_name
 OP_name2 Residue Carbon_name Hydrogen_name
 .....

 Example (CHARMM36):

 beta1 POPC C12 H12A
 beta2 POPC C12 H12B
 alpha1 POPC C11 H11A
 alpha2 POPC C11 H11B
 g3_1 POPC C1 HA
 g3_2 POPC C1 HB
 g2_1 POPC C2 HS
 g1_1 POPC C3 HX
 g1_2 POPC C3 HY


"""

# coding: utf-8

import MDAnalysis as mda
import numpy as np
import math
import os, sys
from optparse import OptionParser
from collections import OrderedDict

bond_len_max=1.5  # in Angstroms, max distance between atoms for reasonable OP calculation (PBC and sanity check)
bond_len_max_sq=bond_len_max**2


class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.

    OP definition consist of:
       - name of the OP
       - residue name
       - involved atoms (exactly 2)
       + extra: mean, std.dev. & err. estimate
                (using standard error of the means from individual residues)
                of the OP (when reading-in an already calculated result)
    """
    def __init__(self, name, resname, atom_A_name, atom_B_name, *args):
        """
        Initialization of an instance of this class.

        it doesn't matter which atom comes first,
        atom A or B, for OP calculation.
        """
        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        # variables for error estimate -- standard error of the mean (STEM)
        self.avg   = None   # average/mean value from all residues
        self.means = None   # list of mean values from each individual residue
        self.std   = None   # standard deviation (sqrt(variance))
        self.stem  = None   # STandard Error of the Mean
        # trajectory as list
        self.traj = []  # for storing OPs
        for field in self.__dict__:
            if not isinstance(field, str):
                raise UserWarning("provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.".format(field))
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.".format(field))
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
        elif len(args) == 3:
            self.avg  = args[0]
            self.std  = args[1]
            self.stem = args[2]
        else:
            if len(args) != 0:
                raise UserWarning("Number of optional positional arguments is {len}, not 3, 2 or 0. Args: {args}\nWrong file format?".format(len=len(args), args=args))


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


    def calc_angle(self, atoms, z_dim=45.0):
        """
        calculates the angle between the vector and z-axis in degrees
        no PBC check!
        assuming a sim-box-centred membrane --> it's centre ~ z_dim/2
        Warning: user has to make sure that correct z_dim is supplied,
                 otherwise - This is a bit DIRTY!!
                 -- this is taken care of in the main trajectory reader in this module
        """
        vec = atoms[1].position - atoms[0].position
        d = math.sqrt(np.square(vec).sum())
        cos = vec[2]/d
        # values for the bottom leaflet are inverted so that
        # they have the same nomenclature as the top leaflet
        cos *= math.copysign(1.0, atoms[0].position[2]-z_dim*0.5)
        try:
            angle = math.degrees(math.acos(cos))
        except ValueError:
            if abs(cos)>=1.0:
                print("Cosine is too large = {} --> truncating it to +/-1.0".format(cos))
                cos = math.copysign(1.0, cos)
                angle = math.degrees(math.acos(cos))
        return angle

    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        This method becomes deprecated after the introduction of error estimation
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))


    @property
    def get_avg_std_stem_OP(self):
        """
        Provides average, stddev and standard error of mean for all OPs in self.traj
        """
        self.means = np.mean(self.traj, axis=0)
        return ( np.mean(self.traj),
                 np.std(self.means),
                 np.std(self.means)/np.sqrt(len(self.means)) )

def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis (mda) Universe instance with topology top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.

    ordPars : list of OrderParameter class instances
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    print(top, trajs)# ; exit()
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use with trajectory
    for op in ordPars.values():
        # selection = pairs of atoms, split-by residues
        #    this selection format preserves the order of the atoms (atA, atB) independent of their order in the topology
        selection = mol.select_atoms("resname {rnm} and name {atA}".format(
                                        rnm=op.resname, atA=op.atAname),
                                     "resname {rnm} and name {atB}".format(
                                        rnm=op.resname, atB=op.atBname)
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
        #Nres=len(op.selection)
        #Nframes=len(mol.trajectory)

    # go through trajectory frame-by-frame
    # and calculate each OP from the list of OPs
    # for each residue separately
    for frame in mol.trajectory:
        for op in ordPars.values():
            # temporary list of order parameters for
            # each individual residue for the given frame
            temp_S = []
            for residue in op.selection:
                DEBUG = False #True
                if DEBUG:
                    for atom in residue:
                        print(atom, atom.position)
                if "vec" in op.name:
                    S = op.calc_angle(residue, z_dim=frame.dimensions[2])
                else:
                    S = op.calc_OP(residue)
                if DEBUG:
                    print("OP: ", S)
                    print()
                temp_S.append(S)
            # resulting S-trajectory will be a list of lists
            # so that individual residues can be easily distinquished
            op.traj.append(temp_S)


def parse_op_input(fname):
    """
    parses input file with Order Parameter definitions
    file format is as follows:
    OP_name    resname    atom1    atom2  +extra: OP_mean  OP_std
    (flexible cols)

    fname : string
        input file name

    returns : dictionary
        with OrderParameters class instances
    """
   # Using ordered dict since it preserves the read-in order. Might come in handy when comparing to experiments.
    ordPars = OrderedDict()
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
    parser.add_option('-x', '--traj', dest='traj_fname', help='beginning of trajectory (xtc, dcd) files names (will use all files beginning with this string).', default="traj")
    parser.add_option('-o', '--out',  dest='out_fname',  help='output (OPs mean&std) file name', default="Headgroup_Glycerol_OPs.dat")
    opts, args = parser.parse_args()

    # dictionary for storing of OrderParameter class instances (name-wise, of course)
    print("\nReading OP definitions ...\n")
    ordPars = parse_op_input(opts.inp_fname)

    # get all parts of trajectories
    trajs = []
    for file_name in os.listdir(os.getcwd()):
        if file_name.startswith(opts.traj_fname):
            trajs.append(file_name)

    # read trajectory and calculate all OPs
    print("Reading trajectories and calculating OPs ...\n")
    read_trajs_calc_OPs(ordPars, opts.top_fname, trajs)


    print("OP Name     mean    std    err.est.")
    print("--------------------------------------------------------------------")
    for op in ordPars.values():
        try:
            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
            print("{:10s} {: 2.4f} {: 2.4f} {: 2.4f}".format(op.name, op.avg, op.std, op.stem))
        except:
            print("{:s} -- problem calculating statistics ".format(op.name))
    print("--------------------------------------------------------------------")


    try:
        with open(opts.out_fname,"w") as f:
            f.write("# OP_name    resname    atom1    atom2    OP_mean   OP_stddev  OP_stem\n\
#--------------------------------------------------------------------\n")
            for op in ordPars.values():
                f.write( "{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f} \n".format(
                         op.name, op.resname, op.atAname, op.atBname,
                         op.avg, op.std, op.stem)
                       )
        print("\nOrderParameters written to >> {fname} <<".format(fname=opts.out_fname))
    except:
        print("ERROR: Problems writing main output file.")


    # this single-line format may become soon deprecated, but
    # it is the format that is used in NMRlipids projects for processing through awk+gnuplot
    try:
        conc_formatted_line = "conc  {b1: 2.6f} {b1e: 2.6f}  {b2: 2.6f} {b2e: 2.6f}    {a1: 2.6f} {a1e: 2.6f}  {a2: 2.6f} {a2e: 2.6f}".format(
                              b1=ordPars['beta1'].avg, b2=ordPars['beta2'].avg,
                              a1=ordPars['alpha1'].avg, a2=ordPars['alpha2'].avg,
                              b1e=ordPars['beta1'].stem, b2e=ordPars['beta2'].stem,
                              a1e=ordPars['alpha1'].stem, a2e=ordPars['alpha2'].stem)
        print()
        print("Single line format:\nconc  beta1 err  beta2 err  alpha1 err  alpha2 err")
        print(conc_formatted_line)
        with open(opts.out_fname+".line","w") as f:
            f.write(conc_formatted_line)
    except:
        print("ERROR: Problems writing the beta-alpha single line format file.")
