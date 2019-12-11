#!/usr/bin/env python3
# coding: utf-8

"""
This script builds hydrogens from a united-atom trajectory and calculate the
order parameter for each C-H bond.

It works in two modes :
  1) A slow mode when an output trajectory (e.g. in xtc format) is requested by
     the user. In this case, the whole trajectory including newly built
     hydrogens are written to this trajectory.
  2) A fast mode without any output trajectory.
For both modes, the order parameter is written to an output file in a format
similar to the code of @jmelcr:
https://github.com/NMRLipids/MATCH/blob/master/scripts/calcOrderParameters.py

This code has been checked against the one from @jmelcr. You might find minor
differences due to rounding errors (in xtc, only 3 digits are written).

The way of building H is largely inspired from a code of Jon Kapla originally
written in fortran :
https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242.

Note: that all coordinates in this script are handled using numpy 1D-arrays
of 3 elements, e.g. atom_coor = np.array((x, y, z)).
Note2: sometimes numpy is slow on small arrays, thus we wrote a few "in-house"
functions for vectorial operations (e.g. cross product).
"""

__authors__ = ("Patrick Fuchs", "Amélie Bâcle", "Hubert Santuz",
               "Pierre Poulain")
__contact__ = ("patrickfuchs", "abacle", "hublot", "pierrepo") # on github
__version__ = "1.0.3"
__copyright__ = "BSD"
__date__ = "2019/05"

# Modules.
import argparse
import copy
import collections
import pickle
import sys
import warnings
import os

import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.coordinates.XTC as XTC

import dic_lipids


# Constants.
# From https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond
LENGTH_CH_BOND = 1.09 # in Angst
# From https://en.wikipedia.org/wiki/Tetrahedron, tetrahedral angle equals
# arccos(-1/3) ~ 1.9106 rad ~ 109.47 deg.
TETRAHEDRAL_ANGLE = np.arccos(-1/3)
# For debugging.
DEBUG = False
# For pickling results (useful for future analyses, e.g. drawing distributions).
PICKLE = False

def isfile(path):
    """Check if path is an existing file.
    If not, raise an error. Else, return the path."""
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path



def calc_OP(C, H):
    """Returns the Order Parameter of a CH bond (OP).

    OP is calculated according to equation:

    S = 1/2 * (3*cos(theta)^2 -1)

    theta is the angle between CH bond and the z(vertical) axis:
    z
    ^  H
    | /
    |/
    C

    Inspired from a function written by @jmelcr.

    Parameters
    ----------
    C : numpy 1D-array
        Coordinates of C atom.
    H : numpy 1D-array
        Coordinates of H atom.

    Returns
    -------
    float
        The normalized vector.
    """
    vec = H - C
    d2 = np.square(vec).sum()
    cos2 = vec[2]**2/d2
    S = 0.5*(3.0*cos2 - 1.0)
    return S


def normalize(vec):
    """Normalizes a vector.

    Parameters
    ----------
    vec : numpy 1D-array

    Returns
    -------
    numpy 1D-array
        The normalized vector.
    """
    return vec / norm(vec)


def norm(vec):
    """Returns the norm of a vector.

    Parameters
    ----------
    vec : numpy 1D-array

    Returns
    -------
    float
        The magniture of the vector.
    """
    return np.sqrt(np.sum(vec**2))


def calc_angle(atom1, atom2, atom3):
    """Calculates the valence angle between atom1, atom2 and atom3.

    Note: atom2 is the central atom.

    Parameters
    ----------
    atom1 : numpy 1D-array.
    atom2 : numpy 1D-array.
    atom3 : numpy 1D-array.

    Returns
    -------
    float
        The calculated angle in radians.
    """
    vec1 = atom1 - atom2
    vec2 = atom3 - atom2
    costheta = np.dot(vec1,vec2)/(norm(vec1)*norm(vec2))
    if costheta > 1.0 or costheta < -1.0:
        raise ValueError("Cosine cannot be larger than 1.0 or less than -1.0")
    return np.arccos(costheta)


def vec2quaternion(vec, theta):
    """Translates a vector of 3 elements and angle theta to a quaternion.

    Parameters
    ----------
    vec : numpy 1D-array
        Vector of the quaternion.
    theta : float
        Angle of the quaternion in radian.

    Returns
    -------
    numpy 1D-array
        The full quaternion (4 elements).
    """
    w = np.cos(theta/2)
    x, y, z = np.sin(theta/2) * normalize(vec)
    q = np.array([w, x, y, z])
    return q


def calc_rotation_matrix(quaternion):
    """Translates a quaternion to a rotation matrix.

    Parameters
    ----------
    quaternion : numpy 1D-array of 4 elements.

    Returns
    -------
    numpy 2D-array (dimension [3, 3])
        The rotation matrix.
    """
    # Initialize rotation matrix.
    matrix = np.zeros([3, 3])
    # Get quaternion elements.
    w, x, y, z = quaternion
    # Compute rotation matrix.
    matrix[0,0] = w**2 + x**2 - y**2 - z**2
    matrix[1,0] = 2 * (x*y + w*z)
    matrix[2,0] = 2 * (x*z - w*y)
    matrix[0,1] = 2 * (x*y - w*z)
    matrix[1,1] = w**2 - x**2 + y**2 - z**2
    matrix[2,1] = 2 * (y*z + w*x)
    matrix[0,2] = 2 * (x*z + w*y)
    matrix[1,2] = 2 * (y*z - w*x)
    matrix[2,2] = w**2 - x**2 - y**2 + z**2
    return matrix


def apply_rotation(vec_to_rotate, rotation_axis, rad_angle):
    """Rotates a vector around an axis by a given angle.

    Note: the rotation axis is a vector of 3 elements.

    Parameters
    ----------
    vec_to_rotate : numpy 1D-array
    rotation_axis : numpy 1D-array
    rad_angle : float

    Returns
    -------
    numpy 1D-array
        The final rotated (normalized) vector.
    """
    # Generate a quaternion of the given angle (in radian).
    quaternion = vec2quaternion(rotation_axis, rad_angle)
    # Generate the rotation matrix.
    rotation_matrix = calc_rotation_matrix(quaternion)
    # Apply the rotation matrix on the vector to rotate.
    vec_rotated = np.dot(rotation_matrix, vec_to_rotate)
    return normalize(vec_rotated)


def pandasdf2pdb(df):
    """Returns a string in PDB format from a pandas dataframe.

    Parameters
    ----------
    df : pandas dataframe with columns "atnum", "atname", "resname", "resnum",
         "x", "y", "z"

    Returns
    -------
    str
        A string representing the PDB.
    """
    s = ""
    chain = ""
    for _, row_atom in df.iterrows():
        atnum, atname, resname, resnum, x, y, z = row_atom
        atnum = int(atnum)
        resnum = int(resnum)
        # See for pdb format:
        # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
        # "alt" means alternate location indicator
        # "code" means code for insertions of residues
    	# "seg" means segment identifier
        # "elt" means element symbol
        if len(atname) == 4:
            s += ("{record_type:6s}{atnum:5d} {atname:<4s}{alt:1s}{resname:>4s}"
                  "{chain:1s}{resnum:>4d}{code:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                  "{occupancy:>6.2f}{temp_fact:>6.2f}          {seg:<2s}{elt:>2s}\n"
                  .format(record_type="ATOM", atnum=atnum, atname=atname, alt="",
                          resname=resname, chain=chain, resnum=resnum, code="",
                          x=x, y=y, z=z, occupancy=1.0, temp_fact=0.0, seg="",
                          elt=atname[0]))
        else:
            s += ("{record_type:6s}{atnum:5d}  {atname:<3s}{alt:1s}{resname:>4s}"
                  "{chain:1s}{resnum:>4d}{code:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                  "{occupancy:>6.2f}{temp_fact:>6.2f}          {seg:<2s}{elt:>2s}\n"
                  .format(record_type="ATOM", atnum=atnum, atname=atname, alt="",
                          resname=resname, chain=chain, resnum=resnum, code="",
                          x=x, y=y, z=z, occupancy=1.0, temp_fact=0.0, seg="",
                          elt=atname[0]))
    return s


def cross_product(A, B):
    """Returns the cross product between vectors A & B.

    Source: http://hyperphysics.phy-astr.gsu.edu/hbase/vvec.html.
    Note: on small vectors (i.e. of 3 elements), computing cross products
          with this functions is faster than np.cross().

    Parameters
    ----------
    A : numpy 1D-array
        A vector of 3 elements.
    B : numpy 1D-array
        Another vector of 3 elements.

    Returns
    -------
    numpy 1D-array
        Cross product of A^B.
    """
    x = (A[1]*B[2]) - (A[2]*B[1])
    y = (A[0]*B[2]) - (A[2]*B[0])
    z = (A[0]*B[1]) - (A[1]*B[0])
    return np.array((x, -y, z))


def get_CH2(atom, helper1, helper2):
    """Reconstructs the 2 hydrogens of a sp3 carbon (methylene group).

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct hydrogens.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom after central atom.

    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the two hydrogens:
        ([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2]).
    """
    # atom->helper1 vector.
    v2 = normalize(helper1 - atom)
    # atom->helper2 vector.
    v3 = normalize(helper2 - atom)
    # Vector orthogonal to the helpers/atom plane.
    #v4 = normalize(np.cross(v3, v2))
    v4 = normalize(cross_product(v3, v2))
    # Rotation axis is atom->helper1 vec minus atom->helper2 vec.
    rotation_axis = normalize(v2 - v3)
    # Vector to be rotated by theta/2, perpendicular to rotation axis and v4.
    #vec_to_rotate = normalize(np.cross(v4, rotation_axis))
    vec_to_rotate = normalize(cross_product(v4, rotation_axis))
    # Reconstruct the two hydrogens.
    unit_vect_H1 = apply_rotation(vec_to_rotate, rotation_axis,
                                 -TETRAHEDRAL_ANGLE/2)
    hcoor_H1 = LENGTH_CH_BOND * unit_vect_H1 + atom
    unit_vect_H2 = apply_rotation(vec_to_rotate, rotation_axis,
                                 TETRAHEDRAL_ANGLE/2)
    hcoor_H2 = LENGTH_CH_BOND * unit_vect_H2 + atom
    return (hcoor_H1, hcoor_H2)


def get_CH(atom, helper1, helper2, helper3):
    """Reconstructs the unique hydrogen of a sp3 carbon.

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct the hydrogen.
    helper1 : numpy 1D-array
        First neighbor of central atom.
    helper2 : numpy 1D-array
        Second neighbor of central atom.
    helper3 : numpy 1D-array
        Third neighbor of central atom.

    Returns
    -------
    numpy 1D-array
        Coordinates of the rebuilt hydrogen: ([x_H, y_H, z_H]).
    """
    helpers = np.array((helper1, helper2, helper3))
    v2 = np.zeros(3)
    for i in range(len(helpers)):
        v2 = v2 + normalize(helpers[i] - atom)
    v2 = v2 / (len(helpers)) + atom
    unit_vect_H = normalize(atom - v2)
    coor_H = LENGTH_CH_BOND * unit_vect_H + atom
    return coor_H


def get_CH_double_bond(atom, helper1, helper2):
    """Reconstructs the hydrogen of a sp2 carbon.

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct the hydrogen.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom after central atom.

    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the rebuilt hydrogen: ([x_H, y_H, z_H]).
    """
    # calc CCC_angle helper1-atom-helper2 (in rad).
    CCC_angle = calc_angle(helper1, atom, helper2)
    # We want to bisect the C-C-C angle ==> we take half of (2pi-CCC_angle).
    # Factorizing yields: pi - CCC_angle/2.
    theta = np.pi - (CCC_angle / 2)
    # atom->helper1 vector.
    v2 = helper1 - atom
    # atom->helper2 vector.
    v3 = helper2 - atom
    # The rotation axis is orthogonal to the atom/helpers plane.
    #rotation_axis = normalize(np.cross(v2, v3))
    rotation_axis = normalize(cross_product(v2, v3))
    # Reconstruct H by rotating v3 by theta.
    unit_vect_H = apply_rotation(v3, rotation_axis, theta)
    coor_H = LENGTH_CH_BOND * unit_vect_H + atom
    return coor_H


def get_CH3(atom, helper1, helper2):
    """Reconstructs the 3 hydrogens of a sp3 carbon (methyl group).

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct hydrogens.
    helper1 : numpy 1D-array
        Heavy atom before central atom.
    helper2 : numpy 1D-array
        Heavy atom before helper1 (two atoms away from central atom).

    Returns
    -------
    tuple of numpy 1D-arrays
        Coordinates of the 3 hydrogens:
        ([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2], [x_H3, y_H3, z_H3]).
    """
    ### Build CH3e.
    theta = TETRAHEDRAL_ANGLE
    # atom->helper1 vector.
    v2 = helper1 - atom
    # atom->helper2 vector.
    v3 = helper2 - atom
    # Rotation axis is perpendicular to the atom/helpers plane.
    #rotation_axis = normalize(np.cross(v3, v2))
    rotation_axis = normalize(cross_product(v3, v2))
    # Rotate v2 by tetrahedral angle. New He will be in the same plane
    # as atom and helpers.
    unit_vect_He = apply_rotation(v2, rotation_axis, theta)
    coor_He = LENGTH_CH_BOND * unit_vect_He + atom
    ### Build CH3r.
    theta = (2/3) * np.pi
    rotation_axis = normalize(helper1 - atom)
    v4 = normalize(coor_He - atom)
    # Now we rotate atom->He bond around atom->helper1 bond by 2pi/3.
    unit_vect_Hr = apply_rotation(v4, rotation_axis, theta)
    coor_Hr = LENGTH_CH_BOND * unit_vect_Hr + atom
    ### Build CH3s.
    theta = -(2/3) * np.pi
    rotation_axis = normalize(helper1 - atom)
    v5 = normalize(coor_He - atom)
    # Last we rotate atom->He bond around atom->helper1 bond by -2pi/3.
    unit_vect_Hs = apply_rotation(v5, rotation_axis, theta)
    coor_Hs = LENGTH_CH_BOND * unit_vect_Hs + atom
    return coor_He, coor_Hr, coor_Hs


###
### The next two functions (buildHs_on_1C() and build_all_Hs_calc_OP())
### build new H, calculate the order parameter and write the new traj with Hs
### to an output file (e.g. .xtc, etc).
### Note: they are slow, they shouldn't be used if the user doesn't want to
###       write the trajectory. Instead, fast_build_all_Hs() should be used.
###
def buildHs_on_1C(atom, dic_lipid):
    """Builds 1, 2 or 3 H on a given carbon.

    This function is a wrapper which gathers the coordinates of the helpers
    and call the function that builds 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described
    in a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    atom : MDAnalysis Atom instance
        This instance contains the carbon on which we want to build Hs.
    dic_lipid : dictionnary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if
        there's only one H that has been rebuilt.
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipid[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipid[atom.name]
    else:
        typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipid[atom.name]
    # Get helper coordinates using atom, which an instance from Atom class.
    # atom.residue.atoms is a list of atoms we can select with
    # method .select_atoms().
    # To avoid too long line, we shorten its name to `sel`.
    sel = atom.residue.atoms.select_atoms
    helper1_coor = sel("name {0}".format(helper1_name))[0].position
    helper2_coor = sel("name {0}".format(helper2_name))[0].position
    if typeofH2build == "CH2":
        H1_coor, H2_coor = get_CH2(atom.position, helper1_coor, helper2_coor)
        return (H1_coor, H2_coor)
    elif typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_coor = sel("name {0}".format(helper3_name))[0].position
        H1_coor = get_CH(atom.position, helper1_coor, helper2_coor,
                         helper3_coor)
        return (H1_coor,)
    elif typeofH2build == "CHdoublebond":
        H1_coor = get_CH_double_bond(atom.position, helper1_coor,
                                     helper2_coor)
        return (H1_coor,)
    elif typeofH2build == "CH3":
        H1_coor, H2_coor, H3_coor = get_CH3(atom.position,
                                            helper1_coor, helper2_coor)
        return (H1_coor, H2_coor, H3_coor)
    else:
        raise UserWarning("Wrong code for typeofH2build, expected 'CH2', 'CH'"
                          ", 'CHdoublebond' or 'CH3', got {}."
                          .format(typeofH2build))


def build_all_Hs_calc_OP(universe_woH, dic_lipid, dic_Cname2Hnames,
                         universe_wH=None, dic_OP=None,
                         dic_corresp_numres_index_dic_OP=None,
                         return_coors=False):
    """Main function that builds all hydrogens and calculates order parameters.

    This function shall be used in two modes :

    1) The first time this function is called, we have to construct a new
    universe with hydrogens. One shall call it like this :

    new_data_frame = build_all_Hs_calc_OP(universe_woH, return_coors=True)

    The boolean return_coors set to True indicates to the function to return
    a pandas dataframe. This latter will be used later to build a new
    universe with H.

    2) For all the other frames, we just need to update the coordinates in
    the universe *with* hydrogens. One shall call it like this :

    build_all_Hs_calc_OP(universe_woH, dic_lipid, dic_Cname2Hnames,
                         universe_wH=universe_wH, dic_OP=dic_OP,
                         dic_corresp_numres_index_dic_OP=dic_corresp_numres_index_dic_OP)

    In this case, the function also calculates the order parameter and returns
    nothing. The coordinates of the universe *with* H are updated in place.
    The order parameter is also added in place (within dic_OP dictionnary).

    NOTE: This function in mode 2 is slow, thus it shall be used when one wants
    to create a trajectory with H (such as .xtc or whatever format).

    NOTE2: This function assumes all possible C-H pairs are present in the .def 
    file (with -d option). They are needed since we want to build an xtc with 
    the whole system. If one is interested in calculating only a subset of OPs,
    please use the function fast_build_all_Hs_calc_OP() instead.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_lipid : dictionnary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionnary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}
    universe_wH : MDAnalysis universe instance (optional)
        This is the universe *with* hydrogens.
    dic_OP : ordered dictionnary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.
    dic_corresp_numres_index_dic_OP : dictionnary
        This dict should contain the correspondance between the numres and
        the corresponding index in dic_OP. For example {..., 15: 14, ...} means
        the residue numbered 15 in the PDB has an index of 14 in dic_OP.
    return_coors : boolean (optional)
        If True, the function will return a pandas dataframe containing the
        system *with* hydrogens.

    Returns
    -------
    pandas dataframe (optional)
        If parameter return_coors is True, this dataframe contains the
        system *with* hydrogens is returned.
    None
        If parameter return_coors is False.
    """
    if universe_wH:
        # We will need the index in the numpy array for updating coordinates
        # in the universe with H.
        row_index_coor_array = 0
    if return_coors:
        # The list newrows will be used to store the new molecule *with* H.
        newrows = []
        # Counter for numbering the new mlcs with H.
        new_atom_num = 1
    # Loop over all atoms in the universe without H.
    for atom in universe_woH.atoms:
        if universe_wH:
            # Update the position of the current atom in the universe with H.
            universe_wH.coord.positions[row_index_coor_array, :] = atom.position
            row_index_coor_array += 1
        if return_coors:
            resnum = atom.resnum
            resname = atom.resname
            name = atom.name
            # Append atom to the new list.
            # 0      1       2        3       4  5  6
            # atnum, atname, resname, resnum, x, y, z
            newrows.append([new_atom_num, name, resname, resnum]
                           + list(atom.position))
            new_atom_num += 1
        # Build new H(s)?
        if (atom.name in dic_lipid and
            atom.residue.resname == dic_lipid["resname"]):
            # Build Hs and store them in a list of numpy 1D-arrays Hs_coor.
            # The "s" in Hs_coor means there can be more than 1 H:
            # For CH2, Hs_coor will contain: [H1_coor, H2_coor].
            # For CH3, Hs_coor will contain: [H1_coor, H2_coor, H3_coor].
            # For CH, Hs_coor will contain: [H1_coor].
            # For CHdoublebond, Hs_coor will contain: [H1_coor].
            Hs_coor = buildHs_on_1C(atom, dic_lipid)
            # To retrieve Hname, we need a counter.
            counter4Hname = 0
            # Loop over Hs_coor (H_coor is a 1D-array with the 3 coors of 1 H).
            for i, H_coor in enumerate(Hs_coor):
                # Retrieve name of newly built H.
                Hname = dic_Cname2Hnames[atom.name][counter4Hname]
                ####
                #### We calculate here the order param on the fly :-D !
                ####
                if dic_OP:
                    if (atom.name, Hname) in dic_OP:
                        op = calc_OP(atom.position, H_coor)
                        # We should get here the index of the residue in dic_OP.
                        # For that we can use dic_corresp_numres_index_dic_OP
                        # (key: resnum in pdb, value: index residue in dic_OP).
                        lipid_ix = dic_corresp_numres_index_dic_OP[atom.resid]
                        # OLD way: dic_OP[(atom.name, Hname)].append(op)
                        if (atom.name, Hname) in dic_OP:
                            dic_OP[(atom.name, Hname)][lipid_ix].append(op)
                        if DEBUG:
                            print(atom.name, H_coor, "OP:", op)
                if return_coors:
                    # Add them to newrows.
                    newrows.append([new_atom_num, Hname, resname, resnum]
                                   + list(H_coor))
                    new_atom_num += 1
                if universe_wH:
                    # Update the position of the current H in the universe with H.
                    universe_wH.coord.positions[row_index_coor_array, :] = H_coor
                    row_index_coor_array += 1
                # Increment counter4Hname for retrieving next H.
                counter4Hname += 1
            if dic_OP and DEBUG:
                print() ; print()
    if dic_OP and DEBUG:
        print("Final dic_OP:", dic_OP)
        print()
    if return_coors:
        # Create a dataframe to store the mlc with added hydrogens.
        new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname",
                                                      "resname", "resnum",
                                                      "x", "y", "z"])
        return new_df_atoms

###
### The next 4 functions (fast_build_all_Hs_calc_OP(), fast_buildHs_on_1C(),
### make_dic_lipids_with_indexes() and get_indexes()) should be used when the
### user doesn't want an output trajectory.
### By using fast indexing to individual Catoms and helpers, they
### are much faster.
###
def fast_buildHs_on_1C(dic_lipids_with_indexes, ts, Cname, ix_first_atom_res):
    """Builds 1, 2 or 3 H on a given carbon using fast indexing.

    This function is a fast wrapper which gathers the coordinates of the
    helpers and call the function that builds 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described
    in a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    dic_lipids_with_indexes : dictionnary
        The dictionnary made in function make_dic_lipids_with_indexes().
    ts : MDAnalysis Timestep instance
        This object contains the actual frame under analysis.
    Cname : str
        The carbon name on which we want to build H(s).
    ix_first_atom_res : int
        The index of the first atom in the lipid under analysis.

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if
        there's only one H that has been rebuilt.
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipids_with_indexes[Cname]) == 6:
        typeofH2build, _, _, Cname_ix, helper1_ix, helper2_ix = dic_lipids_with_indexes[Cname]
    else:
        typeofH2build, _, _, _, Cname_ix, helper1_ix, helper2_ix, helper3_ix = dic_lipids_with_indexes[Cname]
    # Get Cname coordinates.
    Cname_position = ts[Cname_ix+ix_first_atom_res]
    # Get helper coordinates
    helper1_coor = ts[helper1_ix+ix_first_atom_res]
    helper2_coor = ts[helper2_ix+ix_first_atom_res]
    # Build new H(s) and get coordinates.
    if typeofH2build == "CH2":
        H1_coor, H2_coor = get_CH2(Cname_position, helper1_coor, helper2_coor)
        return (H1_coor, H2_coor)
    elif typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_coor = ts[helper3_ix+ix_first_atom_res]
        H1_coor = get_CH(Cname_position, helper1_coor, helper2_coor,
                         helper3_coor)
        return (H1_coor,)
    elif typeofH2build == "CHdoublebond":
        H1_coor = get_CH_double_bond(Cname_position, helper1_coor,
                                     helper2_coor)
        return (H1_coor,)
    elif typeofH2build == "CH3":
        H1_coor, H2_coor, H3_coor = get_CH3(Cname_position,
                                            helper1_coor, helper2_coor)
        return (H1_coor, H2_coor, H3_coor)
    else:
        raise UserWarning("Wrong code for typeofH2build, expected 'CH2', 'CH'"
                          ", 'CHdoublebond' or 'CH3', got {}."
                          .format(typeofH2build))


def get_indexes(atom, universe_woH, dic_lipid):
    """Returns the index of helpers for a given carbon.

    Parameters
    ----------
    atom : MDAnalysis Atom instance
        This is an Atom instance of a carbon on which we want to build Hs.
    universe_woH : MDAnalysis universe instance
        The universe without hydrogens.
    dic_lipid : dictionnary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.

    Returns
    -------
    tuple of 2 or 3 int
        The tuple contains the index of the 2 (or 3) helpers for the atom that
        was passed as argument. (e.g. for atom C37 with index 99, the function
        returns a tuple containing 98 (index of C36 = helper 1) and 100 (index
        of C38=helper2).
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipid[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipid[atom.name]
    else:
        typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipid[atom.name]
    # Get helper coordinates using atom, which an instance from Atom class.
    # atom.residue.atoms is a list of atoms we can select with
    # method .select_atoms().
    # To avoid too long line, we shorten its name to `sel`.
    sel = atom.residue.atoms.select_atoms
    helper1_ix = sel("name {}".format(helper1_name))[0].ix
    helper2_ix = sel("name {}".format(helper2_name))[0].ix
    if typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_ix = sel("name {0}".format(helper3_name))[0].ix
        return (helper1_ix, helper2_ix, helper3_ix)
    else:
        return (helper1_ix, helper2_ix)


def make_dic_lipids_with_indexes(universe_woH, dic_lipid):
    """Expands dic_lipid and adds the index of each atom and helper.

    IMPORTANT: the index of each atom/helper is given with respect to the
               first atom in that residue.
    For example, if we have a POPC where C1 is the first atom, and C50 the
    last one, we want in the end:
    {'C1': ('CH3', 'N4', 'C5', 0, 3, 4), ...,
     'C50': ('CH3', 'C49', 'C48', 49, 48, 47)}
    Where the 3 last int are the index (ix) of the atom, helper1, helper2
    (possibly helper3) with respect to the first atom.
    Thus for C1 : 0 is index of C1, N4 is 3 atoms away from C1 and C5 is 4
    atoms away from C1.
    For C50: C50 is 49 atoms away from C1, C49 is 48 atoms away from C1,
    C48 is 47 atoms away from C1.

    Parameters
    ----------
    universe_woH : MDAnalysis Universe instance
        The universe without hydrogens.
    dic_lipid : dictionnary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.

    Returns
    -------
    dictionnary
        The returned dictionnary as described above in this docstring.
    """
    # Get lipid name.
    resname = dic_lipid["resname"]
    # Get resnum of the 1st lipid encountered in the system whose name
    # is `resname`.
    selection = "resname {}".format(resname)
    first_lipid_residue = universe_woH.select_atoms(selection).residues[0]
    resnum_1st_lipid = first_lipid_residue.resnum
    # Get name of 1st atom of that lipid.
    first_atom_name = first_lipid_residue.atoms[0].name
    # Get index of this atom.
    first_atom_ix = first_lipid_residue.atoms[0].ix
    if DEBUG:
        print("resname: {}, first encountered residue: {},\n"
              "resnum_1st_lipid: {}, first_atom_name: {}, first_atom_ix: {}"
              .format(resname, first_lipid_residue, resnum_1st_lipid,
                      first_atom_name, first_atom_ix))
        print()
    # Keep only carbons on which we want to build Hs.
    carbons2keep = []
    for Cname, Hname in dic_OP:
        if Cname not in carbons2keep:
            carbons2keep.append(Cname)
    dic_lipids_with_indexes = {}
    for Cname in dic_lipid.keys():
        if Cname in carbons2keep:
            dic_lipids_with_indexes[Cname] = dic_lipid[Cname]
    # Now add the helper indexes.
    # The reasonning is over one residue (e.g. POPC). We want to add (to the
    # dict) the index (ix) of each helper of a given carbon with respect to
    # the index of the first atom in that lipid residue.
    # Loop over each carbon on which we want to reconstruct Hs.
    for Cname in dic_lipids_with_indexes.keys():
        # Loop over residues for a given Cname atom.
        selection = "resid {} and name {}".format(resnum_1st_lipid, Cname)
        for Catom in universe_woH.select_atoms(selection):
            # Get the (absolute) index of helpers.
            if dic_lipid[Cname][0] == "CH":
                helper1_ix, helper2_ix, helper3_ix = get_indexes(Catom,
                                                                 universe_woH,
                                                                 dic_lipid)
            else:
                helper1_ix, helper2_ix = get_indexes(Catom, universe_woH,
                                                     dic_lipid)
            # If the first lipid doesn't start at residue 1 we must
            # substract the index of the first atom of that lipid.
            Catom_ix_inres = Catom.ix - first_atom_ix
            helper1_ix_inres = helper1_ix - first_atom_ix
            helper2_ix_inres = helper2_ix - first_atom_ix
            # Then add these indexes to dic_lipids_with_indexes.
            if dic_lipid[Cname][0] == "CH":
                helper3_ix_inres = helper3_ix - first_atom_ix
                tmp_tuple = (Catom_ix_inres, helper1_ix_inres,
                             helper2_ix_inres, helper3_ix_inres)
                dic_lipids_with_indexes[Cname] += tmp_tuple
            else:
                tmp_tuple = (Catom_ix_inres, helper1_ix_inres,
                             helper2_ix_inres)
                dic_lipids_with_indexes[Cname] += tmp_tuple
    if DEBUG:
        print("Everything is based on the following dic_lipids_with_indexes\n{}"
              .format(dic_lipids_with_indexes))
        print()
    return dic_lipids_with_indexes


def fast_build_all_Hs_calc_OP(universe_woH, dic_OP, dic_lipid, dic_Cname2Hnames):
    """Build Hs and calc OP using fast indexing.

    This function uses fast indexing to carbon atoms and helper atoms. It
    should be used when the user doesn't want any output traj with hydrogens.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_OP : ordered dictionnary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.
    dic_lipid : dictionnary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionnary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}

    Returns
    -------
    None
        This function returns nothing, dic_OP is changed *in place*.
    """
    ###
    ### 1) Expand dic_lipids and store there helpers' index.
    ###
    ### We want {'C1': ('CH3', 'N4', 'C5', 0, 3, 4), ...,
    ###          'C50': ('CH3', 'C49', 'C48', 49, 48, 47)}
    ### Where the 3 last int are the index (ix) of the atom, helper1, helper2
    ### (possibly helper3) with respect to the first atom
    ### (e.g. 0 is index of C1, N4 is 3 atoms away from C1, etc).
    ###
    dic_lipids_with_indexes = make_dic_lipids_with_indexes(universe_woH,
                                                           dic_lipid)
    # Get lipid name.
    resname = dic_lipid["resname"]
    # Select first residue of that lipid.
    selection = "resname {}".format(resname)
    first_lipid_residue = universe_woH.select_atoms(selection).residues[0]
    # Get name of 1st atom of that lipid.
    first_atom_name = first_lipid_residue.atoms[0].name
    ###
    ### 2) Now loop over the traj, residues and Catoms.
    ### At each iteration build Hs and calc OP.
    ###
    # Loop over frames (ts is a Timestep instance).
    for ts in universe_woH.trajectory:
        print("Dealing with frame {} at {} ps."
              .format(ts.frame, universe_woH.trajectory.time))
        if DEBUG:
            print("Looping now over residues...")
            print()
        # Loop over the 1st atom of each lipid, which is equiv to loop *over
        # residues* (first_lipid_atom is an Atom instance, lipid_ix is an int
        # that will be used for storing OPs in dic_OP).
        selection = "resname {} and name {}".format(resname, first_atom_name)
        for lipid_ix, first_lipid_atom in enumerate(universe_woH.select_atoms(selection)):
            if DEBUG:
                print("Dealing with Cname", first_lipid_atom)
                print("which is part of residue", first_lipid_atom.residue)
                print("Now looping over atoms of this residue")
                print()
            # Get the index of this first atom.
            ix_first_atom_res = first_lipid_atom.ix
            # Now loop over each carbon on which we want to build Hs
            # (Cname is a string).
            for Cname in dic_lipids_with_indexes.keys():
                # Get Cname coords.
                if len(dic_lipids_with_indexes[Cname]) == 6:
                    _, _, _, Cname_ix, helper1_ix, helper2_ix = dic_lipids_with_indexes[Cname]
                else:
                    _, _, _, _, Cname_ix, helper1_ix, helper2_ix, helper3_ix = dic_lipids_with_indexes[Cname]
                Cname_position = ts[Cname_ix+ix_first_atom_res]
                if DEBUG:
                    print("Dealing with Cname", Cname)
                    sel = first_lipid_atom.residue.atoms.select_atoms
                    Cname_atom = sel("name {}".format(Cname))[0]
                    print(Cname_atom, Cname_atom.position)
                    if len(dic_lipid[Cname]) == 3:
                        _, helper1_name, helper2_name = dic_lipid[Cname]
                    else:
                        _, helper1_name, helper2_name, helper3_name = dic_lipid[Cname]
                    helper1_atom = sel("name {}".format(helper1_name))[0]
                    print("helper1", helper1_atom, helper1_atom.position)
                    helper2_atom = sel("name {}".format(helper2_name))[0]
                    print("helper2", helper2_atom, helper2_atom.position)
                    if len(dic_lipid[Cname]) == 4:
                        helper3_atom = sel("name {}".format(helper3_name))[0]
                        print("helper3", helper3_atom, helper3_atom.position)
                # Get newly built H(s) on that atom.
                Hs_coor = fast_buildHs_on_1C(dic_lipids_with_indexes, ts,
                                             Cname, ix_first_atom_res)
                if DEBUG:
                    print("Cname_position with fast indexing:", Cname_position)
                    print("helper1_position with fast indexing:",
                          ts[helper1_ix+ix_first_atom_res])
                    print("helper2_position with fast indexing:",
                          ts[helper2_ix+ix_first_atom_res])
                    if len(dic_lipid[Cname]) == 4:
                        print("helper3_position with fast indexing:",
                              ts[helper3_ix+ix_first_atom_res])
                # To retrieve Hname, we need a counter.
                counter4Hname = 0
                # Loop over all Hs.
                for i, H_coor in enumerate(Hs_coor):
                    # Retrieve name of newly built H.
                    Hname = dic_Cname2Hnames[Cname][counter4Hname]
                    # Calc and store OP for that couple C-H.
                    Cname_position = ts[Cname_ix+ix_first_atom_res]
                    op = calc_OP(Cname_position, H_coor)
                    # Old way: dic_OP[(Cname, Hname)].append(op)
                    if (Cname, Hname) in dic_OP:
                        dic_OP[(Cname, Hname)][lipid_ix].append(op)
                    if DEBUG:
                        print(Hname, H_coor, "OP:", op)
                    # Increment counter4Hname for retrieving next H.
                    counter4Hname += 1
                if DEBUG:
                    print() ; print()
    if DEBUG:
        print("Final dic_OP:", dic_OP)
        print()


def make_dic_atname2genericname(filename):
    """Make a dict of correspondance between generic H names and PDB names.

    This dict will look like the following: {('C1', 'H11'): 'gamma1_1', ...}.
    Useful for outputing OP with generic names (such as beta1, beta 2, etc.).
    Such files can be found on the NMRlipids MATCH repository:
    https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs.

    Parameters
    ----------
    filename : str
        Filename containing OP definition
        (e.g. `order_parameter_definitions_MODEL_Berger_POPC.def`).

    Returns
    -------
    Ordered dictionnary
        Keys are tuples of (C, H) name, values generic name (as described
        above in this docstring). The use of an ordered dictionnary ensures
        we get always the same order in the output OP.
    """
    dic = collections.OrderedDict()
    try:
        with open(filename, "r") as f:
            for line in f:
                # TODO: This line might have to be changed if the file contains more than
                # 4 columns.
                name, _, C, H = line.split()
                dic[(C, H)] = name
    except:
        raise UserWarning("Can't read order parameter definition in "
                          "file {}".format(filename))
    return dic

def init_dic_OP(universe_woH, dic_atname2genericname):
    """TODO Complete docstring.
    """
    ### To calculate the error, we need to first average over the
    ### trajectory, then over residues.
    ### Thus in dic_OP, we want for each key a list of lists, for example:
    ### OrderedDict([
    ###              (('C1', 'H11'), [[], [], ..., [], []]),
    ###              (('C1', 'H12'), [[], ..., []]),
    ###              ...
    ###              ])
    ### Thus each sublist will contain OPs for one residue.
    ### e.g. ('C1', 'H11'), [[OP res 1 frame1, OP res1 frame2, ...],
    ###                      [OP res 2 frame1, OP res2 frame2, ...], ...]
    dic_OP = collections.OrderedDict()
    # We also need the correspondance between residue number (resnum) and
    # its index in dic_OP.
    dic_corresp_numres_index_dic_OP = {}
    # Create these sublists by looping over each lipid.
    for key in dic_atname2genericname:
        dic_OP[key] = []
        # Get lipid name.
        resname = dic_lipid["resname"]
        selection = "resname {}".format(resname)
        # Loop over each residue on which we want to calculate the OP on.
        for i, residue in enumerate(universe_woH.select_atoms(selection).residues):
            dic_OP[key].append([])
            dic_corresp_numres_index_dic_OP[residue.resid] = i
    if DEBUG:
        print("Initial dic_OP:", dic_OP)
        print("dic_corresp_numres_index_dic_OP:", dic_corresp_numres_index_dic_OP)
    return dic_OP, dic_corresp_numres_index_dic_OP


def make_dic_Cname2Hnames(dic_OP):
    """TODO Complete Docstring.
    """
    dic = {}
    for Cname, Hname in dic_OP.keys():
        if Cname not in dic:
            dic[Cname] = (Hname,)
        else:
            dic[Cname] += (Hname,)
    if DEBUG:
        print("dic_Cname2Hnames contains:", dic)
    return dic


if __name__ == "__main__":
    # 0) Fist ensure Python 3 is used!!!
    major, minor, _, _, _ = sys.version_info
    if major != 3:
        raise UserWarning("buildH only works with Python 3.")
    else:
       if minor < 6:
           warnings.warn("Python version >= 3.6 is recommended with buildH.", UserWarning)
       else:
           print("Python version OK!")

    # 1) Parse arguments.
    # TODO --> Make a function for that.
    message="""This program builds hydrogens and calculate the order
    parameters
    (OP) from a united-atom trajectory. If -opx is requested, pdb and xtc
    output files with hydrogens are created but OP calculation will be slow.
    If no trajectory output is requested (no use of flag -opx), it uses a
    fast procedure to build hydrogens and calculate the OP.
    """
    parser = argparse.ArgumentParser(description=message)
    # Avoid tpr for topology cause there's no .coord there!
    parser.add_argument("topfile", type=isfile,
                        help="Topology file (pdb or gro).")
    parser.add_argument("-x", "--xtc", type=isfile, 
                        help="Input trajectory file in xtc format.")
    parser.add_argument("-l", "--lipid", type=str, required=True,
                        help="Residue name of lipid to calculate the OP on (e.g. POPC).")
    parser.add_argument("-d", "--defop", required=True, type=isfile, 
                        help="Order parameter definition file. Can be found on "
                        "NMRlipids MATCH repository:"
                        "https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs")
    parser.add_argument("-opx", "--opdbxtc", help="Base name for trajectory "
                        "output with hydrogens. File extension will be "
                        "automatically added. For example -opx trajH will "
                        "generate trajH.pdb and trajH.xtc. "
                        "So far only xtc is supported.")
    parser.add_argument("-o", "--out", help="Output base name for storing "
                        "order parameters. Extention \".out\" will be "
                        "automatically added. Default name is OP_buildH.out.",
                        default="OP_buildH.out")
    args = parser.parse_args()

    # Top file is "args.topfile", xtc file is "args.xtc", pdb output file is
    # "args.pdbout", xtc output file is "args.xtcout".
    # Check topology file extension.
    if not args.topfile.endswith("pdb") and not args.topfile.endswith("gro"):
        parser.error("Topology must be given in pdb or gro format")
    # Check residue name validity.
    # Get the dictionnary with helper info using residue name (args.lipid
    # argument). Beware, this dict is then called `dic_lipid` *without s*,
    # while `dic_lipids.py` (with an s) is a module with many different dicts
    # (of different lipids) the user can choose.
    try:
        dic_lipid = getattr(dic_lipids, args.lipid)
    except:
        parser.error("Lipid dictionnary {} doesn't exist in dic_lipids.py".format(args.lipid))

    # 2) Create universe without H.
    print("Constructing the system...")
    if args.xtc:
        try:
            universe_woH = mda.Universe(args.topfile, args.xtc)
        except:
            raise UserWarning("Can't create MDAnalysis universe with files {} "
                              "and {}".format(args.topfile, args.xtc))
    else:
        try:
            universe_woH = mda.Universe(args.topfile)
        except:
            raise UserWarning("Can't create MDAnalysis universe with file {}"
                              .format(args.topfile))
    print("System has {} atoms".format(len(universe_woH.coord)))

    # 2) Initialize dic for storing OP.
    # Init dic of correspondance : {('C1', 'H11'): 'gamma1_1',
    # {('C1', 'H11'): 'gamma1_1', ...}.
    dic_atname2genericname = make_dic_atname2genericname(args.defop)
    # Initialize dic_OP (see function init_dic_OP() for the format).
    dic_OP, dic_corresp_numres_index_dic_OP = init_dic_OP(universe_woH,
                                                          dic_atname2genericname)
    # Initialize dic_Cname2Hnames.
    dic_Cname2Hnames = make_dic_Cname2Hnames(dic_OP)

    # If traj output files are requested.
    # NOTE Here, we need to reconstruct all Hs. Thus the op definition file (passed
    #  with arg -d) needs to contain all possible C-H pairs !!!
    if args.opdbxtc:
        #3) Prepare the system.
        # First check that dic_OP contains all possible C-H pairs.
        # NOTE The user has to take care that .def file has the right atom names !!!
        for atname in dic_lipid.keys():
            if atname != "resname":
                # Check if carbon is present in the definition file.
                if atname not in dic_Cname2Hnames:
                    print("Error: When -opx option is used, the order param "
                          "definition file (passed with -d arg) must contain "
                          "all possible carbons on which we want to rebuild "
                          "hydrogens.")
                    print("Found:", list(dic_Cname2Hnames.keys()))
                    print("Needs:", list(dic_lipid.keys()))
                    raise UserWarning("Order param definition file incomplete.")
                # Check that the 3 Hs are in there for that C.
                nbHs_in_def_file = len(dic_Cname2Hnames[atname])
                tmp_dic = {"CH": 1, "CHdoublebond": 1, "CH2": 2, "CH3": 3}
                correct_nb_of_Hs = tmp_dic[dic_lipid[atname][0]]
                if  correct_nb_of_Hs != nbHs_in_def_file:
                    print("Error: When -opx option is used, the order param "
                          "definition file (passed with -d arg) must contain "
                          "all possible C-H pairs to rebuild.")
                    print("Expected {} hydrogen(s) to rebuild for carbon {}, "
                          "got {} in definition file {}."
                          .format(correct_nb_of_Hs, atname,
                                  dic_Cname2Hnames[atname], args.defop))
                    raise UserWarning("Wrong number of Hs to rebuild.")
        # Create filenames.
        pdbout_filename = args.opdbxtc + ".pdb"
        xtcout_filename = args.opdbxtc + ".xtc"
        # Build a new universe with H.
        # Build a pandas df with H.
        new_df_atoms = build_all_Hs_calc_OP(universe_woH, dic_lipid,
                                            dic_Cname2Hnames,
                                            return_coors=True)
        # Create a new universe with H using that df.
        print("Writing new pdb with hydrogens.")
        # Write pdb with H to disk.
        with open(pdbout_filename, "w") as f:
            f.write(pandasdf2pdb(new_df_atoms))
        # Then create the universe with H from that pdb.
        universe_wH = mda.Universe(pdbout_filename)
        # Create an xtc writer.
        print("Writing trajectory with hydrogens in xtc file.")
        newxtc = XTC.XTCWriter(xtcout_filename, len(universe_wH.atoms))
        # Write 1st frame.
        newxtc.write(universe_wH)

        # 4) Loop over all frames of the traj *without* H, build Hs and
        # calc OP (ts is a Timestep instance).
        for ts in universe_woH.trajectory:
            print("Dealing with frame {} at {} ps."
                  .format(ts.frame, universe_woH.trajectory.time))
            # Build H and update their positions in the universe *with* H (in place).
            # Calculate OPs on the fly while building Hs  (dic_OP changed in place).
            build_all_Hs_calc_OP(universe_woH, dic_lipid, dic_Cname2Hnames,
                                universe_wH=universe_wH, dic_OP=dic_OP,
                                dic_corresp_numres_index_dic_OP=dic_corresp_numres_index_dic_OP)
            # Write new frame to xtc.
            newxtc.write(universe_wH)
        # Close xtc.
        newxtc.close()

    # 6) If no traj output file requested, use fast indexing to speed up OP
    # calculation. The function fast_build_all_Hs() returns nothing, dic_OP
    # is modified in place.
    if not args.opdbxtc:
        fast_build_all_Hs_calc_OP(universe_woH, dic_OP, dic_lipid, dic_Cname2Hnames)

    # 7) Output results.
    # Pickle results? (migth be useful in the future)
    # TODO Implement that option.
    if PICKLE:
        with open("OP.pickle", "wb") as f:
            # Pickle the dic using the highest protocol available.
            pickle.dump(dic_OP, f, pickle.HIGHEST_PROTOCOL)
        #  To unpickle
        #with open("OP.pickle", "rb") as f:
        #    dic_OP = pickle.load(f)
    # Output to a file.
    with open("{}.jmelcr_style.out".format(args.out), "w") as f, \
        open("{}.apineiro_style.out".format(args.out), "w") as f2:
        # J. Melcr output style.
        f.write("# {:18s} {:7s} {:5s} {:5s}  {:7s} {:7s} {:7s}\n"
                .format("OP_name", "resname", "atom1", "atom2", "OP_mean",
                "OP_stddev", "OP_stem"))
        f.write("#-------------------------------"
                "-------------------------------------\n")
        # Loop over each pair (C, H).
        for Cname, Hname in dic_atname2genericname.keys():
            name = dic_atname2genericname[(Cname, Hname)]
            if DEBUG:
                print("Pair ({}, {}):".format(Cname, Hname))
            # Cast list of lists to a 2D-array. It should have dimensions
            # (nb_lipids, nb_frames).
            ### Thus each sublist will contain OPs for one residue.
            ### e.g. ('C1', 'H11'), [[OP res 1 frame1, OP res1 frame2, ...],
            ###                      [OP res 2 frame1, OP res2 frame2, ...],
            ####                     ...]
            a = np.array(dic_OP[(Cname, Hname)])
            if DEBUG:
                print("Final OP array has shape (nb_lipids, nb_frames):",
                      a.shape)
                print()
            # General mean over lipids and over frames (for that (C, H) pair).
            mean = np.mean(a)
            # Average over frames for each (C, H) pair.  Because of how the
            # array is organized (see above), we need to average horizontally
            # (i.e. using axis=1).
            # means is a 1D-array with nb_lipids elements.
            means = np.mean(a, axis=1)
            # Calc standard deviation and STEM (std error of the mean).
            std_dev = np.std(means)
            stem = np.std(means) / np.sqrt(len(means))
            f.write("{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f}\n"
                    .format(name, dic_lipid["resname"], Cname, Hname, mean,
                            std_dev, stem))
        # A. Pineiro output style.
        f2.write("Atom_name  Hydrogen\tOP\t      STD\t   STDmean\n")
        list_unique_Cnames = []
        for Cname, Hname in dic_OP.keys():
            if Cname not in list_unique_Cnames:
                list_unique_Cnames.append(Cname)
        # Order of carbons is similar to that in the PDB.
        list_unique_Cnames_ordered = []
        selection = "resname {}".format(dic_lipid["resname"])
        for atom in universe_woH.select_atoms(selection).residues[0].atoms:
            if atom.name in list_unique_Cnames:
                list_unique_Cnames_ordered.append(atom.name)
        # Now write output.
        for Cname in list_unique_Cnames_ordered:
            cumulative_list_for_that_carbon = []
            for i, Hname in enumerate([H for C, H in dic_OP.keys() if C == Cname]):
                cumulative_list_for_that_carbon += dic_OP[Cname, Hname]
                a = np.array(dic_OP[Cname, Hname])
                mean = np.mean(a)
                means = np.mean(a, axis=1)
                std_dev = np.std(means)
                stem = np.std(means) / np.sqrt(len(means))
                if i == 0:
                    f2.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                             .format(Cname, "HR", mean, std_dev, stem))
                elif i == 1:
                    f2.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                             .format("", "HS", mean, std_dev, stem))
                elif i == 2:
                    f2.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                             .format("", "HT", mean, std_dev, stem))
            a = np.array(cumulative_list_for_that_carbon)
            mean = np.mean(a)
            means = np.mean(a, axis=1)
            std_dev = np.std(means)
            stem = np.std(means) / np.sqrt(len(means))
            f2.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n\n"
                     .format("", "AVG", mean, std_dev, stem))

    print("Results written to {}".format(args.out))

