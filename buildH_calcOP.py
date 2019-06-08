#!/usr/bin/env python3
# coding: utf-8

"""
This script reconstructs hydrogens from a united-atom trajectory.

BLABLABLA TODO

This code is inspired from that of Jon Kapla originally written in fortran
(https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242).

Note, that all coordinates in this script are handled using numpy 1D-arrays
of 3 elements, e.g. atom_coor = np.array((x, y, z)).
"""

__authors__ = ("Patrick Fuchs", "Amélie Bâcle", "Hubert Santuz",
               "Pierre Poulain")
__contact__ = ("patrickfuchs", "abacle", "hublot", "pierrepo") # on github
__version__ = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2019/05"

# Modules.
import argparse
import io
import pickle
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.coordinates.XTC as XTC

import dic_lipids
import copy

# Constants.
# From https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond
LENGTH_CH_BOND = 1.09 # in Angst
# From https://en.wikipedia.org/wiki/Tetrahedron, tetrahedral angle equals
# arccos(-1/3) ~ 1.9106 rad ~ 109.47 deg.
TETRAHEDRAL_ANGLE = np.arccos(-1/3)


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

    This function was initially written by @jmelcr.

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


def pdb2pandasdf(pdb_filename):
    """Reads a PDB file and returns a pandas data frame.

    Arguments
    ---------
    pdb_filename : string

    Returns
    -------
    pandas dataframe
        The col index are: atnum, atname, resname, resnum, x, y, z
    """
    rows = []
    with open(pdb_filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atnum = int(line[6:11])
                atname = line[12:16].strip()
                resname = line[17:21].strip()
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                rows.append((atnum, atname, resname, resnum, x, y, z))
    df_atoms = pd.DataFrame(rows, columns=["atnum", "atname", "resname",
                                           "resnum", "x", "y", "z"])
    return df_atoms


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
        # See for pdb format: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
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
    """Returns the cross product between vectors A & B
    see e.g. http://hyperphysics.phy-astr.gsu.edu/hbase/vvec.html
    """
    x = (A[1]*B[2]) - (A[2]*B[1])
    y = (A[0]*B[2]) - (A[2]*B[0])
    z = (A[0]*B[1]) - (A[1]*B[0])
    return np.array((x, -y, z))

##@profile   
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
    # TODO Rename norm_vec_H1 (this is confusing) --> maybe unit_vect_H1 ? (same for H2)
    norm_vec_H1 = apply_rotation(vec_to_rotate, rotation_axis,
                                 -TETRAHEDRAL_ANGLE/2)
    hcoor_H1 = LENGTH_CH_BOND * norm_vec_H1 + atom
    norm_vec_H2 = apply_rotation(vec_to_rotate, rotation_axis,
                                 TETRAHEDRAL_ANGLE/2)
    hcoor_H2 = LENGTH_CH_BOND * norm_vec_H2 + atom
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
    coor_H = LENGTH_CH_BOND * normalize(atom - v2) + atom
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
    # calc angle theta helper1-atom-helper2 (in rad).
    theta = calc_angle(helper1, atom, helper2)
    # atom->helper1 vector.
    v2 = helper1 - atom
    # atom->helper2 vector.
    v3 = helper2 - atom
    # The rotation axis is orthogonal to the atom/helpers plane.
    #rotation_axis = normalize(np.cross(v2, v3))
    rotation_axis = normalize(cross_product(v2, v3))
    # Reconstruct H by rotating v3 by theta.
    norm_vec_H = apply_rotation(v3, rotation_axis, theta)
    coor_H = LENGTH_CH_BOND * norm_vec_H + atom
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
    norm_vec_He = apply_rotation(v2, rotation_axis, theta)
    coor_He = LENGTH_CH_BOND * norm_vec_He + atom
    ### Build CH3r.
    theta = (2/3) * np.pi
    rotation_axis = normalize(helper1 - atom)
    v4 = normalize(coor_He - atom)
    # Now we rotate atom->He bond around atom->helper1 bond by 2pi/3.
    norm_vec_Hr = apply_rotation(v4, rotation_axis, theta)
    coor_Hr = LENGTH_CH_BOND * norm_vec_Hr + atom
    ### Build CH3s.
    theta = -(2/3) * np.pi
    rotation_axis = normalize(helper1 - atom)
    v5 = normalize(coor_He - atom)
    # Last we rotate atom->He bond around atom->helper1 bond by -2pi/3.
    norm_vec_Hs = apply_rotation(v5, rotation_axis, theta)
    coor_Hs = LENGTH_CH_BOND * norm_vec_Hs + atom
    return coor_He, coor_Hr, coor_Hs

#@profile
def buildHs_on_1C(atom):
    """Builds 1, 2 or 3 H on a given carbon.

    This function is a wrapper which gathers the coordinates of the helpers
    and call the function that builds 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described
    in a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    atom : MDAnalysis Atom instance

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if
        there's only one H that has been rebuilt.
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipids.POPC[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
    else:
        typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipids.POPC[atom.name]
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


def new_buildHs_on_1C(dic_lipids_with_indexes, ts, Cname, ix_first_atom_res):
    """TODO BLABLABLABLA
    Builds 1, 2 or 3 H on a given carbon.

    This function is a wrapper which gathers the coordinates of the helpers
    and call the function that builds 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described
    in a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    atom : MDAnalysis Atom instance

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
        typeofH2build, _, _, Cname_id, helper1_id, helper2_id = dic_lipids_with_indexes[Cname]
    else:
        typeofH2build, _, _, _, Cname_id, helper1_id, helper2_id, helper3_id = dic_lipids_with_indexes[Cname]

    Cname_position = ts[Cname_id+ix_first_atom_res]
    helper1_coor = ts[helper1_id+ix_first_atom_res]
    helper2_coor = ts[helper2_id+ix_first_atom_res]

    if typeofH2build == "CH2":
        H1_coor, H2_coor = get_CH2(Cname_position, helper1_coor, helper2_coor)
        return (H1_coor, H2_coor)
    elif typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_coor = ts[helper3_id+ix_first_atom_res]
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


def get_indexes(atom, universe_woH):

    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipids.POPC[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
    else:
        typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipids.POPC[atom.name]
    # Get helper coordinates using atom, which an instance from Atom class.
    # atom.residue.atoms is a list of atoms we can select with
    # method .select_atoms().
    # To avoid too long line, we shorten its name to `sel`.
    sel = atom.residue.atoms.select_atoms
    helper1_id = sel("name {0}".format(helper1_name))[0].ix
    helper2_id = sel("name {0}".format(helper2_name))[0].ix
    if typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_id = sel("name {0}".format(helper3_name))[0].ix
        return (helper1_id, helper2_id, helper3_id)
    else:
        return (helper1_id, helper2_id)


def build_all_Hs_calc_OP(universe_woH, universe_wH=None, dic_OP=None, return_coors=False):
    """Main function that builds all hydrogens from an MDAnalysis universe and calculate order parameters.

    This function shall be used in two modes :

    1) The first time this function is called, we have to construct a new
    universe with hydrogens. One shall call it like this :

    new_data_frame = build_all_Hs_calc_OP(universe_woH, return_coors=True)

    The boolean return_coors set to True indicates to the function to return
    a pandas dataframe. This latter will be used later to build a new
    universe with H.

    2) For all the other frames, we just need to update the coordinates in
    the universe *with* hydrogens. One shall call it like this :

    build_all_Hs_calc_OP(universe_woH, universe_wH=universe_wH, dic_OP=dic_OP)

    In this case, the function also calculates the order parameter and returns
    nothing. The coordinates of the universe *with* H are update in place.
    The order parameter is also added in place (within dic_OP dictionnary).

    NOTE: This function in mode 2 is slow, thus it shall be used when one wants
    to create a trajectory with H (such as xtc or whatever format).

    Parameters
    ----------
    universe_woH : MDAnalysis universe
        This is the universe *without* hydrogen.
    universe_wH : MDAnalysis universe (optional)
        This is the universe *with* hydrogens.
    dic_OP : dictionnary
        This dictionnary contains all the order parameters. It is structured
        like this: {("C1", "H11"): [val1, val2, ...], ("C1", "H12"): [...], ...}.
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
    # Loop over all atoms in the universe without H..
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
        if (atom.name in dic_lipids.POPC and
            atom.residue.resname == dic_lipids.POPC["resname"]):
            # Build Hs and store them in a list of numpy 1D-arrays Hs_coor.
            # The "s" in Hs_coor means there can be more than 1 H:
            # For CH2, Hs_coor will contain: [H1_coor, H2_coor].
            # For CH3, Hs_coor will contain: [H1_coor, H2_coor, H3_coor].
            # For CH, Hs_coor will contain: [H1_coor].
            # For CHdoublebond, Hs_coor will contain: [H1_coor].
            Hs_coor = buildHs_on_1C(atom)
            # Loop over Hs_coor (H_coor is a 1D-array with the 3 coors of 1 H).
            for i, H_coor in enumerate(Hs_coor):
                # Give a name to newly built H
                # (e.g. if C18 has 3 H, their name will be H181, H182 & H183).
                H_name = atom.name.replace("C", "H") + str(i+1)
                ####
                #### We calculate here the order param on the fly :-D !
                ####
                if dic_OP:
                    op = calc_OP(atom.position, H_coor)
                    dic_OP[(atom.name, H_name)].append(op)
                if return_coors:
                    # Add them to newrows.
                    newrows.append([new_atom_num, H_name, resname, resnum]
                                   + list(H_coor))
                    new_atom_num += 1
                if universe_wH:
                    # Update the position of the current H in the universe with H.
                    universe_wH.coord.positions[row_index_coor_array, :] = H_coor
                    row_index_coor_array += 1
    if return_coors:
        # Create a dataframe to store the mlc with added hydrogens.
        new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname",
                                                      "resname", "resnum",
                                                      "x", "y", "z"])
        return new_df_atoms


####
#### QUICK TEST
####
#@profile
def QUICKTEST_buildHs_on_1C(atom, universe_woH, typeofH2build, index_helper1=None, index_helper2=None, index_helper3=None):
    """BLABLOBLIBLU BLABLOBLIBLU BLABLOBLIBLU TODO
    Builds 1, 2 or 3 H on a given carbon.

    This function is a wrapper which gathers the coordinates of the helpers 
    and call the function that builds 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described
    in a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    atom : MDAnalysis Atom instance

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3 
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if
        there's only one H that has been rebuilt.
    """
    #print("Ds QUICKTEST")
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    #if len(dic_lipids.POPC[atom.name]) == 3:
    #    typeofH2build, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
    #else:
    #    typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipids.POPC[atom.name]
    # Get helper coordinates using atom, which an instance from Atom class.
    # atom.residue.atoms is a list of atoms we can select with
    # method .select_atoms().
    # To avoid too long line, we shorten its name to `sel`.
    #sel = atom.residue.atoms.select_atoms
    #helper1_coor = sel("name {0}".format(helper1_name))[0].position
    #helper2_coor = sel("name {0}".format(helper2_name))[0].position
    #if typeofH2build == "CH":
    #    print("Poueeeeeeeeeet")
    #    helper3_coor = sel("name {0}".format(helper3_name))[0].position
    #print("OLD :-(")
    #print(sel("name {0}".format(helper1_name))[0])
    #print(sel("name {0}".format(helper2_name))[0])
    # if typeofH2build == "CH":
    #     print(sel("name {0}".format(helper3_name))[0])
    #     print(atom.position, helper1_coor, helper2_coor, helper3_coor)
    # else:
    #     print(atom.position, helper1_coor, helper2_coor)
    # print("NEW :-D")
    helper1_coor = universe_woH.coord.positions[index_helper1]
    helper2_coor = universe_woH.coord.positions[index_helper2]
    if index_helper3:
        helper3_coor = universe_woH.coord.positions[index_helper3]
    # print(universe_woH.atoms[index_helper1])
    # print(universe_woH.atoms[index_helper2])
    # if index_helper3:
    #     print(universe_woH.atoms[index_helper3])
    #     print(atom.position, helper1_coor, helper2_coor, helper3_coor)
    #     print()
    # else:
    #     print(atom.position, helper1_coor, helper2_coor)
    #     print()
    if typeofH2build == "CH2":
        H1_coor, H2_coor = get_CH2(atom.position, helper1_coor, helper2_coor)
        return (H1_coor, H2_coor)
    elif typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        #helper3_coor = sel("name {0}".format(helper3_name))[0].position
        universe_woH.coord.positions[index_helper3]
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


# Quick try to make a function that only loops over carbons on which we want to
# build new H and calc OP.
def fast_build_all_Hs(universe_woH, dic_OP):
    """BLABLABLA
    """
    # Get lipid name.
    resname = dic_lipids.POPC["resname"]
    # Get resnum of the 1st lipid whose name is `resname`.
    resnum_1st_lipid = universe_woH.select_atoms("resname {}".format(resname)).residues[0].resnum
    # Get name of 1st atom of lipid.
    first_atom_name = universe_woH.residues[0].atoms[0].name
    # Now expand dic_lipids.
    # We want {..., 'C1': ('CH3', 'N4', 'C5', 0, 3, 4), ...} 
    # Where the 3 last int are the index of the atom, helper1, helper2 
    # (possibly helper3) with respect to the first atom 
    # (e.g. 0 is index of C1, N4 is 3 atoms away from C1, etc)
    dic_lipids_with_indexes = copy.deepcopy(dic_lipids.POPC)
    
    ###
    ### !!! TODO !!!
    ### Remove key resname in dic_lipids_with_indexes
    ###

    # #Retrieve once the index of the helpers
    # # Loop over each couple C-H.
    for Cname in dic_lipids.POPC.keys():
        if Cname != "resname":
            # Loop over residues for a given Cname atom.
            for Catom in universe_woH.select_atoms("resid {} and name {}".format(resnum_1st_lipid, Cname)):
                helper_ids = get_indexes(Catom, universe_woH)
                dic_lipids_with_indexes[Cname] = dic_lipids_with_indexes[Cname] + (Catom.ix,) + helper_ids

    #current_resid = -1
    #idx_first_atom_res = -1
    DEBUG = False
    # Loop over frames. ts is a Timestep instance.
    for ts in universe_woH.trajectory:
        print("Dealing with frame {} at {} ps."
              .format(ts.frame, universe_woH.trajectory.time))
        #Loop over the 1st atom of each lipid, which is equiv to loop over residues.
        for first_lipid_atom in universe_woH.select_atoms("resname {} and name {}".format(resname, first_atom_name)):
            if DEBUG:
                print("Dealing with Cname", first_lipid_atom)
                print("    residue is", first_lipid_atom.residue)
            # Get the index of this first atom.
            ix_first_atom_res = first_lipid_atom.ix
            # Now loop over each carbon on which we want to build H.
            for Cname in dic_lipids.POPC.keys():
                if Cname != "resname":
                    # Get Cname coords.
                    if len(dic_lipids_with_indexes[Cname]) == 6:
                        _, _, _, Cname_id, _, _ = dic_lipids_with_indexes[Cname]
                    else:
                        _, _, _, _, Cname_id, _, _, _ = dic_lipids_with_indexes[Cname]
                    Cname_position = ts[Cname_id+ix_first_atom_res]
                    if DEBUG:
                        print("Dealing with Cname", Cname)
                    # Get newly built H on that atom.
                    Hs_coor = new_buildHs_on_1C(dic_lipids_with_indexes, ts, Cname, ix_first_atom_res)
                    #Hs_coor = buildHs_on_1C(Catom)
                    # Loop over all Hs.
                    if DEBUG:
                        print("Cname_position:", Cname_position)
                    for i, H_coor in enumerate(Hs_coor):
                        # Give a name to newly built H
                        # (e.g. if C18 has 3 H, their name will be H181,H182 & H183).
                        H_name = Cname.replace("C", "H") + str(i+1)
                        Cname_position = ts[Cname_id+ix_first_atom_res]
                        op = calc_OP(Cname_position, H_coor)
                        dic_OP[(Cname, H_name)].append(op)
                        if DEBUG:
                            print(H_name, H_coor, op)
                    if DEBUG:
                        print() ; print()



# For now make a quick dic (to be removed later).
def quick_dic():
    dic = {}
    with open("order_parameter_definitions_MODEL_Berger_POPC.def", "r") as f:
        for line in f:
            name, _, C, H = line.split()
            dic[(C, H)] = name
    return dic


if __name__ == "__main__":
    # 1) Parse arguments.
    # TODO --> Make a function for that.
    parser = argparse.ArgumentParser(description="Reconstruct hydrogens and calculate order parameter from a united-atom trajectory.")
    # Avoid tpr for topology cause there's no .coord there!
    parser.add_argument("topfile", type=str, help="topology file (pdb or gro)")
    parser.add_argument("--xtc", help="input trajectory file in xtc format")
    parser.add_argument("--pdbout", help="output pdb file name with hydrogens "
                        "(default takes topology name + \"H\")")
    parser.add_argument("--xtcout", help="output xtc file name with hydrogens "
                        "(default takes topology name + \"H\")")
    args = parser.parse_args()
    # Top file is "args.topfile", xtc file is "args.xtc", pdb output file is
    # "args.pdbout", xtc output file is "args.xtcout".
    # Check topology file extension.
    if not args.topfile.endswith("pdb") and not args.topfile.endswith("gro"):
        raise argparse.ArgumentTypeError("Topology must be given in pdb"
                                         " or gro format")
    # Check other extensions.
    if args.pdbout:
        if not args.pdbout.endswith("pdb"):
            raise argparse.ArgumentTypeError("pdbout must have a pdb extension")
    if args.xtcout:
        if not args.xtcout.endswith("xtc"):
            raise argparse.ArgumentTypeError("xtcout must have an xtc extension")

    # 2) Create universe without H.
    print("Constructing the system...")
    if args.xtc:
        universe_woH = mda.Universe(args.topfile, args.xtc)
    else:
        universe_woH = mda.Universe(args.topfile)
    print("System has {} atoms".format(len(universe_woH.coord)))


    # QUICK TRY to build H only on carbons (thus lines below are
    # temporarilly not executed).
    QUICK = True
    if not QUICK:
        # 3) Build a new universe with H.
        # Build a pandas df with H.
        new_df_atoms = build_all_Hs_calc_OP(universe_woH, return_coors=True)
        # Create a new universe with H using that df.
        if args.pdbout:
            print("Writing new pdb with hydrogens.")
            # Write pdb with H to disk.
            with open(args.pdbout, "w") as f:
                f.write(pandasdf2pdb(new_df_atoms))
            # Then create the universe with H from that pdb.
            universe_wH = mda.Universe(args.pdbout)
        else:
            ###
            ### !!!FIX ME !!! So far when using StringIO stream, an exception is
            ### raised which is not neat.
            ### See https://github.com/MDAnalysis/mdanalysis/issues/2089.
            ###
            # In this else we don't want to create a pdb file, use a stream instead.
            pdb_s = pandasdf2pdb(new_df_atoms)
            universe_wH = mda.Universe(io.StringIO(pdb_s), format="pdb")
        if args.xtcout:
            # Create an xtc writer.
            newxtc = XTC.XTCWriter(args.xtcout, len(universe_wH.atoms))
            # Write 1st frame.
            newxtc.write(universe_wH)

    # 4) Initialize dic for storing OP.
    # Init dic of correspondance : {('C1', 'H11'): 'gamma1_1',
    # {('C1', 'H11'): 'gamma1_1', ...}.
    # TODO --> Add arguments for passing file name with OP definition.
    # TODO --> Make a class for storing all this stuff!
    dic_atname2genericname = quick_dic()
    dic_OP = {}
    for key in dic_atname2genericname:
        dic_OP[key] = []

    # Quick try to loop only over C on which we want to build Hs.
    if QUICK:
        fast_build_all_Hs(universe_woH, dic_OP)

    # 5) Loop over all frames of the traj *without* H, build H and calc OP.
    # (ts is a Timestep instance).
    if not QUICK:
        for ts in universe_woH.trajectory:
            print("Dealing with frame {} at {} ps."
                .format(ts.frame, universe_woH.trajectory.time))
            # Build H and update their positions in the universe *with* H (in place).
            build_all_Hs_calc_OP(universe_woH, universe_wH=universe_wH, dic_OP=dic_OP)
            if args.xtcout:
                # Write new frame to xtc.
                newxtc.write(universe_wH)
        if args.xtcout:
            # Close xtc.
            newxtc.close()

    # 6) Output results.
    # Pickle results?
    PICKLE = False
    if PICKLE:
        with open("OP.pickle", "wb") as f:
            # Pickle the dic using the highest protocol available.
            pickle.dump(dic_OP, f, pickle.HIGHEST_PROTOCOL)
        #  To unpickle
        #with open("OP.pickle", "rb") as f:
        #    dic_OP = pickle.load(f)
    # Output to a file.
    with open("OUT.buildH", "w") as f:
        f.write("# OP_name    resname    atom1    atom2    OP_mean   OP_stddev  OP_stem\n"
                "#--------------------------------------------------------------------\n")
        for key in dic_atname2genericname.keys():
            name = dic_atname2genericname[key]
            at1, at2 = key
            a = np.array(dic_OP[key])
            #print("{:15s} {:4s} {:4s} {:10.6f} +/- {:10.6f}".format(name, at1, at2, a.mean(), a.std()))
            f.write("{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f}\n"
                    .format(name, "POPC", at1, at2, a.mean(), a.std(), 0.0))

