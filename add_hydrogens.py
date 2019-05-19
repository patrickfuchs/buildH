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
import numpy as np
import pandas as pd
import MDAnalysis as mda

import dic_lipids

# Constants.
LENGTH_CH_BOND = 1.0 # in Angst
# From https://en.wikipedia.org/wiki/Tetrahedron, tetrahedral angle equals
# arccos(-1/3) ~ 1.9106 rad or 109.47 deg.
TETRAHEDRAL_ANGLE = np.arccos(-1/3)

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
    return vec / magnitude(vec)


def magnitude(vec):
    """Returns the magnitude of a vector.

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
    costheta = np.dot(vec1,vec2)/(magnitude(vec1)*magnitude(vec2))
    if costheta > 1.0 or costheta < -1.0:
        raise(ValueError, "Cosine cannot be larger than 1.0 or less than -1.0")
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
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                rows.append((atnum, atname, resname, resnum, x, y, z))
    df_atoms = pd.DataFrame(rows, columns=["atnum", "atname", "resname",
                                           "resnum", "x", "y", "z"])
    return df_atoms


def pdb2list_pandasdf_residues(pdb_filename):
    """Reads a PDB file and returns a list of pandas dataframes for each residue.

    Arguments
    ---------
    pdb_filename : string

    Returns
    -------
    list of pandas dataframe representing a residue
        Each dataframe has the folowing columns: atnum, atname, resname, resnum,
        x, y, z, typeofH2build, helper1_name, helper2_name.
    """
    # Put PDB rows in a list.
    rows = []
    # This list will be used for indexing rows with atom names.
    all_atom_names = []
    with open(pdb_filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atnum = int(line[6:11])
                atname = line[12:16].strip()
                all_atom_names.append(atname)
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if atname in dic_lipids.POPC:
                    typeofH2build, helper1_name, helper2_name = dic_lipids.POPC[atname]
                else:
                    typeofH2build, helper1_name, helper2_name = None, None, None
                rows.append((atnum, atname, resname, resnum, x, y, z,
                             typeofH2build, helper1_name, helper2_name))
    # Make a first dataframe with those rows.
    df_atoms = pd.DataFrame(rows, index=all_atom_names,
                            columns=["atnum", "atname", "resname",
                                     "resnum", "x", "y", "z",
                                     "typeofH2build", "helper1_name",
                                     "helper2_name"])
    # Make a list of dataframes (each df is a residue).
    list_df_residues = []
    for res_num in df_atoms.resnum.unique():
        list_df_residues.append( df_atoms[ df_atoms["resnum"] == res_num ] )
    return list_df_residues


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
    for i, row_atom in df.iterrows():
        atnum, atname, resname, resnum, x, y, z = row_atom
        atnum = int(atnum)
        resnum = int(resnum)
        s += ("{:6s}{:5d} {:>4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}"
              "{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
              .format("ATOM", atnum, atname, "", resname, chain, resnum, "",
                      x, y, z, 1.0, 0.0, "", ""))
    return s
 
    
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
    v4 = normalize(np.cross(v3, v2))
    # Rotation axis is atom->helper1 vec minus atom->helper2 vec.
    rotation_axis = normalize(v2 - v3)
    # Vector to be rotated by theta/2, perpendicular to rotation axis and v4.
    vec_to_rotate = normalize(np.cross(v4, rotation_axis))
    # Reconstruct the two hydrogens.
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
    rotation_axis = normalize(np.cross(v2, v3))
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
    rotation_axis = normalize(np.cross(v3, v2))
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


def get_name_H(name_carbon, nb_of_H):
    """Returns the name of newly built hydrogens.

    Parameters
    ----------
    name_carbon : string
        The name of the carbon atom holding the hydrogen(s).
    nb_of_H : int
        Number of names to generate.

    Returns
    -------
    str
        The name of the unique H if nb_of_H equals 1.
    or tuple of str
        The names of the 2 or 3 rebuilt H if nb_of_H > 1.
    """
    name_H1 = name_carbon.replace("C", "H") + "1"
    name_H2 = name_carbon.replace("C", "H") + "2"
    name_H3 = name_carbon.replace("C", "H") + "3"
    if nb_of_H == 1:
        return name_H1
    elif nb_of_H == 2:
        return name_H1, name_H2
    elif nb_of_H == 3:
        return name_H1, name_H2, name_H3


def buildH_wMDanalysis(pdb_filename, return_coors=False):
    """Builds hydrogens from a united atom frame.
    
    TODO This function gets big, divide it in 2: actual func reads the traj 
         frame by frame and then calls another one which builds H.
    TODO2 Implement same stuff with a trajectory instead of single PDB frame.
    TODO3 Implement order parameter calculation.
    TODO4 BLABLABLA.
    """
    # load PDB
    universe = mda.Universe(pdb_filename)
    if return_coors:
        # The list newrows will be used to store the new molecule *with* H.
        newrows = []
        # Counter for numbering the new mlcs with H.
        new_atom_num = 1
    # Loop over all atoms.
    for atom in universe.atoms:
        if return_coors:
            resnum = atom.resnum
            # beware, resname must be 3 letters long in my routine
            resname = atom.resname[:-1]
            name = atom.name
            # Append atom to the new list.
            # 0      1       2        3       4  5  6
            # atnum, atname, resname, resnum, x, y, z
            newrows.append([new_atom_num, name, resname, resnum]
                           + list(atom.position))
            new_atom_num += 1
        if atom.name in dic_lipids.POPC:
            typeofH2build = dic_lipids.POPC[atom.name][0]
            if typeofH2build == "CH2":
                _, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
                # atom is a Atom object.
                # atom.residue.atoms is a list of atoms we can select with
                # method .select_atoms().
                # To avoid too long line, we shorten its name to `sel`.
                sel = atom.residue.atoms.select_atoms
                # [0] is because select_atoms returns a AtomGroup which
                # contains only 1 atom.
                helper1_coor = sel("name {0}".format(helper1_name))[0].position
                helper2_coor = sel("name {0}".format(helper2_name))[0].position
                # Build H(s).
                H1_coor, H2_coor = get_CH2(atom.position, helper1_coor,
                                           helper2_coor)
                ####
                #### We could calculate here the order param on the fly :-D !
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name, H2_name = get_name_H(atom.name, 2)
                    newrows.append([new_atom_num, H1_name, resname, resnum]
                                   + list(H1_coor))
                    new_atom_num += 1
                    newrows.append([new_atom_num, H2_name, resname, resnum]
                                   + list(H2_coor) )
                    new_atom_num += 1
            elif typeofH2build == "CH":
                _, helper1_name, helper2_name, helper3_name = (dic_lipids
                                                               .POPC[atom.name])
                sel = atom.residue.atoms.select_atoms
                helper1_coor = sel("name {0}".format(helper1_name))[0].position
                helper2_coor = sel("name {0}".format(helper2_name))[0].position
                helper3_coor = sel("name {0}".format(helper3_name))[0].position
                H1_coor = get_CH(atom.position, helper1_coor, helper2_coor,
                                 helper3_coor)
                ####
                #### We could calculate here the order param on the fly :-D !
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name = get_name_H(atom.name, 1)
                    newrows.append([new_atom_num, H1_name, resname, resnum]
                                   + list(H1_coor))
                    new_atom_num += 1
            elif typeofH2build == "CHdoublebond":
                _, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
                sel = atom.residue.atoms.select_atoms
                helper1_coor = sel("name {0}".format(helper1_name))[0].position
                helper2_coor = sel("name {0}".format(helper2_name))[0].position
                H1_coor = get_CH_double_bond(atom.position, helper1_coor,
                                             helper2_coor)
                ####
                #### We could calculate here the order param on the fly :-D !
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name = get_name_H(atom.name, 1)
                    newrows.append([new_atom_num, H1_name, resname, resnum]
                                   + list(H1_coor))
                    new_atom_num += 1
            elif typeofH2build == "CH3":
                _, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
                sel = atom.residue.atoms.select_atoms
                helper1_coor = sel("name {0}".format(helper1_name))[0].position
                helper2_coor = sel("name {0}".format(helper2_name))[0].position
                H1_coor, H2_coor, H3_coor = get_CH3(atom.position,
                                                    helper1_coor, helper2_coor)
                ####
                #### We could calculate here the order param on the fly :-D !
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name, H2_name, H3_name = get_name_H(atom.name, 3)
                    newrows.append([new_atom_num, H1_name, resname, resnum]
                                   + list(H1_coor))
                    new_atom_num += 1
                    newrows.append([new_atom_num, H2_name, resname, resnum]
                                   + list(H2_coor))
                    new_atom_num += 1
                    newrows.append([new_atom_num, H3_name, resname, resnum]
                                   + list(H3_coor))
                    new_atom_num += 1
    if return_coors:
        # Create a dataframe to store the mlc with added hydrogens.
        new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname",
                                                      "resname", "resnum",
                                                      "x", "y", "z"])
        return new_df_atoms


if __name__ == "__main__":
    use_pandas = False
    pdb_filename = "1POPC.pdb" #"POPC_only.pdb"
    if use_pandas:
        # read coordinates in a pandas dataframe
        list_df_residues = pdb2list_pandasdf_residues(pdb_filename)
        new_df_atoms = reconstruct_hydrogens3(list_df_residues)
        print(pandasdf2pdb(new_df_atoms))
    else:
        new_df_atoms = buildH_wMDanalysis(pdb_filename, return_coors=True)
        print(pandasdf2pdb(new_df_atoms))
