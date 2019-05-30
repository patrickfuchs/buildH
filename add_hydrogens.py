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
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.coordinates.XTC as XTC

import dic_lipids

# Constants.
LENGTH_CH_BOND = 1.0 # in Angst
# From https://en.wikipedia.org/wiki/Tetrahedron, tetrahedral angle equals
# arccos(-1/3) ~ 1.9106 rad ~ 109.47 deg.
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


def buildH_on_1C(atom):
    """Reconstructs 1, 2 or 3 H on a given carbon.

    This function is a wrapper which gathers the coordinates of the helpers and
    call the function that build 1, 2 or 3 H.

    The name of the helpers as well as the type of H to build are described in
    a dictionnary stored in dic_lipids.py.

    Parameters
    ----------
    atom : MDAnalysis Atom instance
        (see https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=atom%20class#MDAnalysis.core.groups.Atom
        for class definition)

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3 
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if there's
        only one H that has been rebuilt.
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipids.POPC[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
    else:
        (typeofH2build, helper1_name, helper2_name,
         helper3_name) = dic_lipids.POPC[atom.name]
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
        H1_coor = get_CH_double_bond(atom.position, helper1_coor, helper2_coor)
        return (H1_coor,)
    elif typeofH2build == "CH3":
        H1_coor, H2_coor, H3_coor = get_CH3(atom.position,
                                            helper1_coor, helper2_coor)
        return (H1_coor, H2_coor, H3_coor)
    else:
        raise UserWarning("Wrong code for typeofH2build, expect 'CH2', 'CH', "
                          "'CHdoublebond' or 'CH3', got {}"
                          .format(typeofH2build))


def build_all_H(universe_woH, universe_wH=None, return_coors=False):
    """Main function that reconstructs hyddrogens.

    This function shall be used in two modes :

    1) The first time this function is called, we have to construct a new 
    universe with hydrogens. One shall call it like this :

    new_data_frame = build_all_H(universe_woH, universe_wH=None, return_coors=False)

    This dataframe will be used to write a pdb with H, which will allow to build
    a new universe with H.

    2) For all the other frames, we just need to update the coordinates in the 
    universe *with* hydrogens. One shall call it like this :

    build_all_H(universe_woH, universe_wH=universe_wH)

    The function returns nothing and the coordinates of the universe *with* H
    are changed in place.

    Parameters
    ----------
    universe_woH : MDAnalysis universe
        This is the universe *without* hydrogen.
    universe_wH : MDAnalysis universe (optional)
        This is the universe *with* hydrogens.
    return_coors : boolean (optional)
        If True, the function will return a pandas dataframe containing the system *with* hydrogens.

    Returns
    -------
    pandas dataframe (optional)
        If parameter return_coors is True, this dataframe  containing the 
        system *with* hydrogens is returned.
    None
        If parameter return_coors is False, the function returns nothing.
    """
    if universe_wH:
        # We will need the index in the numpy array for updating coordinates
        # in the universe with H.
        atom_index_in_nparray = 0
    if return_coors:
        # The list newrows will be used to store the new molecule *with* H.
        newrows = []
        # Counter for numbering the new mlcs with H.
        new_atom_num = 1
    # Loop over all atoms.
    for atom in universe_woH.atoms:
        if universe_wH:
            # Update the position of the current atom in the universe with H.
            universe_wH.coord.positions[atom_index_in_nparray, :] = atom.position
            atom_index_in_nparray += 1
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
        # Build new H(s)?
        if atom.name in dic_lipids.POPC:
            # Build Hs and store them in a list of numpy 1D-arrays Hs_coor.
            # The "s" in Hs_coor means there can be more than 1 H:
            # For CH2, Hs_coor will contain: [H1_coor, H2_coor].
            # For CH3, Hs_coor will contain: [H1_coor, H2_coor, H3_coor].
            # For CH, Hs_coor will contain: [H1_coor].
            # For CHdoublebond, Hs_coor will contain: [H1_coor].
            Hs_coor = buildH_on_1C(atom)
            # Loop over Hs_coor (H_coor is a 1D-array with the 3 coors of 1 H).
            for i, H_coor in enumerate(Hs_coor):
                # Give a name to newly built H
                # (e.g. if C18 has 3 H, their name will be H181, H182 & H183).
                H_name = atom.name.replace("C", "H") + str(i+1)
                ####
                #### We calculate here the order param on the fly :-D !
                ####
                if return_coors:
                    # Add them to newrows.
                    newrows.append([new_atom_num, H_name, resname, resnum]
                                   + list(H_coor))
                    new_atom_num += 1
                if universe_wH:
                    # Update the position of the current H in the universe with H.
                    universe_wH.coord.positions[atom_index_in_nparray, :] = H_coor
                    atom_index_in_nparray += 1
    if return_coors:
        # Create a dataframe to store the mlc with added hydrogens.
        new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname",
                                                      "resname", "resnum",
                                                      "x", "y", "z"])
        return new_df_atoms

    

if __name__ == "__main__":
    # Parse arguments.
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
    # Create universe.
    print("Constructing the system...")
    if args.xtc:
        universe_woH = mda.Universe(args.topfile, args.xtc)
    else:
        universe_woH = mda.Universe(args.topfile)
    print("System has {} atoms".format(len(universe_woH.coord)))
    # Build a pandas df with H.
    new_df_atoms = build_all_H(universe_woH, return_coors=True)
    # Create a new universe with H.
    if args.pdbout:
        print("Writing first frame.")
        # Write pdb with H to disk.
        with open(args.pdbout, "w") as f:
            f.write(pandasdf2pdb(new_df_atoms))
        # Then create the universe with H from that pdb.
        universe_wH = mda.Universe(args.pdbout)
    else:
        exit("For now please specify --pdbout argument.")
        ###
        ### !!!FIX ME !!! So far when using StringIO stream, it complains.
        ###
        # We don't want to create a pdb file, use a stream instead.
        #pdb_s = pandasdf2pdb(new_df_atoms)
        #universe_wH = mda.Universe(io.StringIO(pdb_s), format="pdb")
    if args.xtcout:
        # Create an xtc writer.
        newxtc = XTC.XTCWriter(args.xtcout, len(universe_wH.atoms))
        # Write 1st frame.
        newxtc.write(universe_wH)
    # Loop over all frames of the traj *without* H.
    # (ts is a timestep object).
    for ts in universe_woH.trajectory:
        print("Dealing with frame {} at {} ps."
              .format(ts.frame, universe_woH.trajectory.time))
        # Build H and update positions in the universe *with* H.
        build_all_H(universe_woH, universe_wH=universe_wH)
        if args.xtcout:
            # Write new frame to xtc.
            newxtc.write(universe_wH)
    if args.xtcout:
        # Close xtc.
        newxtc.close()
