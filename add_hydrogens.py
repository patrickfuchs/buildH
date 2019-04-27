#!/usr/bin/env python

import numpy as np
import pandas as pd

"""
This script reconstructs hydrogens from BLABLABLA...
TODO
"""

def normalize(vec):
    """Normalizes a vector.

    Parameters
    ----------
    vec : numpy 1D-array.

    Returns
    -------
    numpy 1D-array
        The normalized vector.
    """
    nvec = vec / np.sqrt(np.sum(vec**2))
    return nvec


def v2q(vec, theta):
    """Translates a 3-D vector and angle theta in a quaternion.

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
    x, y, z = np.sin(theta/2)*normalize(vec)
    q = np.array([w, x, y, z])
    return q


def rotational_matrix(quaternion):
    """Translates a quaternion to a rotational matrix.

    Parameters
    ----------
    quaternion : numpy 1D-array of 4 elements.

    Returns
    -------
    numpy 2D-array (dimension [3, 3])
        The rotational matrix.
    """
    #init mat_rot np.array of dim 3,3 default values set to 0
    mat_rot = np.zeros([3, 3])
    w, x, y, z = quaternion
    mat_rot[0,0] = w**2 + x**2 - y**2 - z**2
    mat_rot[1,0] = 2*(x*y + w*z)
    mat_rot[2,0] = 2*(x*z - w*y)
    mat_rot[0,1] = 2*(x*y - w*z)
    mat_rot[1,1] = w**2 - x**2 + y**2 - z**2
    mat_rot[2,1] = 2*(y*z + w*x)
    mat_rot[0,2] = 2*(x*z + w*y)
    mat_rot[1,2] = 2*(y*z - w*x)
    mat_rot[2,2] = w**2 - x**2 - y**2 + z**2
    return mat_rot


def apply_rotation(vec_to_rotate, rotational_vector, rad_angle):
    """Rotates a vector around another vector by a given angle.

    Parameters
    ----------
    vec_to_rotate : numpy 1D-array.
    vector to rotate around : numpy 1D-array.
    rad_angle : float.

    Returns
    -------
    numpy 1D-array
        The final rotated (normalized) vector.
    """
    # Generate a quaternion of the given angle (in radian)
    quat_rot = v2q(rotational_vector, rad_angle)
    # Generate the rotational matrix 
    rot_mat_quat = rotational_matrix(quat_rot)
    # Use the rotational matrix on the vector to rotate
    vec_rotated = np.dot(rot_mat_quat, vec_to_rotate)
    norm_vec = normalize(vec_rotated)
    return norm_vec


def read_pdb(pdb_filename):
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


def write_PDB(atom_num, atom_type, coor):
    """Prints atom coordinates in PDB format.

    Parameters
    ----------
    atom_num : int
    atom_type : string
    coor : numpy 1D-array

    Returns
    -------
    None
    """
    x, y, z = coor
    # pdb format (source: http://cupnet.net/pdb-format/)
    print("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}"
          "{:6.2f}{:6.2f}          {:>2s}{:2s}"
          .format("ATOM", atom_num, atom_type, "", "LIP", "A", 1,"",  x, y, z,
                  1.0, 0.0, "", ""))
   
    
def get_SP2_H(atom, helper1, helper2):
    """Reconstructs the 2 hydrogens of a SP2 carbon.

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
    tuple of numpy 1D-array
        Coordinates of the two hydrogens, e.g. ([x_H1, y_H1, z_H1],
        [x_H2, y_H2, z_H2]).
    """
    #atom - helper1 vector
    v2 = normalize(helper1 - atom)
    #atom - helper2 vector
    v3 = normalize(helper2 - atom)
    #####case(3) !CH2####
    #Perpendicular to the helpers - atom plane
    v4 = normalize(np.cross(v3, v2))
    #Rotational vector
    rot_vec = normalize(v2 - v3)
    #Vector to be rotated by theta/2, perpendicular to rot_vec and v4
    vec_to_rotate = normalize(np.cross(v4, rot_vec))
    norm_vec_H1 = apply_rotation(vec_to_rotate, rot_vec, -1.911/2)
    hcoor_H1 = 1 * norm_vec_H1 + atom
    norm_vec_H2 = apply_rotation(vec_to_rotate, rot_vec, 1.911/2)
    hcoor_H2 = 1 * norm_vec_H2 + atom
    return (hcoor_H1, hcoor_H2)


if __name__ == "__main__":
    # read coordinates in a pandas dataframe
    df_atoms = read_pdb("POPC.pdb")
    #print(df_atoms)
    # select only atoms N4
    #print(df_atoms[ (df_atoms["resname"] == "POP") &
    #                (df_atoms["atname"] == "N4") ])
    # select coor of atoms N4
    N4 = df_atoms[ (df_atoms["resname"] == "POP") &
                   (df_atoms["atname"] == "N4") ]
    # select coor only of N4
    N4_coor_only = N4[["x", "y", "z"]]
    #print(N4_coor_only)
    # convert N4_coor_only dataframe to an np 2D-array
    N4_2Darray = np.array(N4_coor_only.values.tolist())
    #print(N4_2Darray)
    # do the same on C5 and O11
    C5_2Darray = np.array(df_atoms[ (df_atoms["resname"] == "POP") &
                                    (df_atoms["atname"] == "C5") ] \
                          [["x", "y", "z"]].values.tolist())
    O11_2Darray = np.array(df_atoms[ (df_atoms["resname"] == "POP") &
                                     (df_atoms["atname"] == "C6") ] \
                           [["x", "y", "z"]].values.tolist())
    index = 1
    for i in range(len(N4_2Darray)):
        #print(i, C5_2Darray[i],
        #      N4_2Darray[i],
        #      O11_2Darray[i])
        coor_H1, coor_H2 = get_SP2_H(C5_2Darray[i],
                                     N4_2Darray[i],
                                     O11_2Darray[i])
        write_PDB(index, "C", C5_2Darray[i])
        index += 1 
        write_PDB(index, "H", coor_H1)
        index += 1 
        write_PDB(index, "H", coor_H2)
