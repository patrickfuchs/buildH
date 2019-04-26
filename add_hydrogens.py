#!/usr/bin/env python

import numpy as np 

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
    #Traduction case(3) fortran subroutine add_hydrogen
    ##Tests##
    #N4, C5, O11, C13
    helper1 = np.array([[43.89, 74.55, 21.67],
                        [44.38, 75.92, 21.85],
                        [46.58, 78.14, 25.37],
                        [45.10, 79.68, 26.62]])
    #C5, C6, C12, C32
    atom = np.array([[44.38, 75.92, 21.85],
                     [45.16 , 76.27 , 23.12],
                     [45.77 , 79.32 , 25.29],
                     [44.82 , 78.27,  27.15]])
    #C6, O7, C13, O33
    helper2 = np.array([[45.16, 76.27, 23.12],
                        [45.59  , 77.63,  23.07],
                        [45.10 , 79.68 ,  26.62],
                        [43.94 , 78.22 , 28.27]])
    #Tests of the news functions on simple cases 
    index = 1
    for i in range(len(atom)):
        (coor_H1, coor_H2) = get_SP2_H(atom[i], helper1[i], helper2[i])
        write_PDB(index, "C", atom[i])
        index += 1 
        write_PDB(index, "H", coor_H1)
        index += 1 
        write_PDB(index, "H", coor_H2)
    


exit()

#####case(4) !CH3e
v1 = C5
#ref vectors
#N4-C5 vector
v2 = normalize(N4 - v1)
#C6-C5 vector
v3 = normalize(C6 - v1)

u = normalize(np.cross(v3, v2))

quat_N4_C5 = v2q(u, -1.911/2)

#generate the rotational matrix 
rot_mat_quat_N4_C5 = rotational_matrix(quat_N4_C5)

#Use the rotational matrix on the C5-C6 vector
vec_H51 = np.dot(rot_mat_quat_N4_C5, v2)

norm_vec_H51 = normalize(vec_H51)
hcoor = 0.1 * norm_vec_H51 + C5



