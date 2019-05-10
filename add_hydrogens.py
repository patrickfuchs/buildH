#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

import dic_lipids

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
                #dictmp = {}
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
          .format("ATOM", atom_num, atom_type, "", "POP", "", 1,"",  x, y, z,
                  1.0, 0.0, "", ""))
   
def pandasdf2pdb(df):
    """Returns a string from a pandas dataframe.
    TODO
    """
    s = ""
    chain = ""
    for i, row_atom in df.iterrows():
        atnum, atname, resname, resnum, x, y, z = row_atom
        atnum = int(atnum)
        resnum = int(resnum)
        s += ("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}"
             "{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
             .format("ATOM", atnum, atname, "", resname, chain, resnum, "",  x, y, z,
                     1.0, 0.0, "", ""))
    return s
 
    
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
        Coordinates of the two hydrogens: 
        ([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2]).
    """
    # Atom -> helper1 vector.
    v2 = normalize(helper1 - atom)
    # Atom -> helper2 vector.
    v3 = normalize(helper2 - atom)
    # Vector perpendicular to the helpers - atom plane.
    v4 = normalize(np.cross(v3, v2))
    # Rotational axis vector.
    rot_vec = normalize(v2 - v3)
    # Vector to be rotated by theta/2, perpendicular to rot_vec and v4
    vec_to_rotate = normalize(np.cross(v4, rot_vec))
    # Reconstruct the two hydrogens.
    norm_vec_H1 = apply_rotation(vec_to_rotate, rot_vec, -1.911/2)
    hcoor_H1 = 1 * norm_vec_H1 + atom
    norm_vec_H2 = apply_rotation(vec_to_rotate, rot_vec, 1.911/2)
    hcoor_H2 = 1 * norm_vec_H2 + atom
    return (hcoor_H1, hcoor_H2)


def get_name_H(name_carbon):
    name_H1 = name_carbon.replace("C", "H") + "1"
    name_H2 = name_carbon.replace("C", "H") + "2"
    return name_H1, name_H2


if __name__ == "__main__":
    # read coordinates in a pandas dataframe
    df_atoms = read_pdb("1POPC.pdb")
    # create an empty data frame to store the new mlc with added hydrogens
    new_df_atoms = pd.DataFrame(columns=["atnum", "atname", "resname",
                                 "resnum", "x", "y", "z"])
    new_atom_num = 1
    # loop over all existing atoms (the iter var row_atom is a pandas Series)
    for i, row_atom in df_atoms.iterrows():
        # add that atom to the new dataframe
        new_atom = row_atom
        new_atom["atnum"] = new_atom_num
        new_atom.name = new_atom_num
        new_df_atoms = new_df_atoms.append(new_atom)
        new_atom_num += 1
        # check whether it needs hydrogen(s) to be reconstructer
        atom_name = row_atom["atname"]
        if atom_name in dic_lipids.POPC:
            # get atom info
            atom_coor = np.array(row_atom[["x", "y", "z"]].values, dtype=float) # force to float
            res_name, res_num = row_atom[["resname", "resnum"]]
            # get helper atom names
            helper1_name, helper2_name = dic_lipids.POPC[atom_name]
            # get helper coords (needs [0] because it comes from a dataframe)
            helper1_coor = df_atoms [ (df_atoms["resnum"] == res_num) &
                                      (df_atoms["atname"] == helper1_name) ] \
                                      [["x", "y", "z"]].values[0]
            helper2_coor = df_atoms [ (df_atoms["resnum"] == res_num) &
                                      (df_atoms["atname"] == helper2_name) ] \
                                      [["x", "y", "z"]].values[0]
            # construct H1 & H2
            H1_coor, H2_coor = get_SP2_H(atom_coor, helper1_coor, helper2_coor)
            H1_name, H2_name = get_name_H(atom_name)
            # now create H1 & H2 as a pandas Series and append them into the new dataframe
            H1_atom = pd.Series([new_atom_num, H1_name, res_name, res_num]
                                + list(H1_coor),
                                name = new_atom_num,
                                index=["atnum", "atname", "resname", "resnum",
                                       "x", "y", "z"])
            new_df_atoms = new_df_atoms.append(H1_atom)
            new_atom_num += 1   
            H2_atom = pd.Series([new_atom_num, H2_name, res_name, res_num]
                                + list(H2_coor),
                                name = new_atom_num,
                                index=["atnum", "atname", "resname", "resnum",
                                       "x", "y", "z"])
            new_df_atoms = new_df_atoms.append(H2_atom)
            new_atom_num += 1
        #if i % 100 == 0:
        #    print("Atom {:5d} done!".format(i))
    #print(new_df_atoms)
    print(pandasdf2pdb(new_df_atoms))

exit()
#########
# Note by P@t: 28/04/2019
# Below is Amelie's stuff I didn't touch to
#########

######case(1) CH 
C13 = np.array([ 45.10 , 79.68,  26.62])
#C12 
helper1 = np.array([45.77, 79.32, 25.29])
#C32
helper2 = np.array([44.82,  78.27, 27.15])
#O14
helper3 = np.array([45.81,   80.59,  27.47])
helpers = np.array([[45.77, 79.32, 25.29],[44.82,  78.27, 27.15], [45.81,   80.59,  27.47]])
v2 = 0.0
for i in range(len(helpers)):
    v2 = v2 + normalize(helpers[i] - C13)

v2 = v2 / (len(helpers)) + C13 
coor_H = 1 * normalize(C13 - v2)+C13


write_PDB(1, "C", C13)
write_PDB(2, "C", helper1)
write_PDB(3, "C", helper2)
write_PDB(4, "O", helper3)
write_PDB(5, "H", coor_H)

######case(2) CH double bond
#Fonctions qui seront surement remplacées par une fonction issues de MDAnalysis
def vect_AB(A,B):
    """Returns a vector from point A to point B
    """
    return [B[0]-A[0],B[1]-A[1],B[2]-A[2]]
    
def scalar(A,B):
    """Returns the scalar (or inner) product between vectors A & B
    """
    return (A[0]*B[0]) + (A[1]*B[1]) + (A[2]*B[2])

def magnitude(A):
    """Returns the magnitude of vector A
    """
    return math.sqrt(A[0]**2+A[1]**2+A[2]**2)

def rad2deg(ang):
    """Convert an angle in radians to degrees
    """
    return ang*(180/math.pi)


def angle(A,B,C):
    """Returns the angle IN RAD between 3 points A, B & C
    """
    # compute vector BA
    vectBA = vect_AB(B,A) # [(A[0]-B[0]),(A[1]-B[1]),(A[2]-B[2])]
    # compute vector BC
    vectBC = vect_AB(B,C) # [(C[0]-B[0]),(C[1]-B[1]),(C[2]-B[2])]
    # compute the cosine of angle ABC ( BA * BC = ba.bc.cos(theta) )
    costheta = scalar(vectBA,vectBC)/(magnitude(vectBA)*magnitude(vectBC))
    # compute the angle ABC
    theta = math.acos(costheta)
    return rad2deg(theta)



#C24 
C24 = np.array([05.82, 31.07,  33.03])
#C23 
helper1 = np.array([06.66,  31.71,  31.93])
v2 = helper1 - C24
#C25
helper2 = np.array([06.29,  30.67,  34.27])
v3 = helper2 - C24

#thetal is the angle 2pi - C-C-C devided by 2
#to ensure equal (C-C-H) angles from both directions
angle_Cs = angle(helper1, C24, helper2)
#Dans le code fortran angle_Cs est divise par 180 ? passage en rad ? 
theta = math.pi * (2 - angle_Cs/180.) /2 
u = normalize(np.cross(v2, v3))
norm_vec_H = apply_rotation(v3, u, theta)
coor_H = 1 * norm_vec_H + C24

write_PDB(1, "C", helper1)
write_PDB(2, "C", C24)
write_PDB(3, "H", coor_H)
write_PDB(4, "C", helper2)

#C25 
C25 = np.array([06.29,   30.67,  34.273])
#C24 
helper1 = np.array([05.82,   31.07,   33.03])
v2 = helper1 - C25
#C26
helper2 = np.array([07.80,   30.74,   34.55])
v3 = helper2 - C25

#thetal is the angle 2pi - C-C-C devided by 2
#to ensure equal (C-C-H) angles from both directions
angle_Cs = angle(helper1, C25, helper2)
#Dans le code fortran angle_Cs est divise par 180 ? passage en rad ? 
theta = math.pi * (2 - angle_Cs/180.) /2 
u = normalize(np.cross(v2, v3))
norm_vec_H = apply_rotation(v3, u, theta)
coor_H = 1 * norm_vec_H + C25

write_PDB(1, "C", helper1)
write_PDB(2, "C", C25)
write_PDB(3, "H", coor_H)
write_PDB(4, "C", helper2)



#####case(4) !CH3e
theta = 1.911
C50 = np.array([25.26 ,  08.20,   24.58])
C49 = np.array([25.99 ,  07.97 ,  25.91])
v2 = C49 - C50
C48 = np.array([24.97,   07.55 ,  26.97])
v3 = C48- C50

u = normalize(np.cross(v3, v2))
norm_vec_H = apply_rotation(v2, u, theta)
coor_H = 1 * norm_vec_H + C50

write_PDB(1, "C", C48)
write_PDB(2, "C", C49)
write_PDB(3, "C", C50)
write_PDB(4, "H", coor_H)

####case(5) !CH3r
theta = 120 * math.pi / 180.
C50 = np.array([25.26 ,  08.20,   24.58])
C49 = np.array([25.99 ,  07.97 ,  25.91])
H3 = np.array([ 25.91960191,   8.47515027,  23.88055904])
u = normalize(C49 - C50)
v4 = normalize(H3 - C50) 
norm_vec_H = apply_rotation(v4, u, theta)
coor_H = 1 * norm_vec_H + C50 

write_PDB(5, "H", coor_H)

####case(5) !CH3s
theta = -120 * math.pi / 180.
C50 = np.array([25.26 ,  08.20,   24.58])
C49 = np.array([25.99 ,  07.97 ,  25.91])
H3 = np.array([ 25.91960191,   8.47515027,  23.88055904])
u = normalize(C49 - C50)
v4 = normalize(H3 - C50) 
norm_vec_H = apply_rotation(v4, u, theta)
coor_H = 1 * norm_vec_H + C50 

write_PDB(6, "H", coor_H)

#####case(4) !CH3e
theta = 1.911
CA2 = np.array([07.52,   23.79,   38.89])
CA1 = np.array([06.85,  24.93,   38.12])
v2 = CA1 - CA2
C31 = np.array([ 07.99,   25.93,   37.91])
v3 = C31- CA2

u = normalize(np.cross(v3, v2))
norm_vec_H = apply_rotation(v2, u, theta)
coor_H = 1 * norm_vec_H + CA2

write_PDB(1, "C", C31)
write_PDB(2, "C", CA1)
write_PDB(3, "C", CA2)
write_PDB(4, "H", coor_H)

####case(5) !CH3r
theta = 120 * math.pi / 180.
CA2 = np.array([07.52,   23.79,   38.89])
CA1 = np.array([06.85,  24.93,   38.12])
H3 = np.array([ -4.758,  -2.206,  26.417])
u = normalize(CA1 - CA2)
v4 = normalize(H3 - CA2) 
norm_vec_H = apply_rotation(v4, u, theta)
coor_H = 1 * norm_vec_H + CA2 

write_PDB(5, "H", coor_H)

####case(5) !CH3s
theta = -120 * math.pi / 180.
CA2 = np.array([07.52,   23.79,   38.89])
CA1 = np.array([06.85,  24.93,   38.12])
H3 = np.array([ -4.758,  -2.206,  26.417])
u = normalize(CA1 - CA2)
v4 = normalize(H3 - CA2) 
norm_vec_H = apply_rotation(v4, u, theta)
coor_H = 1 * norm_vec_H + CA2 

write_PDB(6, "H", coor_H)

