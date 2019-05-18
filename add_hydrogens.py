#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pandas as pd
import MDAnalysis

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
    return vec / magnitude(vec)


def magnitude(vec):
    """Returns the magnitude of a vector.

    Parameters
    ----------
    vec : numpy 1D-array.

    Returns
    -------
    float
        The magniture of the vector.
    """
    return np.sqrt(np.sum(vec**2))


def calc_angle(atom1, atom2, atom3):
    """Calculate the valence angle between atom1, atom2 and atom3.

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


def apply_rotation(vec_to_rotate, rotational_axis, rad_angle):
    """Rotates a vector around another vector by a given angle.

    Parameters
    ----------
    vec_to_rotate : numpy 1D-array.
    rotational_axis : numpy 1D-array.
    rad_angle : float.

    Returns
    -------
    numpy 1D-array
        The final rotated (normalized) vector.
    """
    # Generate a quaternion of the given angle (in radian).
    quaternion = v2q(rotational_axis, rad_angle)
    # Generate the rotation matrix.
    rotation_matrix = rotational_matrix(quaternion)
    # Apply the rotational matrix on the vector to rotate.
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
        s += ("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}"
              "{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
              .format("ATOM", atnum, atname, "", resname, chain, resnum, "",
                      x, y, z, 1.0, 0.0, "", ""))
    return s
 
    
def get_CH2(atom, helper1, helper2):
    """Reconstructs the 2 hydrogens of a sp3 carbon.

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
    # atom -> helper1 vector.
    v2 = normalize(helper1 - atom)
    # atom -> helper2 vector.
    v3 = normalize(helper2 - atom)
    # Vector perpendicular to the helpers/atom plane.
    v4 = normalize(np.cross(v3, v2))
    # Rotational axis.
    rotation_axis = normalize(v2 - v3)
    # Vector to be rotated by theta/2, perpendicular to rot_vec and v4.
    vec_to_rotate = normalize(np.cross(v4, rotation_axis))
    # Reconstruct the two hydrogens.
    norm_vec_H1 = apply_rotation(vec_to_rotate, rotation_axis, -1.911/2)
    hcoor_H1 = 1 * norm_vec_H1 + atom
    norm_vec_H2 = apply_rotation(vec_to_rotate, rotation_axis, 1.911/2)
    hcoor_H2 = 1 * norm_vec_H2 + atom
    return (hcoor_H1, hcoor_H2)


def get_CH(atom, helper1, helper2, helper3):
    helpers = np.array((helper1, helper2, helper3))
    v2 = np.zeros(3)
    for i in range(len(helpers)):
        v2 = v2 + normalize(helpers[i] - atom)
    v2 = v2 / (len(helpers)) + atom
    coor_H = 1 * normalize(atom - v2) + atom
    return coor_H


def get_CH_double_bond(atom, helper1, helper2):
    # calc angle theta helper1-atom-helper2 (in rad).
    theta = calc_angle(helper1, atom, helper2)
    # atom -> helper1 vector.
    v2 = helper1 - atom
    # atom -> helper2 vector.
    v3 = helper2 - atom
    # The rotation axis is orthogonal to the atom/helpers plane.
    rotation_axis = normalize(np.cross(v2, v3))
    # Reconstruct H by rotating v3 by theta.
    norm_vec_H = apply_rotation(v3, rotation_axis, theta)
    coor_H = 1 * norm_vec_H + atom
    return coor_H


def get_name_H(name_carbon, nb_of_H):
    if nb_of_H == 1:
        name_H1 = name_carbon.replace("C", "H") + "1"
        return name_H1
    elif nb_of_H == 2:
        name_H1 = name_carbon.replace("C", "H") + "1"
        name_H2 = name_carbon.replace("C", "H") + "2"
        return name_H1, name_H2


#@profile
def reconstruct_hydrogens3(list_df_residues):
    # The list newrows will be used to store the new molecule *with* H.
    newrows = []
    # Counter for numbering the new mlcs with H.
    new_atom_num = 1
    # Loop over all residues.
    for df_residue in list_df_residues:
        # Loop over all existing atoms within that residue.
        # (iter var row_atom is a pandas Series)
        for i, row_atom in df_residue.iterrows():
            # Renum atom (since we'll have additional H, atom num will change).
            row_atom["atnum"] = new_atom_num
            # Append row_atom to the new list (columns 0 to 6 included,
            # we don't need the following ones).
            # 0      1       2        3       4  5  6
            # atnum, atname, resname, resnum, x, y, z
            newrows.append(list(row_atom[0:7]))
            new_atom_num += 1
            # Check whether the atom needs hydrogen(s) to be reconstructed onto.
            if row_atom["typeofH2build"]:
                # Get atom info (force float for the coordinates!).
                atom_coor = np.array(row_atom[["x", "y", "z"]].values,
                                     dtype=float)
                #####
                ##### TRY TO USE .to_numpy() method
                ##### 
                res_name, res_num = row_atom[["resname", "resnum"]]
                # Get name of helper atoms.
                helper1_name, helper2_name = (row_atom["helper1_name"],
                                              row_atom["helper2_name"])
                # Get helper coords.
                #print(row_atom) ; exit()
                #helper1_coor = np.array(df_residue.loc[helper1_name]
                #                        [["x", "y", "z"]].values, dtype=float)
                #helper2_coor = np.array(df_residue.loc[helper2_name]
                #                        [["x", "y", "z"]].values, dtype=float)
                helper1_coor = np.array(df_residue.loc[helper1_name,
                                                       ["x", "y", "z"]]
                                        .values, dtype=float)
                helper2_coor = np.array(df_residue.loc[helper2_name,
                                                       ["x", "y", "z"]]
                                        .values, dtype=float)
                # Build H(s).
                H1_coor, H2_coor = get_CH2(atom_coor, helper1_coor,
                                             helper2_coor)
                # Add new H(s) to the newrows list.
                H1_name, H2_name = get_name_H(row_atom["atname"], 2)
                newrows.append([new_atom_num, H1_name, res_name, res_num]
                               + list(H1_coor))
                new_atom_num += 1
                newrows.append([new_atom_num, H2_name, res_name, res_num]
                               + list(H2_coor))
                new_atom_num += 1
    # Create a dataframe to store the mlc with added hydrogens.
    new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname", "resname",
                                                  "resnum", "x", "y", "z"])
    return new_df_atoms


def reconstruct_hydrogens_wMDanalysis(pdb_filename, return_coors=False):
    # load PDB
    u = MDAnalysis.Universe(pdb_filename)
    if return_coors:
        # The list newrows will be used to store the new molecule *with* H.
        newrows = []
        # Counter for numbering the new mlcs with H.
        new_atom_num = 1
    # Loop over all atoms.
    for atom in u.atoms:
        if return_coors:
            resnum = atom.resnum
            resname = atom.resname[:-1] # beware, resname must be 3 letters long in my routine
            name = atom.name
            # Append atom to the new list.
            # 0      1       2        3       4  5  6
            # atnum, atname, resname, resnum, x, y, z
            newrows.append( [new_atom_num, name, resname, resnum] + list(atom.position) )
            new_atom_num += 1
        if atom.name in dic_lipids.POPC:
            typeofH2build = dic_lipids.POPC[atom.name][0]
            if typeofH2build == "CH2":
                _, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
                # atom is a Atom object : atom.residue.atoms is a list of its atom we can select.
                # [0] is because select_atoms return a AtomGroup which contains only 1 atom.
                helper1_coor = atom.residue.atoms.select_atoms("name {0}".format(helper1_name))[0].position
                helper2_coor = atom.residue.atoms.select_atoms("name {0}".format(helper2_name))[0].position
                # Build H(s).
                H1_coor, H2_coor = get_CH2(atom.position, helper1_coor, helper2_coor)
                ####
                #### We could calculate here the order parameter on the fly :-D !
                #### call_routine_op()
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name, H2_name = get_name_H(atom.name, 2)
                    newrows.append( [new_atom_num, H1_name, resname, resnum] + list(H1_coor) )
                    new_atom_num += 1
                    newrows.append( [new_atom_num, H2_name, resname, resnum] + list(H2_coor) )
                    new_atom_num += 1
            elif typeofH2build == "CH":
                _, helper1_name, helper2_name, helper3_name = dic_lipids.POPC[atom.name]
                helper1_coor = atom.residue.atoms.select_atoms("name {0}".format(helper1_name))[0].position
                helper2_coor = atom.residue.atoms.select_atoms("name {0}".format(helper2_name))[0].position
                helper3_coor = atom.residue.atoms.select_atoms("name {0}".format(helper3_name))[0].position
                H1_coor = get_CH(atom.position, helper1_coor, helper2_coor, helper3_coor)
                ####
                #### We could calculate here the order parameter on the fly :-D !
                #### call_routine_op()
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name = get_name_H(atom.name, 1)
                    newrows.append( [new_atom_num, H1_name, resname, resnum] + list(H1_coor) )
                    new_atom_num += 1
            elif typeofH2build == "CHdoublebond":
                _, helper1_name, helper2_name = dic_lipids.POPC[atom.name]
                helper1_coor = atom.residue.atoms.select_atoms("name {0}".format(helper1_name))[0].position
                helper2_coor = atom.residue.atoms.select_atoms("name {0}".format(helper2_name))[0].position
                H1_coor = get_CH_double_bond(atom.position, helper1_coor, helper2_coor)
                ####
                #### We could calculate here the order parameter on the fly :-D !
                #### call_routine_op()
                ####
                if return_coors:
                    # Add them to newrows.
                    H1_name = get_name_H(atom.name, 1)
                    newrows.append( [new_atom_num, H1_name, resname, resnum] + list(H1_coor) )
                    new_atom_num += 1
    if return_coors:
        # Create a dataframe to store the mlc with added hydrogens.
        new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname", "resname",
                                                      "resnum", "x", "y", "z"])
        return new_df_atoms


if __name__ == "__main__":
    use_pandas = False
    pdb_filename = "POPC_only.pdb"#"1POPC.pdb"
    if use_pandas:
        # read coordinates in a pandas dataframe
        list_df_residues = pdb2list_pandasdf_residues(pdb_filename)
        new_df_atoms = reconstruct_hydrogens3(list_df_residues)
        print(pandasdf2pdb(new_df_atoms))
    else:
        new_df_atoms = reconstruct_hydrogens_wMDanalysis(pdb_filename, return_coors=True)
        print(pandasdf2pdb(new_df_atoms))

exit()
#########
# Note by P@t: 28/04/2019
# Below is Amelie's stuff I didn't touch to
#########

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


def oldmagnitude(A):
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
    costheta = scalar(vectBA,vectBC)/(oldmagnitude(vectBA)*oldmagnitude(vectBC))
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

