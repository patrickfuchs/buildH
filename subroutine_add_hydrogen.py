#!/usr/bin/env python


import numpy as np 

def normalize(vec):
    """Function that normalize a 3-D vector
	Take in a numpy array of dimension 3
	Return a normalized numpy array of dimension 3"""

    nvec = vec / np.sqrt(np.sum(vec**2))
    return nvec

def v2q(vec, theta):
    """Function that translate a 3-D vector and angle theta in a quaternion.
	Takes in a numpy array of dimension 3 and a theta angle (in radian)
	Returns a quaternion; numpy array of dimension 4"""

    w = np.cos(theta/2)
    (x, y, z) = np.sin(theta/2)*normalize(vec)
    q = np.array([w, x, y, z])
    return q


def rotational_matrix(quaternion):
    """Function that translate de quaternion into a rotational matrix
	Takes in a numpy array of dimension 4 i.e a quaternion
	Returns a rotational matrix; numpey array dimension [3,3]"""
    #init mat_rot np.array of dim 3,3 default values set to 0
    mat_rot = np.zeros([3,3])
    w = quaternion[0]
    x = quaternion[1]
    y = quaternion[2]
    z = quaternion[3]
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

#    1DOPC    N4    4   4.389   7.455   2.167
#    1DOPC    C5    5   4.438   7.592   2.185
#    1DOPC   H51    6   4.357   7.651   2.181
#    1DOPC   H52    7   4.497   7.611   2.106
#    1DOPC    C6    8   4.516   7.627   2.312

def apply_rotation(vec_to_rotate, rotational_vector, rad_angle):
    """Function that take a vector to rotate around an other vector 
    by a given angle. It returns the final vector normalised"""
    #generate a quaternion of the given angle (in radian)
    quat_rot = v2q(rotational_vector, rad_angle)
    #generate the rotational matrix 
    rot_mat_quat = rotational_matrix(quat_rot)
    #Use the rotational matrix on the vector to rotate
    vec_rotated = np.dot(rot_mat_quat, vec_to_rotate)
    norm_vec = normalize(vec_rotated)
    return norm_vec
    

def write_PDB(id, atom_type, coor):
    """Function that print information in a PDB format
    Helps check the result with VMD """
    x = coor[0]
    y = coor[1]
    z = coor[2]
    #From cupnet website ! 
    print "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format("ATOM",id, atom_type, "", "LIP", "A", 1,"",  x, y, z, 1.0, 0.0, "", "")
    
    
def get_SP2_H(atom, helper1, helper2):
    """Function that reconstruct the 2 hydrogens of a SP2 carbon
    It needs the carbon coordinate and 2 helpers, the previous atom and 
    the next one.
    It returns the coordinates of the hydrogens in a tuple of 3D list 
    ([x_H1, y_H1, z_H1], [x_H2, y_H2, z_H2])
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
    
    

#Traduction case(3) fortran subroutine add_hydrogen
##Tests##
#N4, C5, O11, C13
helper1 = np.array([[43.89, 74.55, 21.67], [44.38, 75.92, 21.85], [46.58  ,78.14   ,25.37], [45.10 , 79.68 ,  26.62]])
#C5, C6, C12, C32
atom = np.array([[44.38, 75.92, 21.85],[45.16 , 76.27 , 23.12], [45.77 , 79.32 , 25.29], [44.82 , 78.27,  27.15]])
#C6, O7, C13, O33
helper2 = np.array([[45.16, 76.27, 23.12], [45.59  , 77.63,  23.07], [45.10 , 79.68 ,  26.62], [43.94 , 78.22 , 28.27]])

#Tests of the news functions on simple cases 
index = 1
for i in range(len(atom)):
    (coor_H1, coor_H2) = get_SP2_H(atom[i], helper1[i], helper2[i])
    write_PDB(index, "C", atom[i])
    index += 1 
    write_PDB(index, "H", coor_H1)
    index += 1 
    write_PDB(index, "H", coor_H2)
    




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



