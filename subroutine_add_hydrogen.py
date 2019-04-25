#!/usr/bin/env python


import numpy as np 

'''Function that normalize a 3-D vector
	Take in a numpy array of dimension 3
	Return a normalized numpy array of dimension 3'''
def normalize(vec):
    nvec = vec / np.sqrt(np.sum(vec**2))
    return nvec

'''Function that translate a 3-D vector and angle theta in a quaternion.
	Takes in a numpy array of dimension 3 and a theta angle (in radian)
	Returns a quaternion; numpy array of dimension 4'''
def v2q(vec, theta):
    w = np.cos(theta/2)
    (x, y, z) = np.sin(theta/2)*normalize(vec)
    q = np.array([w, x, y, z])
    return q

'''Function that translate de quaternion into a rotational matrix
	Takes in a numpy array of dimension 4 i.e a quaternion
	Returns a rotational matrix; numpey array dimension [3,3]'''
def rotational_matrix(quaternion):
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


#Traduction case(3) fortran subroutine add_hydrogen
##Tests##
N4 = np.array([43.89, 74.55, 21.67])
#v1
C5 = np.array([44.38, 75.92, 21.85])
C6 = np.array([45.16, 76.27, 23.12])

v1 = C5
#ref vectors
#N4-C5 vector
v2 = normalize(N4 - v1)
#C6-C5 vector
v3 = normalize(C6 - v1)


#####case(3) !CH2####
#Perpendicular to the N4-C5-C6 plane
v4 = normalize(np.cross(v3, v2))
#Rotational vector
u = normalize(v2 - v3)
#Vector to be rotated by theta/2, perpendicular to u and v4
#Comm' Amelie : why theta/2 ?? 
v4 = normalize(np.cross(v4, u))

## Part : Rotate the new v4 around u by theta
#quaternion of 109.5 degrees aka 1.911rad around v_N4_C5 normalized vector
quat_N4_C5 = v2q(u, -1.911/2)

#generate the rotational matrix 
rot_mat_quat_N4_C5 = rotational_matrix(quat_N4_C5)

#Use the rotational matrix on the C5-C6 vector
vec_H51 = np.dot(rot_mat_quat_N4_C5, v4)

norm_vec_H51 = normalize(vec_H51)
hcoor = 1 * norm_vec_H51 + C5
##End part Rotate the new v4 around u by theta

#+ 109.5°
quat_N4_C5 = v2q(u, 1.911/2)

#generate the rotational matrix 
rot_mat_quat_N4_C5 = rotational_matrix(quat_N4_C5)

#Use the rotational matrix on the C5-C6 vector
vec_H52 = np.dot(rot_mat_quat_N4_C5, v4)

norm_vec_H52 = normalize(vec_H52)
hcoor = 1 * norm_vec_H52 + C5



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



