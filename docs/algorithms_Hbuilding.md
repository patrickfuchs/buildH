# Algorithm for building hydrogens

The way of building H is largely inspired from [a code of Jon Kapla](https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242) originally written in fortran.

## Building CH2

**TODO**: to be updated.

Below is an example of a reconstruction of 2 hydrogens (*H51* and *H52*) attached to a carbon *C5* with the help of 2 others atom *C6* and *N4*.

![Vectors](vectors.png)

First, we compute the cross product between the vector C5-C6 and C5-N4 (red).

We determine a rotational axis determined by the vector N4-C6. (green)

We compute the cross product between the red one and green one. (orange)

This orange vector is the rotational vector to construct the hydrogens.

For this case, 2 hydrogens are constructed (yellow) : we apply a rotation of 109.47 deg for one and -109.47 deg for the other one.

## TODO

How the other hydrogens are reconstructed (CH3, CH double bond, CH).
