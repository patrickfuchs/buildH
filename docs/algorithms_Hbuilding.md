# Algorithm for building hydrogens

**buildH** builds hydrogens using general geometric rules which are explained in this document. All the Python functions implementing these reconstructions are written in `hydrogens.py`. These functions are largely inspired from [a code of Jon Kapla](https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242) originally written in fortran. All mathematical functions (vector operations, rotations, etc.) are written in `geometry.py` and accelarated using [Numba](https://numba.pydata.org/). In this page, we use the following conventions:

- all represented vectors are unit vectors;
- $\theta$ is the [tetrahedral bond angle](https://en.wikipedia.org/wiki/Tetrahedron) which equals ~109.5°;
- $l_{CH}$ is the [carbon-hydrogen bond length](https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond) which equals ~ 1.09 Å.

All images in this page were generated using [VMD](http://www.ks.uiuc.edu/Research/vmd/).

## Building CH3

Building a methyl on a primary carbon requires two helpers: i) helper1 is connected to that carbon, ii) helper2 is connected to helper1. We start with the reconstructruction of the first hydrogen as explained in the figure below.

![CH3_building](img/how_CH3_building1.png)

`vect1` (red) is first computed as the vector product between vectors "carbon -> helper2" and "carbon -> helper1", it will be our rotation axis in the next step. `vect2` (blue) is then computed by rotating vector "carbon -> helper1" about `vect1` by $\theta$. The first H will be obtained by translating a point located at the "carbon" along `vect2` of $l_{CH}$ Å. Note that the newly built H is in a *trans* configuration with respect to helper2.

We then go on with the reconstruction of the two other Hs as explained in the figure below.

![CH_building](img/how_CH3_building2.png)

`vect3` (green) is obtained by rotating `vect2` of $+\frac{2\pi}{3}$ about vector "carbon -> helper1". 
`vect4` (magenta) is obtained by rotating `vect2` of $-\frac{2\pi}{3}$ about vector "carbon -> helper1". 

The second and third H will be obtained by translating of $l_{CH}$ Å a point located at the "carbon" along `vect3` and `vect4` respectively.

## Building CH2

The building of 2 hydrogens on a secondary carbon involves a few geometrical procedures that are explained in the figure below.

![CH2_building](img/how_CH2_building.png)

We start with the 3 atoms, the central carbon on which we want to reconstruct hydrogens (`C26`), helper1 (`C25`) and helper2 (`C27`) which are connected to the central carbon. The two helpers will help us build the new hydrogens following standard [tetrahedral geometry](https://en.wikipedia.org/wiki/Tetrahedral_molecular_geometry). 

On the left panel, we first show how to construct 3 vectors:

- `vect1` (red) is normal to the plane of the 3 atoms. It is calculated as the cross product between vectors "central carbon -> helper2" and "central carbon -> helper1".
- `vect2` (blue) will be our **rotation axis** used later. It is calculated as vector "central carbon -> helper1" minus vector "central carbon -> helper2".
- `vect3` (green) is a vector that will be rotated in the next step. It is the cross product between `vect1` and `vect2`.

On the right panel, we go on to construct 2 other vectors:

- `vect4` (magenta) is obtained by rotating `vect3` of $\frac{\theta}{2}$ about `vect2`. The first H will be obtained by translating a point located at the central carbon along `vect4` of $l_{CH}$ Å.
- `vect5` (orange) is obtained by rotating `vect3` of $-\frac{\theta}{2}$ about `vect2`. The second H will be obtained by translating a point located at the central carbon along vect5 of $l_{CH}$ Å.

## Building CH

The building of 1 hydrogen on a tertiary carbon is quite simple and explained in the figure below.

![CH_building](img/how_CH_building.png)

We first compute the red vector `vect1` as the sum of the 3 vectors "central carbon -> helper1" + "central carbon -> helper2" + "central carbon -> helper3". We see that this `vect1` defines a [median](https://en.wikipedia.org/wiki/Median_(geometry)#Tetrahedron) of the tetrahedron. The blue vector `vect2` is merely the opposite of `vect1` and gives the direction of the C-H bond. The new H will be obtained by translating a point located at the central carbon along `vect2` of $l_{CH}$ Å.

## Building CH on a double bond

For the H to reconstruct on a carbon involved in a double bond, **buildH** uses the following strategy as explained in the figure below.

![CH_building](img/how_CHdoublebond_building.png)

First we compute the angle $\gamma$ between atoms "helper1-central carbon-helper2". `vect1` is next calculated as the cross product between vectors "central carbon -> helper1" and "central carbon -> helper2". `vect1` will be used in the next step as a rotation axis. `vect2` is finally obtained by rotating vector "helper2 -> atom" of $\pi-\frac{\gamma}{2}$ rad about `vect1`. Why this angle value? In fact, **buildH** uses here the bisection strategy. `vect2` is along the same axis as the vector bisecting the angle "helper1-central carbon-helper2" but on the opposite direction. To obtain this bisecting vector, we would need to rotate vector "helper2 -> atom" of $-\frac{gamma}{2}. Since `vect2` is on the opposite direction, we simply add $\pi$ to that value obtaining thus $\pi-\frac{\gamma}{2}$ rad.

The double bond case is somewhat more complicated as largely discussed in an article published in JCTC from [Piggot at al.](https://doi.org/10.1021/acs.jctc.7b00643). One of the problem is that the ideal CCC angle (helper1-carbon-helper2) of 120° can vary according to the functional groups and some details in the force field. In **buildH** we decided to use the bisection strategy (as described above) so that the H reconstruction will adapt itself to any value of ${\gamma}$. Comparing the resulst with H reconstructed like this vs the real ones from an all-atom CHARMM36 snapshot of 256 POPC yielded a acceptable difference in $S_{CH}$ of 0.008 and 0.016 (see [here](https://zenodo.org/record/4715962)). This difference was even reduced upon averaging over a trajectory. 
