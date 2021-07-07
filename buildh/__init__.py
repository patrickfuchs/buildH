"""buildH builds hydrogen atoms from united-atom molecular dynamics of lipids.

It also calculates the order parameters.

It works in two modes :
  1) A slow mode when an output trajectory (e.g. in xtc format) is requested by
     the user. In this case, the whole trajectory including newly built
     hydrogens are written to this trajectory.
  2) A fast mode without any output trajectory.
For both modes, the order parameter is written to an output file in a format
similar to the code of @jmelcr:
https://github.com/NMRLipids/MATCH/blob/master/scripts/calcOrderParameters.py

This code has been checked against the one from @jmelcr. You might find minor
differences due to rounding errors (in xtc, only 3 digits are written).

The way of building H is largely inspired from a code of Jon Kapla originally
written in fortran :
https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242.

Note: that all coordinates in this script are handled using numpy 1D-arrays
of 3 elements, e.g. atom_coor = np.array((x, y, z)).
Note2: sometimes numpy is slow on small arrays, thus we wrote a few "in-house"
functions for vectorial operations (e.g. cross product).
"""

__authors__ = ("Hubert Santuz", "Amélie Bâcle",
               "Pierre Poulain", "Patrick Fuchs")
__email__ = "patrick.fuchs@u-paris.fr"
__version__ = "1.5.0"
__license__ = "BSD 3-Clause License"

# Silent this unharmful warning from MDAnalysis
import warnings
warnings.filterwarnings("ignore", message="Element information is absent or missing for a few ")

from .UI import BuildHError, launch
