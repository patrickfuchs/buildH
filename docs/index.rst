==================================
Welcome to buildH's documentation!
==================================

**Version** |release|

    *Build hydrogen atoms from a united-atom MD of lipids and calculate the order parameters.*

**buildH** is a software that reads a united-atom (UA) trajectory of lipids, builds the hydrogens on it and calculates the order parameter on each C-H bond. **buildH** also allows to output the trajectory with the new reconstructed hydrogens.

**buildH** works in two modes:

1. A slow mode when an output trajectory is requested by the user. In this case, the whole trajectory including newly built hydrogens is written to this trajectory file. If other molecules are present (e.g. water, ions, etc.), they will just be copied to the output trajectory with the same coordinates. So far, only the xtc format is supported.
2. A fast mode without any output trajectory.

In both modes, the order parameters are calculated.

It is possible to select only a part of the lipid on which **buildH** will do his job (e.g. the polar head, the sn-1 aliphatic chain, etc) thanks to the `def file <def_format.md>`_.


**buildH** has been carefully validated as explained in :doc:`Validation of buildH <CHARMM36_POPC_validation/validation>`.
The algorithms used to reconstruct hydrogens are detailed in :doc:`Algorithms for building hydrogens <algorithms_Hbuilding>`
and the formulas for computing the order parameters in :doc:`Order parameters and statistics <order_parameter>`.

All basic geometrical operations in **buildH** are accelerated using `Numba <https://numba.pydata.org>`_. **buildH** is hosted on `Github <https://github.com/patrickfuchs/buildH>`_.

Motivation
==========

The initial motivation comes from the `NMRlipids <https://nmrlipids.blogspot.com/>`_ project.
As stated in this `post <https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html>`_,
there was a lack of suitable program for reconstructing hydrogens.
In the past, we used to use `g_protonate` in GROMACS 3 but this program has been removed in recent versions.

Our idea was to build our own implementation in Python using libraries such as ``MDAnalysis``, ``Numpy`` and ``Pandas``.

**buildH** is used actively in the recent projects of NMRlipids such as `NMRlipidsIVPEandPG <https://github.com/NMRLipids/NMRlipidsIVPEandPG>`_ or `Databank <https://github.com/NMRLipids/Databank>`_.
**buildH** can also be used by anyone willing to analyze the order parameters from a UA trajectory, or if one needs to have explicit hydrogens for some further analyzes.

Citations
=========

If you use buildH, please cite:

```
Santuz et al., (2021). buildH: Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters. Journal of Open Source Software, 6(65), 3521, https://doi.org/10.21105/joss.03521
```

License
=======

buildH is licensed under `BSD 3-Clause <https://github.com/patrickfuchs/buildH/blob/master/LICENSE.txt>`_.


Content
=======

User manual
-----------
.. toctree::
    :maxdepth: 2

    installation
    usage
    command_line_options
    def_format
    json_format

Tutorials
---------
.. toctree::
   :maxdepth: 2

   tutorials

Algorithms, OP calculations & validation
----------------------------------------

.. toctree::
   :maxdepth: 2

   algorithms


API documentation
-----------------
.. toctree::
   :maxdepth: 1

   api_reference

Changelog
---------

.. toctree::
   :maxdepth: 1

   changelog
