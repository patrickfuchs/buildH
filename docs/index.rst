Welcome to buildH's documentation!
============================================

**Version** |release|

    Build hydrogens from a united-atom MD of lipids and calculate the order parameter.

Installation
============

Install buildH with `pip`:

.. code-block:: bash

    $ python3 -m pip install buildh

Installation within a conda environment:

.. code-block:: bash

   $ conda create -n env_buildH python pip
   $ conda activate env_buildH
   $ pip install buildh
   
User manual
===========
.. toctree::
    :maxdepth: 1

    buildh
    algorithms_Hbuilding
    command_line_options
    def_format
    json_format
    
Tutorial
========
.. toctree::
    :maxdepth: 1

    notebooks/README.md

Reference manual
================
.. toctree::
    :maxdepth: 2

    api/core
    api/geometry
    api/hydrogens
    api/init_dics
    api/lipids
    api/utils
    api/writers


Index
=====

* :ref:`genindex`


Changelog
=========

.. include:: ../CHANGELOG.md
