**Dev**

**1.6.0**

- Avoid output trajectory rewind when writing box dimensions
- Switch to MDAnalysis 2.0
- Add support of Python 3.9
- Improve docstrings

**1.5.0**

- Write box dimensions in the requested trajectory output
- Fix write duplicate 1st frame when a trajectory output is requested
- Avoid using universe.trajectory.time on a single pdb
- Limit Python version >= 3.6 <=3.8 (for MDAnalysis compatibility)
- Add support for: Berger DOPC/DPPC/POPS, GROMOS-CKP POPC/POPS, GROMOS-53A6L DPPC, CHARMM36UA
- Force Python 3.8 for doc building

**1.4.0**

- Add -v / --version option
- Reorganize doc
- Add Notebook04 (launch buildH as a module)
- Support Berger cholesterol
- Add Notebook05 (mixture POPC / cholesterol)
- Create buildH logo and add it to doc
- Add paper for JOSS
- Add community guidelines

**1.3.1**

- Fix setup.cfg to include json files in python package archive

**1.3.0**

- Complete documentation
- Accelerate functions within geometry.py with Numba
- Implement the use of buildH as a module
- Simplify calculation of CH on an sp3 carbon
- Use MyST parser for documentation (handles latex equations)
- Clarify some error messages
- Fix residue number exceeding 9999
- Add POPE def and json files
- Add Notebook01 (basic buildH analysis on a Berger traj)
- Add Notebook02 (+trajectory output)
- Add Notebook03 (analysis on a mixture POPC/POPE)
- Move CHARMM36 POPC validation to Zenodo

**1.2.0**

- Build docs
- Rename '-x/--xtc' flag to -t/--traj' one to be more generic
- Replace mandatory topology argument to '-c/--coord' flag
- Improve performance of control functions.
- Move misc functions to a module utils.py
- Improve Exception handling & add proper exits
- Improve PEP8 & PEP257 compliance
- Improve test coverage
- Fix bug when a trajectory was written when only a pdb was provided.
- Add sanity checks for the various input files
- Use json files instead of python module to read lipid topologies.
- Optimize package for better performance

**1.1.0**

- Create Python package structure
- Create conda environment
- Fix tests
- Separate entry point
- Update README for dev version installation
- Handle version with bump2version
