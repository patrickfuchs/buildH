**Dev**

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
