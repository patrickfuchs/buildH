**Dev**
- Add first notebook (basic buildH analysis on a Berger traj)
- Complete documentation

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
