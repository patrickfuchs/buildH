# Validation of buildH on a CHARMM36 POPC trajectory

This directory contains a [report](report_buildH.pdf) describing a validation of **buildH** made in August 2019.

**buildH** reconstructs hydrogens from a united-atom trajectory and calculates the order parameter on each reconstructed C-H bond. To validate whether buildH works, we took an all-atom trajectory (generated with the CHARMM36 force field), removed the hydrogens and reconstructed them with buildH. Then we compared the H reconstruction and the order parameter values calculated with buildH to the real ones from the all-atom trajectory.

All the files used for making this validation have been deposited on [Zenodo](https://zenodo.org/record/4715962) with the following DOI: 10.5281/zenodo.4715962.

The file [report_buildH.pdf](report_buildH.pdf) (taken from the Zenodo archive) describes the validation.

The file [CHARMM_vs_buildH.gif](CHARMM_vs_buildH.gif) is an animated gif showing the difference between the hydrogens in a CHARMM structure vs those reconstructed by buildH.

Note that this validation was made in August 2019. You can retrieve the corresponding old version of buildH [here](https://github.com/patrickfuchs/buildH/tree/7cf8a331b1758abffd03ebb9737704dee3f12a88).
