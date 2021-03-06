In this directory, we demonstrate the validity of order parameter calculation with buildH, by comparing it with the main program used in NMRlipids (written by J. Melcr). The detailed protocol is described in `launch.sh`.

Practically, we launch buildH to reconstruct hydrogens on the Berger trajectory present in the parent directory.

Then we launch JMelcr Prog (`calcOrderParameters.py` retrieved on Jan 18 2020 from [here](https://raw.githubusercontent.com/NMRLipids/MATCH/master/scripts/calcOrderParameters.py) and converted to Python 3 by P. Fuchs) to calculate the order parameters on this traj with Hs.

Eventually, we compare the order parameters between the two programs. The results are plotted in `comp_JMelcrProg.pdf`. The mean order parameter difference averaged over all C-Hs is 0.0002. This very small difference is probably due to rounding errors, but it is by far negligible. The OP_stdev and OP_stem are also very close to builH.

In summary, the difference between both programs is negligible. It thus validates the part of buildH that calculates the order parameter.

