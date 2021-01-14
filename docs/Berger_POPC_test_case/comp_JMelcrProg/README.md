In this dir, we demonstrate the validity of order parameter calculation with buildH, by comparing with the main program used in NMRlipids (written by J. Melcr). The detailed protocol is desribed in `launch.sh`.

Practically, we launch buildH to reconstruct hydrogens on the Berger trajectory present in the parent directory.

Then we launch JMelcr Prog to calculate the order parameters on this traj with Hs.

Eventually, we compare the order parameters between the two programs. The results are plotted in `comp_JMelcrProg.pdf`. The mean difference averaged all order parameters is 0.0002. This very small difference is probably due to rounding errors, but it is by far completely negligible.

In summary, the difference between both programs is negligible. It thus validates the part of buildH that calculates the order parameter.

