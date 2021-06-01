# Order parameters and statistics

The mean order parameter of bond $CH_j$ is calculated using the standard formula:

$$\overline{S_{CH_j}} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the $CH_j$ bond and the normal to the membrane (usually the *z* axis), <...> means averaging over molecules and frames. $S_{CH}$ can be measured by NMR which is useful to validate simulation results, as largely described in the [NMRlipids project](http://nmrlipids.blogspot.com).

The order parameter output of buildH (default name `OP_buildH.out`) looks like this:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
gamma1_1             POPC    C1    H11    0.01304  0.12090  0.01069
gamma1_2             POPC    C1    H12    0.00666  0.09279  0.00820
gamma1_3             POPC    C1    H13   -0.01531  0.09141  0.00808
[...]
```

Each line corresponds to a given CH. The 4 first columns contain the generic name, residue name, carbon and hydrogen names respectively. The other columns contains different statistics on order parameters (OP):

- `OP_mean` is the mean OP of bond $CH_j$ averaged over all lipids and all frames of the trajectory, we shall write it $\overline{S_{CH_j}}$.
- `OP_stddev` is the standard deviation of the OP, we shall write it $\sigma(S_{CH_j})$; first we average each OP of bond $CH_j$ (e.g. the C-H of beta1) of residue $i$ (i.e. lipid $i$) over the whole trajectory:

$$ \overline{S_{CH_j}(i)} = \frac{1}{nframes} \sum_{t=0}^{t=nframes} S_{CH_j}(i)(t) $$

where $nframes$ is the total number of frames, then we calculate the standard deviation of those means over all residues:

$$ \sigma(S_{CH_j}) =
\sqrt{
\frac{1}{nres} \sum_{i=1}^{i=nres} (\overline{S_{CH_j}(i)} - \overline{S_{CH_j}})^2
}$$

where $nres$ is the total number of residues (i.e. lipids).
- `OP_stem` is the standard error of the mean averaged in the same spirit, let's call it $err(S_{CH_j})$:

$$err(S_{CH_j}) = \frac{\sigma(S_{CH_j})}{\sqrt{nres}}$$
