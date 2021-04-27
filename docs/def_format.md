# Lipid def file format

**TODO**: to be completed!!!

`order_parameters_definitions_MODEL_X_Y.def`

This file is a mapping file created for the [NMRlipids](https://nmrlipids.blogspot.com/) project.

It aims at giving a unique name for each order parameter value along the lipid regardless the model of lipid used. This `.def` files can be found on the [MATCH repository](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs). Two examples of such files are present in directories `docs/Berger_POPC_test_case` and `docs/CHARMM36_POPC_validation`.

If an output trajecotry (option `-opx`) is requested, this `.def` file **must** contain all possible pairs of C-H to reconstruct (since the whole trajectory with Hs will be reconstructed). This option is slow, we do not recommend it if an output xtc file is not wanted.

If no option `-opx` is used buildH uses fast indexing. In this case the `.def` file can contain any subset of all possible C-H pairs. For example, if one wants to get OPs for the polar head only (Berger POPC), the `.def` file could look like the following:

```
beta1 POPC C5  H51
beta2 POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1 POPC C12 H121
g3_2 POPC C12 H122
g2_1 POPC C13 H131
g1_1 POPC C32 H321
g1_2 POPC C32 H322
```
