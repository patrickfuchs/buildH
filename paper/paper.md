---
title: 'buildH: Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters'
tags:
- python
- molecular-dynamics-simulation
- order-parameters
- lipids
- united-atom
authors:
- name: Hubert Santuz
  orcid: 0000-0001-6149-9480
  affiliation: "1, 2"
- name: Amélie Bacle
  orcid: 0000-0002-3317-9110
  affiliation: 3
- name: Pierre Poulain
  orcid: 0000-0003-4177-3619
  affiliation: 4
- name: Patrick F.J. Fuchs^[corresponding author]
  orcid: 0000-0001-7117-994X
  affiliation: "5, 6"
affiliations:
- name: CNRS, Université de Paris, UPR 9080, Laboratoire de Biochimie Théorique, 13 Rue Pierre et Marie Curie, F-75005 Paris, France
  index: 1
- name: Institut de Biologie Physico-Chimique–Fondation Edmond de Rothschild, PSL Research University, Paris, France
  index: 2
- name: Laboratoire Coopératif "Lipotoxicity and Channelopathies - ConicMeds", Université de Poitiers, F-86000 Poitiers, France
  index: 3
- name: Université de Paris, CNRS, Institut Jacques Monod, F-75006, Paris, France
  index: 4
- name: Sorbonne Université, Ecole Normale Supérieure, PSL Research University, CNRS, Laboratoire des Biomolécules (LBM), F-75005 Paris, France
  index: 5
- name: Université de Paris, UFR Sciences du Vivant, F-75013 Paris, France
  index: 6
date: 27 May 2021
bibliography: paper.bib
---

# Background

Molecular dynamics (MD) simulations of lipids are widely used to understand the complex structure and dynamics of biological or model membranes [@Tieleman1997; @Feller2000; @Lyubartsev2011]. They are very complementary to biophysical experiments and have thus become important to get insights at the microscopic scale. Many of them are performed using all-atom (AA) or united-atom (UA) representations. In AA force fields (such as CHARMM36 [@Klauda2010]), all the atoms are considered whereas in UA force fields (such as Berger [@Berger1997]) the aliphatic hydrogen atoms (Hs) are merged with their parent carbon into a larger particle representing a CH, CH2 or CH3 (e.g. a methyl group is represented by a single CH3 particle). In simulations of phospholipids, the use of UA representations allows one to reduce the number of particles that are simulated to almost a third of the true number of atoms. This is because lipid molecules contain many aliphatic Hs. This simplification thus reduces the computational cost without losing important chemical details.

MD simulations of lipids are usually validated against experimental data [@Klauda2010] or used to help interpret experiments [@Feller2007]. One type of experiment which is often used for that is $^2H$ NMR. In this type of experiment, aliphatic Hs are replaced by deuterons. $^2H$ NMR allows one to measure the order parameter of a given C-H bond (where the H is replaced by a deuteron):

$$S_{CH} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the C-H bond and the magnetic field. The symbol $\langle ... \rangle$ means averaging over time and molecules.

This order parameter is useful because it is directly related to the flexibility of the given C-H bond. It describes the amount of possible orientations visited by the C-H bond, ranges between -0.5 and 1 (inclusive), and is unitless. Values close to 0 indicate high mobility, while values with higher magnitudes (either positive or negative) indicate that the C-H bond is less mobile. Although it varies theoretically between -0.5 and 1, its absolute value $\lvert S_{CH} \rvert$ is often reported because its sign is usually difficult to measure experimentally [@Ollila2014]. In MD simulations, $\theta$ is the angle between the C-H bond and the normal to the membrane (usually the $z$ axis). For AA simulations, $S_{CH}$ is trivial to calculate. However, it is more difficult for UA simulations since the aliphatic Hs are not present. There are two main strategies to compute $S_{CH}$ from UA simulations [@Piggot2017]: i) expressing $S_{CH}$ as a function of the coordinates of other atoms [@Douliez1995], ii) reconstructing hydrogen coordinates and calculating $S_{CH}$ as in AA simulations. The trend in the recent years has been towards strategy ii), such as in the NMRlipids project [@Botan2015; @Catte2016; @Antila2019; @Bacle2021]. NMRlipids is an open science project which uses MD simulations and experimental $S_{CH}$ with the goal of improving lipid force fields or conducting fundamental research on the structure and dynamics of lipid membranes.

# Statement of Need

Reconstructing Hs from the heavy atom coordinates can be done using standard geometric rules respecting stereochemistry. In the first NMRlipids project [@Botan2015], H reconstruction was performed using a program called `g_protonate` from the software GROMACS version 3.* [@Berendsen1995]. However, `g_protonate` has been removed from GROMACS version 4.* or higher. So, there was a need to find other solutions. Currently, there are many programs in the field of chemoinformatics that are able to build Hs such as OpenBabel [@OBoyle2011] or other proprietary softwares. It is also possible to use `pdb2gmx` from the GROMACS software [@Abraham2015] to build Hs, but it is not initially intended for that since its main purpose is to build a topology. Many of these solutions remain workarounds where one uses a sophisticated software for just doing one basic task. Using these open or proprietary software usually has a slow learning curve. Moreover, it is not always easy to use them in the Unix command line when they have complicated graphical interfaces, thus complicating their use in automatic pipelines of trajectory analyses.
Here, we propose the `buildH` software to meet this need. `buildH` is very light and usable in the Unix command line or as a Python module, making it a tool of choice to integrate in analysis pipelines or Jupyter notebooks. It can be easily extended or customized if needed. `buildH` is currently widely used in the NMRlipids project IVb dealing with the conformational plasticity of lipid headgroups in cellular membranes and protein-lipid complexes [@Bacle2021]. In addition, it is planned to use `buildH` in the next NMRlipids project dealing with a databank containing MD trajectories of lipids [@Kiirikki2021].

# Overview

`buildH` is a Python software that enables automatic analyses of order parameter calculations from UA trajectories of lipids. The software has the following features:

- It reads a single structure or a trajectory of lipids in a UA representation.
- It reconstructs the aliphatic Hs using standard geometric rules.
- From the reconstructed Hs, it calculates and outputs the order parameters on each requested C-H bond.
- Optionally, it outputs a structure (in pdb format) and a trajectory (in xtc format) with all reconstructed Hs.

Beyond order parameter calculations, the trajectory with Hs can be used for any further analyses (e.g. precise molecular volume calculation).

`buildH` has been natively developed for a use in the Unix command line. It possesses a minimum number of options, making it easy to use. It is also possible to utilize it as a Python module, which may be convenient in some cases.

To reconstruct H atoms, `buildH` uses standard geometric rules. These rules require so-called *helper* atoms. For example, the reconstruction of the two Hs of a CH2 on carbon $C_i$, requires two helpers which are $C_{i-1}$ and $C_{i+1}$, that is, the two neighbors of $C_i$ in the chain (note that helpers can also be other heavy atoms such as oxygen or nitrogen). The list of helpers used for the reconstruction of each H is written in a json file. Many json files are already present on the `buildH` repository representing the major lipids: Phoshphatidylcholine (PC), Phoshphatidylethanolamine (PE), Phoshphatidylglycerol (PG) for the polar heads and palmitoyl, myristoyl, oleoyl for the aliphatic chains, as well as cholesterol. Major UA force fields are also represented (Berger [@Berger1997], GROMOS-CKP [@Piggot2012], CHARMM-UA [@Lee2014]). In case a user wants to analyze a lipid which is not present in `buildH`, a step-by-step documentation guides the user in the process of creating and supplying his/her lipid description as a json file.

All structure and trajectory read / write operations are handled by the MDAnalysis module [@Michaud-Agrawal2011; @Gowers2016]. Mathematical vector operations are performed by Numpy [@Harris2020] and accelerated with Numba [@Lam2015], leading to very decent performances. For example, the reconstruction of all Hs and order parameter calculation on a trajectory of 2500 frames with 128 POPC molecules can be handled in approximately 7 minutes using a single core Xeon processor at 3.60 GHz, whereas it was almost 30 minutes long without Numba.

`buildH` has been implemented with good practices of software development in mind [@Jimenez2017; @Taschuk2017]:

- version controlled repository on GitHub [https://github.com/patrickfuchs/buildH](https://github.com/patrickfuchs/buildH),
- open-source license (BSD-3-Clause),
- continuous integration through tests,
- and documentation [https://buildh.readthedocs.io/](https://buildh.readthedocs.io/).

Some notebooks are provided in the GitHub repository to explain how `buildH` works and how to analyze the data produced. In case of trouble, any user can post an issue on GitHub.

`buildH` is available in the Python Package Index (PyPI) as well as in the Bioconda repository. The current version 1.5.0 of `buildH` is archived in the Zenodo repository ([https://zenodo.org/record/5080126](https://zenodo.org/record/5080126)) and in the Software Heritage archive ([swh:1:dir:eb46a03f6e6188f93bd0b39ab78b4640e777ecd1](https://archive.softwareheritage.org/swh:1:dir:eb46a03f6e6188f93bd0b39ab78b4640e777ecd1;origin=https://github.com/patrickfuchs/buildH;visit=swh:1:snp:ce3e8efd06ea951364d1f5ee26cbd5776c117dc3;anchor=swh:1:rel:980b3a2f1332d172230eaf74a05bbfc9145b1266)).

# Acknowledgements

The authors thank the community of [NMRlipids](http://nmrlipids.blogspot.com/) for useful discussions, especially Samuli Ollila.

# References

