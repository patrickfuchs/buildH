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
  orcid: 
  affiliation: 1
- name: Amélie Bacle
  orcid: 
  affiliation: 2
- name: Pierre Poulain
  orcid: 0000-0003-4177-3619
  affiliation: 3
- name: Patrick F.J. Fuchs^[corresponding author]
  orcid: 0000-0001-7117-994X
  affiliation: "4, 5"
affiliations:
- Laboratoire de Biochimie Théorique (LBT), CNRS, F-75005 Paris, France
  index: 1
- Laboratoire Coopératif "Lipotoxicity and Channelopathies - ConicMeds", Université de Poitiers, F-86000 Poitiers, France
  index: 2
- name: Mitochondria, Metals and Oxidative Stress group, Institut Jacques Monod, UMR 7592, Université de Paris, CNRS, F-75013 Paris, France.
  index: 3
- Sorbonne Université, Ecole Normale Supérieure, PSL Research University, CNRS, Laboratoire des Biomolécules (LBM), F-75005 Paris, France
  index: 4
- Université de Paris, UFR Sciences du Vivant, F-75013 Paris, France
  index: 5
date: 27 May 2021
bibliography: paper.bib
---

# Statement of need

Molecular dynamics (MD) simulations of lipids are widely used to understand the complex structure and dynamics of biological or model membranes [Tieleman1997; @Feller2000; @Lyubartsev2011]. They are very complementary to biophysical experiments and have thus become important to get insights at the microscopic scale. Many of them are performed using all-atom (AA) or united-atom (UA) representations. In AA force fields (such as CHARMM36 [@Klauda2010]), all the atoms are considered, whereas in UA force fields (such as Berger [@Berger1997]) the aliphatic hydrogen atoms (Hs) are merged to their parent carbon into a larger particule representing a CH, CH2 or CH3 (e.g. a methyl group is represented by a single CH3 particle). In simulations of phospholipids, the use of UA representations allows to divide by almost 3 the number of atoms to simulate because these molecules contain many aliphatic Hs.
MD simulations of lipids are usually validated against experimental data [@Klauda2010] or used to help interpret experiments [@Feller2007]. One type of experiment which is widely for that is $^2H$ NMR. In this type of experiment, aliphatic Hs atoms are replaced by deuterons. $^2H$ NMR allows one to measure the order parameter of a given C-H bond (where the H is replaced by a deuteron):

$$S_{CH} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the C-H bond and the magnetic field. The symbol $\langle ... \rangle$ means averaging over time and molecules. This order parameter is useful because it is directly related to the flexibility of the given C-H bond. In fact, it describes the amount of possible orientations visited by the C-H bond. It varies between 1 and -0.5 and is unitless. Values close to 0 mean high mobility, when it goes away from 0 (towards negative or positive values), it means the C-H bond gets less mobile. Although it varies theoretically between 1 and -0.5, its absolute value $\lvert S_{CH} \rvert$ is often reported because its sign is usually difficult to measure experimentally [@Ollila2014]. In MD simulations, $\theta$ is the angle between the C-H bond and the normal to the membrane (usually the $z$ axis). For AA simulations, $S_{CH}$ is trivial to calculate. However, it is more difficult for UA simulations since the aliphatic Hs are not present. There are two main strategies to compute $S_{CH}$ from UA simulations [@Piggot2017]: i) expressing $S_{CH}$ as a function of the coordinates of other atoms (for example the neighboring carbon atoms $i-1$ and $i+1$ of carbon $i$ [@Douliez1995], ii) reconstructing Hs coordinates and calculate $S_{CH}$ as in AA simulations. The trend in the last years was more to use strategy ii), such as in the NMRlipids project [@Botan2015; @Catte2016; @Antilla2019]. NMRlipids is an open science project comparing MD simulations to experimental $S_{CH}$ with the goal of improving force fields for lipids. 
It is pretty easy to reconstruct Hs from the heavy atoms using standard geometric rules respecting stereochemistry. In the first NMRlipids project [@Botan2015], H reconstruction was performed using a program called `g_protonate` from the software GROMACS version 3.* [@Berendsen1995]. However, `g_protonate`  has been removed from GROMACS version 4.* or higher. So there was a need to find other solutions. Currently, there are programs in the field of chemoinformatics that are able to build Hs such as OpenBabel [@O'Boyle2011] or other proprietary softwares. It is also doable to use `pdb2gmx` from the GROMACS software [@Abraham2015] to build Hs, but it is not initially intended for that, its purpose is to build a topology. Many of these solutions remain tricks where one uses a complicated software for just doing one easy task. Using these open or proprietary softwares usually has a slow learning curve. Moreover, it is not always easy to use them in the Unix command line when they have complicated graphical interfaces, thus complicating their use in automatic pipelines of trajectory analyses.
Here, we propose a software `buildH` to fill this need. `buildH` is very light and usable in the Unix command line or as a Python module making it a tool of choice to integrate in pipeline analyses or Jupyter notebooks. `buildH` has been and is currently widely used in the NMRlipids project [@Botan2015; @Catte2016; @Antilla2019].

# Overview

`buildH` is a Python software that enables automatic analyses of order parameter calculations from UA trajectories of lipids. The software has the following features:

- It reads a single structure or a trajectory of lipids in a UA representation.
- It reconstructs the aliphatic Hs using standard geometric rules.
- From the reconstructed Hs, it calculates and outputs the order parameters on each requested C-H bond.
- Optionnaly, it outputs a structure (in pdb format) and a trajectory (in xtc format) with all reconstructed Hs.

Beyond order parameter calculations, the trajectory with Hs can be be used for any further analyses (e.g. precise molecular volume calculation).

`buildH` has been natively developed for a use in the Unix command line. It possesses a minimum number of options so that it is very easy to use and to learn. It is also possible to use it as a Python module which may be convenient in some cases.

To reconstruct H atoms, `buildH` uses standard geometric rules. These rules require so-called *helper* atoms. For example, the reconstruction of the two Hs of a CH2 on carbon $C_i$, requires two helpers which are $C_{i-1}$ and $C_{i+1}$, that is, the two neighbours of $C_i$ in the chain (note that helpers can also be other heavy atoms such as oxygen or nitrogen). The list of helpers used for the reconstruction of each H is written in a json file. Many json files are already present on the `buildH` repository representing the major lipids: Phoshphatidylcholine (PC), Phoshphatidylethanolamine (PE), Phoshphatidylglycerol (PG) for the polar heads, palmitoyl, myristoyl, oleoyl for the aliphatic chains, as well as cholesterol. Major UA force fields are also represented (Berger [@Berger1997], GROMOS-CPK [@Piggot2012], CHARMM-UA [@Lee2014]). In case a user wants to analyze a lipid which is not present in `buildH`, it is possible to supply his/her own json file with a full documentation on how to do that.

All structure and trajectory reading / writing are dealt with the MDAnalysis module [Michaud-Agrawal2011; Gowers2016]. Mathematical vector operations are dealt with Numpy [@Harris2020] and accelerated with Numba [@Lam2015]. This allows to get very decent performances although the software is written in pure Python. For example, the reconstruction of all Hs and order parameter calculation on a trajectory of 2500 frames with 128 POPC molecules can be handled in approximately 7 minutes using a single core Xeon processor at 3.60 GHz, whereas it was almost 30 minutes long without Numba.

`buildH` has been implemented with good practices of software development in mind [@Jimenez2017; @Taschuk2017]:

- version control repository on GitHub (https://github.com/patrickfuchs/buildH),
- open-source license (BSD-3-Clause),
- continuous integration through tests,
- and documentation (https://buildh.readthedocs.io/).

Some notebooks are present on the the github repository which explain how `buildH` works and how to analyze data produced by it. In case of trouble, any user can post an issue.

`buildH` is available in the Python Package Index (PyPI) as well as in the Bioconda repository. All versions of the software are archived in the Zenodo repository (https://zenodo.org/record/4792443) and in the Software Heritage archive (???).

# Acknowledgements

The authors thank the community of [NMRlipids](http://nmrlipids.blogspot.com/) for useful discussions, especially Samuli Ollila.

# References

