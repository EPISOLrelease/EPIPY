# Epipy: A Python Package for 3DRISM Solvation Calculations of Chemical and Biological Molecules
>Please citate the following papers for Epipy:
>1. "A Python Tutorial for 3DRISM Solvation Calculations of Chemical and Biological Molecules", Swanson, P., Cao, S., Huang, X, https://chemrxiv.org/engage/chemrxiv/article-details/68a6903c728bf9025e6c91ed
>2. “EPISOL: A Software Package with Expanded Functions to Perform 3D-RISM Calculations for the Solvation of Chemical and Biological Molecules“, Cao, S.; Kalin, M.L.; Huang, X., J. Comput. Chem., 44, 1536-1549, (2023)

## Abstract

The 3-Dimensional Reference Interaction Site Model (3DRISM) provides a powerful grid-based solvation model for chemical and biological solutes, which balances the calculation accuracy and efficiency. We previously developed EPISOL (Expanded Package of Integral Equation Theory-Based Solvation) to enable efficient 3DRISM calculations. EPISOL implements 22 different closures and several variations of 3DRISM. EPISOL is compatible with both AMBER and GROMACS simulation packages. The original EPISOL was written in C++ and includes a kernel library that allows integration of EPISOL routines into other software. In this work, we introduce EPIPY, a Python-based package that leverages the EPISOL kernel library to streamline 3DRISM calculations. We first provide an overview of 3DRISM, then present a step-by-step tutorial on running and analyzing 3DRISM calculations using EPIPY. Our tutorial examples demonstrate how to generate water distributions around chemical compounds and ions, identify the most probable water coordinates near a protein, and compute the solvation free energies of small organic molecules. 



> [!WARNING]
> You must have the EPISOL kernel downloaded and installed ( see [kernel install instrcutions](./EPISOL_kernel_install_instrcutions.md) for more)
> 
> In the future this will be included with pip, however currently you must do a 2-step installation

Includes:
* setting up and running 3DRISM commands
* Placement of waters 
* post-processing and selection tools  
>[!TIP]
>see wiki for more

For installation with PIP run 
```
pip install episol
```
epipy is contained within the episol module
```python
from episol import epipy
```
---

# tutorials
We have created several jupyter notebook tutorials exhibiting some of the many tools available in EPIPY. The tutorials come in several flavors, notebooks with topology file generating capabilities, those which use pre-written topologies (and require no 3rd party packages), and those for offline non-colab use. The notebooks with topology generating capabilites will work for essentailly any system and we recomend the users to upload their own molcules/proteins. 

## Google colab - with topology generation  
Contains third-party packages to generate topologies. This will work for _essentially_ anything you throw at it so feel free to get creative.
* Introduction tutorial and walkthrough of methane and nitrous ion [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/episol_methane_and_nitrous_tutorial_release.ipynb)
  
* Protein calculation, and water placement [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/episol_protein_and_water_placement_release.ipynb)

* High throughput small molecule generation and solvation free energy calculation [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/episol_small_molecule_high_throughput_release.ipynb)


## Google colab - no external packages (aside from epipy) 
for the 'plain' tutorials we use pre-crafted PDBs and topology files.

* Introduction tutorial and walkthrough of methane and nitrous ion [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/plain_methane_and_nitrous_tutorial_notebook.ipynb)
  
* Protein calculation, and water placement [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/plain_protein_and_water_placement_notebook.ipynb)

* High throughput small molecule generation and solvation free energy calculation  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/EPISOLrelease/EPIPY/blob/main/Colab_tutorials/plain_episol_small_molecule_high_throughput.ipynb)

---
Notebooks for local use (and their corresponding environment.yml file) are found in the [tutorial folder](./tutorials). 

> [!WARNING]
> epipy is in its beta phase and the code will be continually improved
>
> thank you for your feedback :)
