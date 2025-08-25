# Epipy: A Python Package for 3DRISM Solvation Calculations of Chemical and Biological Molecules
> Citation: "A Python Tutorial for 3DRISM Solvation Calculations of Chemical and Biological Molecules", Swanson, P., Cao, S., Huang, X, https://chemrxiv.org/engage/chemrxiv/article-details/68a6903c728bf9025e6c91ed

The 3-Dimensional Reference Interaction Site Model (3DRISM) provides a powerful grid-based solvation model for chemical and biological solutes, which balances the calculation accuracy and efficiency. We previously developed EPISOL (Expanded Package of Integral Equation Theory-Based Solvation) to enable efficient 3DRISM calculations. EPISOL implements 22 different closures and several variations of 3DRISM. EPISOL is compatible with both AMBER and GROMACS simulation packages. The original EPISOL was written in C++ and includes a kernel library that allows integration of EPISOL routines into other software. In this work, we introduce EPIPY, a Python-based package that leverages the EPISOL kernel library to streamline 3DRISM calculations. We first provide an overview of 3DRISM, then present a step-by-step tutorial on running and analyzing 3DRISM calculations using EPIPY. Our tutorial examples demonstrate how to generate water distributions around chemical compounds and ions, identify the most probable water coordinates near a protein, and compute the solvation free energies of small organic molecules. The EPIPY source code and accompanying tutorial files are available at: https://github.com/EPISOLrelease/EPIPY.


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

The code in the tutorials will work for _essentially_ anything you throw at it
so feel free to get creative.

# Google colab tutorials
* [Introduction tutorial and walkthrough of methane and nitrous ion](https://colab.research.google.com/drive/1EUE9_B168XrVjyvnOackqDtIH13vcLl-#offline=true&sandboxMode=true)
* [For protein calculation, and water placement](https://colab.research.google.com/drive/1uGGPi8CAEIJNg1_mkDwbyBzpy7WiqbnQ#offline=true&sandboxMode=true)
* [For high throughput small molecule generation and solvation free energy calculation](https://colab.research.google.com/drive/175uh_Fh0YWnoxpCPDr3Qav_mQCQ0YT0Z#offline=true&sandboxMode=true)

  
> [!WARNING]
> epipy is in its beta phase and the code will be continually improved
>
> thank you for your feedback :)
