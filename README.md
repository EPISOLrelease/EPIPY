# Epipy: Python interface for episol
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

  
> [!WARNING]
> epipy is in its beta phase and the code will be continually improved
>
> thank you for your feedback :)
