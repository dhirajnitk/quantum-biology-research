# Superradiance in Tryptophan Architectures
* This repository contains the code used to generate Figs. 7-12 in the Quantum Biology Laboratory’s recent work: **Patwa, H., Babcock, N.S., Kurian, P. "Quantum-enhanced photoprotection in neuroprotein architectures emerges from collective light-matter interactions." *Frontiers in Physics* (2024).** For further details, please visit https://www.quantumbiolab.com.
* Contains code that calculates the complex eigenvalues and/or eigenvectors of a non-Hermitian Hamiltonian that describes the interaction of a network of tryptophan chromophores with the electromagnetic field.
* It includes a pipeline that can calculate the thermal average (over a given number of realizations/trials) of the quantum yield with varying thermal disorder and system size values.
* It can be easily generalized to calculate the thermal average of any observable quantity that depends on the eigenvalues or eigenvectors for varying system sizes and static disorder values.

## Workflow
### Pipeline Input
* **DIR_TO_SCAN**: For each architecture/structure size, a file must be created that contains the positions and transition dipole vectors of each tryptophan molecule in the architecture (call this the *Position-Dipole* file). All Position-Dipole files should be put in one directory, and the full directory path should be passed in.
* **TRP_PER_SPIRAL**: Number of tryptophan molecules for each "subunit." This number is only for plotting purposes, the results of the calculations don't depend on it. If the notion of a "subunit" is not well-defined for a particular type of structure, it can just be set to 1.
* **NOISE_VALUES**: List of static disorder (W; units of inverse cm) values to run for, separated by spaces. Include W=0 in the list. 
* **NUM_REALIZATIONS**: Number of trials/realizations for a fixed Psoition-Dipole file and static disorder value. 
### Pipeline Output (artifacts)
* The outputs will be a collection of plots, both in .PDF and .PNG format. Figs. 7-10 in our work are direct outputs of the pipeline.
### Pipeline Output (files/directories created)
* All eigenvalues can be stored on the machine on which the pipeline is run (depending on user inputs), and all non-Hermitian Hamiltonians are stored as NumPy files. NumPy (.npy extension) files are binary representations of NumPy arrays, and can be saved and loaded through NumPy functions save() and load().

to be continued ...
