# Quantum block Krylov subspace projector algorithm

Release 1.0.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17552224.svg)](https://doi.org/10.5281/zenodo.17552224)

## About 
This repository contains the code utilized in the article **Quantum block Krylov subspace projector algorithm for computing low-lying eigenenergies** by M. G. J. Oliveira and N. Glaser in PRA, 2025. Therefore, if you use our code, please adequately cite the article.

Contributors: M. G. J. Oliveira and  N. Glaser, NNF Quantum Computing Programme, Niels Bohr Institute, University of Copenhagen, Denmark

## Repository structure
This repository is organized as follows:

* Models:
    - *Heisenberg.py* - Contains the code to generate the Heisenberg model;
    - *molecules.py* - Contains the code to generate diatomic molecules operators;

* *data_generator.py* - Contains the code to generate all data needed to reproduce the plots in the article;
* *utils.py* - Contains some auxiliary functions;
* *qbksp.py* - Contains the quantum block Krylov subspace projector algorithm optimized for a real and symmetric Hamiltonian and real initial references;
* *Heisenberg.ipynb* - Contains the code to reproduce the Heisenberg model plots;
* *molecules.ipynb* - Contains the code to reproduce the molecular plots.

-------------------------------------------------------------------------------------------

## How to use this code

Please use **Python version 3.12.3** and the package requirements as listed in *requirements.txt*.
These can be installed via `pip install -r requirements.txt`.

To generate the necessary data, you should run the script *data_generator.py*. Then, you can use each notebook to plot and analyze the data. 

In all files, please replace the variable `Data_path` with the location where you would like to store the generated data.
