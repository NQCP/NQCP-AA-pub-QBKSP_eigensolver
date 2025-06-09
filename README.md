# Quantum block Krylov subspace projector algorithm

## About 
This repository contains the code utilized in the article **Quantum block Krylov subspace projector algorithm for computing low-lying eigenenergies**. Therefore, if you use our code, please adequately cite the article.

Contributors: M. G. J. Oliveira and  N. A. Glaser, NNF Quantum Computing Programme, Niels Bohr Institute, University of Copenhagen, Denmark

## Repository structure
This repository is organized as follows:

* Models:
    - *Heisenberg.py* - Contains the code to generate the Heisenberg model;
    - *molecules.py* - Contains the code to generate diatomic molecules operators;

* *data_generator.py* - Contains the code to generate all data needed to reproduce the plots in the article;
* *utils.py* - Contains some auxiliary functions;
* *qbksp.py* - Contains the quantum block krylov subspace projector algorithm;
* *Heisenberg.ipynb* - Contains the code to reproduce the Heisenberg model plots;
* *molecules.ipynb* - Contains the code to reproduce the molecular plots.

-------------------------------------------------------------------------------------------

## How to use this code

Please use **Python version 3.12.3** and the package requirements as listed in *requirements.txt*.
These can be installed via `pip install -r requirements.txt`.

To generate the necessary data, you should run the script *data_generator.py*. Then, you can use each notebook to plot and analyze the data. 

Please replace the data path with your data path.
