[![Build Status](https://travis-ci.org/hachmannlab/chemlg.svg?branch=master)](https://travis-ci.org/hachmannlab/chemlg)
[![codecov](https://codecov.io/gh/hachmannlab/chemlg/branch/master/graph/badge.svg)](https://codecov.io/gh/hachmannlab/chemlg)
# ChemLG
ChemLG is a smart and massive parallel molecular library generator for chemical and materials sciences.


## Code Design:
ChemLG is developed in the Python 3 programming language and uses OpenBabel and its Python extension, Pybel for handling molecules. The development follows a strictly modular and object-oriented design to make the overall code as flexible and versatile as possible.

## Installation and Dependencies:
The dependencies for ChemLG are OpenBabel and MPI4Py. It is recommended that these two dependencies are installed in a virtual environment prior to installing ChemLG. They can be installed via the conda installer:


    conda create --name my_chemlg_env python=3.6
    source activate my_chemlg_env
    conda install -c openbabel openbabel
    conda install -c anaconda mpi4py
    
You can download ChemML from Python Package Index (PyPI) via pip. 

    pip install chemlg


You can test the installation with:

    pytest -v


Following additional (optional) packages are required (within the environment) to unlock some of the other functionalities of ChemLG:

ChemLG graphical config file builder:

    conda install –c conda-forge rdkit
    conda install ipykernel
    conda install -c conda-forge ipywidgets
    python -m ipykernel install --name my_chemlg_env --display-name "chemlg"


ChemLG genetic algorithm module:

    pip install deap





## Contributors:

- Mohammad Atif Faiz Afzal, CBE Department, SUNY Buffalo
- Gaurav Vishwakarma, CBE Department, SUNY Buffalo
- Janhavi Dudwadkar, CBE Department, SUNY Buffalo
- Mojtaba Haghighatlari, CBE Department, SUNY Buffalo
- Johannes Hachmann, CBE Department, SUNY Buffalo

- We encourage any contributions and feedback. Feel free to fork and make pull-request to the "development" branch.



## Citation:
Please cite the use of ChemLG as:


    Afzal, M. A. F., Hachmann J. (2018) "ChemLG: A smart and massively parallel molecular library generator" https://github.com/hachmannlab/chemlg


## License:
ChemLG is copyright (C) 2014-2018 Johannes Hachmann, all rights reserved.
ChemLG is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause).
