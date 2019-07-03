.. ChemLG documentation master file, created by
   sphinx-quickstart on Wed Feb  6 14:34:02 2019.

ChemLG – A Library Generator for the Exploration and Enumeration of Chemical and Materials Space
==========================================================================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



ChemLG is a smart and massive parallel molecular library generator for chemical and materials sciences.


Github Repository: https://github.com/hachmannlab/chemlg

Program Version: 0.3

Release Date: July 4, 2019




With contributions by:
Janhavi Abhay Dudwadkar (UB): Jupyter GUI


Code Design
+++++++++++

ChemLG is developed in the Python 3 programming language and uses OpenBabel and its Python extension, Pybel for handling molecules. The development follows a strictly modular and object-oriented design to make the overall code as flexible and versatile as possible. ChemLG can be run on a single core or in parallel on multiple cores. For the parallel execution, MPI4Py is also required along with OpenBabel as dependencies of ChemLG. 

Installation and Dependencies
+++++++++++++++++++++++++++++

It is highly recommended that a virtual environment is used to run ChemLG. The virtual environment and ChemLG and its dependencies can be installed as:

.. code:: bash

    conda create --name my_chemlg_env python=3.6
    source activate my_chemlg_env
    conda install -c openbabel openbabel
    conda install -c anaconda mpi4py
    pip install chemlg


You can test the installation with:

.. code:: bash

    pytest -v



.. toctree::
   :maxdepth: 4
   :caption: ChemLG Guide

   chemlg



Citation
++++++++


Please cite ChemLG as follows:

- Afzal, M. A. F.; Vishwakarma, G.; Dudwadkar, J. A.;Haghighatlari, M.; Hachmann, J. ChemLG– A Library Generator for the Exploration and Enumeration of Chemical and Materials Spaces. 2019; https://github.com/hachmannlab/chemlg

- M.A.F. Afzal, G. Vishwakarma, J. Hachmann, ChemLG – ChemLG– A Library Generator for the Exploration and Enumeration of Chemical and Materials Spaces. Available from: https://hachmannlab.github.io/chemlg

- J. Hachmann, M.A.F. Afzal, M. Haghighatlari, Y. Pal, Building and Deploying a Cyberinfrastructure for the Data-Driven Design of Chemical Systems and the Exploration of Chemical Space, Mol. Simul. 44 (2018), 921-929. DOI: 10.1080/08927022.2018.1471692




Acknowledgement
+++++++++++++++

ChemLG is based upon work supported by the U.S. National Science Foundation under grant #OAC-1751161. It was also supported by start-up funds provided by UB's School of Engineering and Applied Science and UB's Department of Chemical and Biological Engineering, the New York State Center of Excellence in Materials Informatics through seed grant #1140384-8-75163, and the U.S. Department of Energy under grant #DE-SC0017193.



License and Copyright
+++++++++++++++++++++

ChemLG is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause).

ChemLG is copyright (C) 2015-2018 Johannes Hachmann and Mohammad Atif Faiz Afzal, all rights reserved. 
(C) 2015-2019 Johannes Hachmann, Mohammad Atif Faiz Afzal
University at Buffalo - The State University of New York (UB)

Contact: hachmann@buffalo.edu, m27@buffalo.edu
http://hachmannlab.cbe.buffalo.edu
