Additional module requirements
==============================

The following additional packages should be installed in the environment.

.. code:: bash

    conda install –c conda-forge rdkit
    conda install ipykernel
    conda install -c conda-forge ipywidgets
    python -m ipykernel install --name my_chemlg_env --display-name "chemlg"

To build the chemlg input files, launch jupyter notebook and enter the following:

..code:: bash

    from chemlg.templates.main import config_builder


Building Blocks File:

User has an option to either upload an existing file or create a new one. For an existing file, all SMILES should start on a new line. To create a new file, enter the SMILES individually or as comma-separated values. The SMILES can be verified to correspond to the correct structure from the visualize option. The 'add building blocks' option can be used repeatedly until all the SMILES are added following which the 'create building blocks' file option creates the required building blocks file in the specified directory.

Generation Rules:

For creating the config file, follow the on-screen instructions and for more details, refer to:
 
.. toctree::

    chemlg.inputs