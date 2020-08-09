from chemlg.libgen import library_generator
import pytest
from openbabel import pybel
import os
import pkg_resources
import tempfile
import shutil
import random

def objective(mol_ob):
    mol2 = pybel.readstring("smi", "c1ccccc1")
    tanimoto = mol_ob.calcfp()|mol2.calcfp() 
    return random.random()

@pytest.fixture()
def setup_teardown():
    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    # return test directory to save outputs
    yield test_dir
    # Remove the directory after the test
    shutil.rmtree(test_dir)

def test_genetic_algorithm(setup_teardown):
    ga_test = library_generator(config_file='tests/libgen/input_files/config.dat', building_blocks_file='tests/libgen/input_files/building_blocks.dat', output_dir=setup_teardown, genetic_algorithm_config='tests/libgen/input_files/genetic_algorithm_config.dat', cost_function=objective)
    
    assert len(ga_test.population) >= 10

def test_batch(setup_teardown):
    ga_test = library_generator(config_file='tests/libgen/input_files/config.dat', building_blocks_file='tests/libgen/input_files/building_blocks.dat', output_dir=setup_teardown, genetic_algorithm_config='tests/libgen/input_files/ga_batch.dat', fitnesses_list=[])

    assert ga_test.population is None   # since batch mode does not update the self.population variable
