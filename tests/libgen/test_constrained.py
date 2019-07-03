from chemlg.constrained import GeneticAlgorithm
import pytest
import pybel
import os
import pkg_resources
import tempfile
import shutil

def objective(mol_ob):
    mol2 = pybel.readstring("smi", "c1ccccc1")
    tanimoto = mol_ob.calcfp()|mol2.calcfp() 
    return tanimoto

@pytest.fixture()
def cf_path():
    return pkg_resources.resource_filename(
        'chemlg', os.path.join('templates', 'config.dat'))

@pytest.fixture()
def bb_path():
    return pkg_resources.resource_filename(
        'chemlg', os.path.join('templates', 'building_blocks.dat'))

@pytest.fixture()
def setup_teardown():
    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    # return test directory to save outputs
    yield test_dir
    # Remove the directory after the test
    shutil.rmtree(test_dir)

def test_genetic_algorithm(bb_path, cf_path, setup_teardown):
    bb_loc = bb_path
    cf_loc = cf_path
    ga_test = GeneticAlgorithm(evaluate=objective,
                                fitness = (('max', 0.01),),
                                bb_file=bb_loc,
                                config_file=cf_loc,
                                output_dir=setup_teardown,
                                crossover_size=5,
                                mutation_size=5,
                                algorithm=2)
    
    ga_test.search(n_generations=1)
    assert len(ga_test.population) >= 1
