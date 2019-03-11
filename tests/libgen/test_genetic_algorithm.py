from chemlg.genetic_algorithm import GeneticAlgorithm
import pytest
import pybel
import os
import pkg_resources


def objective(mol_ob):
    mol2 = pybel.readstring("smi", "c1ccccc1")
    tanimoto = mol_ob.calcfp()|mol2.calcfp() 
    return tanimoto

def cf_path():
    return pkg_resources.resource_filename(
        'chemlg', os.path.join('templates', 'config.dat'))


def bb_path():
    return pkg_resources.resource_filename(
        'chemlg', os.path.join('templates', 'building_blocks.dat'))


def test_genetic_algorithm():
    bb_loc = bb_path()
    cf_loc = cf_path()
    ga_test = GeneticAlgorithm(evaluate=objective,
                                fitness = (('max', 0.01),),
                                bb_file=bb_loc,
                                config_file=cf_loc,
                                output_dir='./output',
                                crossover_size=5,
                                mutation_size=5,
                                algorithm=2)
    
    ga_test.search(n_generations=20)
    print(len(ga_test.population))
    assert len(ga_test.population) >= 10
