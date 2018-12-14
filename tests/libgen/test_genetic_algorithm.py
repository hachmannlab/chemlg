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
                                bb_file=bb_loc,
                                config_file=cf_loc,
                                output_dir='./',
                                max_indi_len = 20,
                                weights=(-1.0, ), 
                                pop_size=30, 
                                n_generations=20)
    ga_test.fit()
    best_ind_df, best_individual = ga_test.search()
    assert best_ind_df is not None
    assert best_individual is not None
