import pytest
from chemlg.libgen import create_fused, molecule


def test_substitutions():
    mol1 = molecule('C1CCCCC1', 'F1-F2')
    mol2 = molecule('N1CC(Cl)CC1', 'F3')
    lib = create_fused(mol1, mol2, {'bb_final_lib':False})
    assert isinstance(lib, list)
    assert len(lib) >= 4
    
def test_general():
    mol1 = molecule('CCc1ccccc1', 'F1-F2')
    mol2 = molecule('o1cccc1', 'F3')
    lib = create_fused(mol1, mol2, {'bb_final_lib':False})
    assert isinstance(lib, list)
    assert len(lib) >= 4
