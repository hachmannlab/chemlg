import pytest
from chemlg.libgen import create_link, molecule


def test_general():
    mol1 = molecule('CCc1ccccc1', 'F1-F2-F5')
    mol2 = molecule('CC', 'F3')
    lib = create_link(mol1, mol2, {})
    assert isinstance(lib, list)
    assert len(lib) == 5
    
