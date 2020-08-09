import pytest
from chemlg.libgen import lipinski
from openbabel import pybel

def test_lipinski():
    mol = pybel.readstring("smi", 'CCC')
    descr = lipinski(mol)
    assert isinstance(descr, dict)

