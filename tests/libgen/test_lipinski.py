import pytest
from chemlg.libgen import lipinski
import pybel

def test_lipinski():
    mol = pybel.readstring("smi", 'CCC')
    descr = lipinski(mol)
    assert isinstance(descr, dict)

