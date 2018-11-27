"""
This is a test file for the molecule function in the libgen module.
"""
import pytest

from chemlg.libgen import molecule

@pytest.fixture()
def smiles_code():
    return ('CcC', 'F1-F2')


def test_general(smiles_code):
    mol = molecule(*smiles_code)
    assert isinstance(mol, dict)
    assert mol['can_smiles'] == 'C[CH]C\t\n'
    assert mol['code'] == 'F1-F2'
    assert mol['smiles'] == 'CcC'
    assert mol['reverse_smiles'] == 'C[CH]C'

def test_exceptions(smiles_code):
    with pytest.raises(OSError):
        molecule('asdf', 'F1-F2')