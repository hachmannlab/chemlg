import pytest
from chemlg.libgen import if_add


def test_general():

    val = if_add( 'CC', {'include_bb': ['F2']}, 'F1-F3-F2')
    assert val

    val = if_add( 'CCCC', {'2': (0, 5)}, 'F1')
    assert val

    val = if_add( 'CC', {'3': (0, 2)}, 'F1')
    assert val

    val = if_add( 'CO', {'4': (25, 35)}, 'F1')
    assert val 

    val = if_add( 'C1CCC1', {'5': (0, 2)}, 'F1')
    assert val 

    val = if_add( 'c1ccc1', {'6': (0, 2)}, 'F1')
    assert val 

    val = if_add( 'C1CCC1', {'7': (0, 2)}, 'F1')
    assert val 

    val = if_add( 'CC', {'8': (0, 1)}, 'F1')
    assert val 

    val = if_add( 'CC=CC', {'9': (0, 2)}, 'F1')
    assert val 

    val = if_add( 'CC#CC', {'10': (0, 2)}, 'F1')
    assert val 

    val = if_add( 'CCOCOCCC', {'heteroatoms': [('O', 2)]}, 'F1')
    assert val

    val = if_add( 'CCCC', {'lipinski': True}, 'F1')
    # assert val

    val = if_add( 'CCCC', {'fingerprint': ['c1ccccc1', 0.1]}, 'F1')
    # assert val

    val = if_add( 'CCCC', {'14': ['CCCCCC']}, 'F1')
    # assert val

    val = if_add( 'CCCC', {'15': ['CCCCCCC']}, 'F1')
    # assert not val



def test_exceptions():

    val = if_add( '', {}, 'F1')
    assert not val

    val = if_add( 'smiles', {}, 'F1')
    assert not val

    val = if_add( 'CC', {'include_bb': ['F2']}, 'F1-F3')
    assert not val 

    val = if_add( 'CC', {'2': (10, 20)}, 'F1')
    assert not val

    val = if_add( 'CC', {'3': (10, 20)}, 'F1')
    assert not val 

    val = if_add( 'CCCCCC', {'4': (10, 20)}, 'F1')
    assert not val 

    val = if_add( 'C1CCC1', {'5': (10, 20)}, 'F1')
    assert not val

    val = if_add( 'C1CCC1', {'6': (10, 20)}, 'F1')
    assert not val 

    val = if_add( 'c1ccc1', {'7': (10, 20)}, 'F1')
    assert not val 

    val = if_add( 'CCCC', {'8': (40, 50)}, 'F1')
    assert not val 

    val = if_add( 'CC=CC', {'9': (40, 50)}, 'F1')
    assert not val 

    val = if_add( 'CC#CC', {'10': (40, 50)}, 'F1')
    assert not val 

    val = if_add( 'CCOCOCCCO', {'heteroatoms': [('O', 2)]}, 'F1')
    assert not val 

    val = if_add( 'CCCC', {'lipinski': True}, 'F1')
    # assert not val 

    val = if_add( 'CCCC', {'fingerprint': ['c1ccccc1', 0.1]}, 'F1')
    # assert not val 

    val = if_add( 'CCCC', {'14': ['c1ccccc1']}, 'F1')
    # assert not val 

    val = if_add( 'CCCC', {'15': ['c1ccccc1']}, 'F1')
    # assert val

    val = if_add( 'CCCC', {'bb_final_lib': False}, 'F1')
    assert not val 
