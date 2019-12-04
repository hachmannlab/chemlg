import os
import pkg_resources
import shutil
import tempfile
import pytest
import pandas as pd

from chemlg.libgen import library_generator

@pytest.fixture()
def setup_teardown():
    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    # return test directory to save outputs
    yield test_dir
    # Remove the directory after the test
    shutil.rmtree(test_dir)


def test_library_generator(setup_teardown):
    output_dir = setup_teardown
    a = library_generator(config_file='tests/libgen/input_files/config.dat', building_blocks_file='tests/libgen/input_files/building_blocks.dat', output_dir=output_dir)
    assert a is None