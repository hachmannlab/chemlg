import os
import pkg_resources
import shutil
import tempfile
import pytest
import pandas as pd

from chemlg.libgen import library_generator

@pytest.fixture()
def data_path():
    return pkg_resources.resource_filename(
        'chemlg', os.path.join('templates', 'config.dat'))

@pytest.fixture()
def setup_teardown():
    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    # return test directory to save figures
    yield test_dir
    # Remove the directory after the test
    shutil.rmtree(test_dir)


def test_library_generator(data_path, setup_teardown):
    output_dir = setup_teardown
    library_generator(config_file=data_path, output_dir=output_dir)
    a = pd.read_csv(os.path.join(output_dir, 'final_library.csv'))
    # test properly what should be in the first line of file, length of the file , ...
    assert a is not None


