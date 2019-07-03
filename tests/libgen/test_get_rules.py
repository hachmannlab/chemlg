import os
import pkg_resources
import shutil
import tempfile
import pytest
import pandas as pd

from chemlg.libgen import get_rules

@pytest.fixture()
def data_path():
    config_path = pkg_resources.resource_filename('chemlg', os.path.join('templates', 'config.dat'))
    bb_path = pkg_resources.resource_filename('chemlg', os.path.join('templates', 'building_blocks.dat'))
    return config_path, bb_path

@pytest.fixture()
def setup_teardown():
    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    # return test directory to save outputs
    yield test_dir
    # Remove the directory after the test
    shutil.rmtree(test_dir)


def test_get_rules(data_path, setup_teardown):
    output_dir = setup_teardown
    config_path, bb_path = data_path
    config_file = open(config_path)
    rules_dict, args = get_rules(config_file, output_dir, {'rank': 0})
    assert rules_dict['include_bb'] == ['C']


