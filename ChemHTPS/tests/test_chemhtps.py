#!/usr/bin/env python

_MODULE_NAME = "test_chemhtps"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2015-08-18"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "These are the tests"

# Version history timeline:
# v0.1.0 (2015-08-18): basic implementation

###################################################################################################
# TASKS OF THIS MODULE:
# -Test different portion of the project to make sure it works
###################################################################################################

###################################################################################################
# TODO:
#
###################################################################################################

import pytest
import os
import shutil
from chemhtps.project_setup import setup_project
from chemhtps.job_generator import generate_jobs


def test_gonna_fail():
    assert True == False  # Going to fail here on line 4


def test_setup_project():
    setup_project('testing')
    dir_list = ['/archive', '/db', '/jobpool', '/job_templates', '/lost+found', '/screeninglib/geometrylib',
                '/screeninglib/structurelib']
    file_list = ['/testing.config', '/screeninglib/buildingblocks.dat', '/screeninglib/generatorrules.dat']
    for dir in dir_list:
        assert True == os.path.isdir('testing' + dir)
    for file in file_list:
        assert True == os.path.isfile('testing' + file)


def test_generate_jobs():
    template = ['!DFT PBE0 Def2-SVP smallPRINT PRINTBASIS PRINTGAP\n', '!opt']
    with open('testing/job_templates/PBE0_def2svp.inp', 'w') as test_template:
        test_template.writelines(template)
    testgeo = ['C 0.000000 0.000000 0.000000\n', 'C 0.000000 0.000000 1.400000\n', 'C 1.212436 0.000000 2.100000\n',
               'C 2.424871 0.000000 1.400000\n', 'C 2.424871 0.000000 0.000000\n', 'C 1.212436 0.000000 -0.700000\n',
               'H -0.943102 0.000000 1.944500\n', 'H 1.212436 0.000000 3.189000\n', 'H 3.367973 0.000000 1.944500\n',
               'H 3.367973 0.000000 -0.544500\n', 'H 1.212436 0.000000 -1.789000\n', 'H -0.943102 0.000000 -0.544500\n']
    with open('testing/screeninglib/geometrylib/testing.xyz', 'w') as test_geo:
        test_geo.writelines(testgeo)
    cwd = os.getcwd()
    os.chdir(cwd + '/testing')
    generate_jobs()
    os.chdir('..')
    assert True == os.path.isfile('testing/jobpool/short/testing/testing.xyz')
    assert True == os.path.isfile('testing/jobpool/short/testing/testing.inp')
    assert True == os.path.isfile('testing/screeninglib/geometrylib/testing.xyz')


def test_wrap_up():
    shutil.rmtree('testing')
    assert False == os.path.isdir('testing')

if __name__ == '__main__':
    pytest.main()
