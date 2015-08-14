#!/usr/bin/env python

_MODULE_NAME = "project_setup"
_MODULE_VERSION = "v1.0.0"
_REVISION_DATE = "2015-06-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module sets up the project, including the file and directory structure."

# Version history timeline:
# v1.0.0 (2015-06-24): basic implementation

###################################################################################################
# TASKS OF THIS MODULE:
# -set up the project, including the file and directory structure
###################################################################################################

###################################################################################################
# TODO:
# 
###################################################################################################

import os

###################################################################################################

def setup_project(project_name):
    """(setup_project):
        This function sets up the project, including the file and directory structure.
    """
    dir_list = ['/archive', '/db', '/jobpool', '/job_templates', '/lost+found', '/screeninglib/geometrylib',
                '/screeninglib/structurelib']
    cwd = os.getcwd()  # just in case we need this

    for dir in dir_list:
        os.makedirs(project_name + dir, 0755)

    with open(project_name + '/' + project_name + '.config', 'w') as config:
        config.write('project_name = ' + project_name)

    with open(project_name + '/screeninglib/buildingblocks.dat', 'w') as build:
        build.close()

    with open(project_name + '/screeninglib/generatorrules.dat', 'w') as gener:
        gener.close()
