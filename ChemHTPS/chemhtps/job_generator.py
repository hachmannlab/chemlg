#!/usr/bin/env python

_MODULE_NAME = "job_generator"
_MODULE_VERSION = "v1.0.0"
_REVISION_DATE = "2015-06-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module generates the computational chemistry jobs."

# Version history timeline:
# v1.0.0 (2015-06-24): basic implementation

###################################################################################################
# TASKS OF THIS MODULE:
# -generate the computational chemistry jobs
# -prioritize the jobs pool
###################################################################################################

###################################################################################################
# TODO:
#
###################################################################################################

import sys
import os

from misc import (chk_mkdir)

###################################################################################################

def generate_jobs():
    """(generate_jobs):
        This function generates the computational chemistry jobs.
    """
    cwd = os.getcwd()
    geom = os.listdir(cwd + '/screeninglib/geometrylib')
    struct = os.listdir(cwd + '/screeninglib/structurelib')
    if geom == [] and struct == []:
        sys.exit('There is no candidate library of molecules in ' + cwd)

    # Not sure this is the best way to get the template should maybe ask which template to use
    template = 'PBE0_def2svp.inp'
    with open(cwd + '/job_templates/' + template, 'r', 0) as jt:
        job_template = jt.readlines()
    for geo in geom:
        job_file = cwd + '/screeninglib/geometrylib/' + geo
        job_dir = cwd + '/jobpool/' + geo.split('.')[0]
        chk_mkdir(job_dir)
        os.rename(job_file, job_dir + '/' + geo)  # maybe should use shutil.move in case multiple file systems
        job_template.append('* xyz 0 1 ' + geo)
        with open(job_dir + '/' + geo.split('.')[0] + '.inp', 'w+') as tmp:
            tmp.writelines(job_template)


def prioritize_pool():
    """(prioritize_pool):
        This function prioritizes the jobs pool.
    """

# TODO: write this function creating all the stuff in the dummy_project, including the config file containing the project name
