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
import shutil
import curses

from misc import (chk_mkdir)
from template_generator import generate_template, runmenu

###################################################################################################

def generate_jobs():
    """(generate_jobs):
        This function generates the computational chemistry jobs.
    """
    cwd = os.getcwd()

    menu = curses.initscr()
    curses.noecho()
    curses.cbreak()
    curses.curs_set(0)
    curses.start_color()
    menu.keypad(1)

    libraries = ['Available libraries', 'Please choose a library: ', 'Exit Program']
    files = []
    for root, directories, filenames in os.walk(cwd + '/screeninglib'):
        for directory in directories:
            if os.listdir(os.path.join(root, directory)):
                files.append(os.listdir(os.path.join(root, directory))[0:10])
                libraries.insert(-1, os.path.join(root, directory))

    pos = runmenu(menu, libraries)
    library = libraries[pos]
    lib = os.listdir(library)

    curses.nocbreak()
    curses.echo()
    curses.curs_set(1)
    menu.keypad(0)
    curses.endwin()

    # Not sure this is the best way to get the template should maybe ask which template to use
    template = generate_template()
    with open(template, 'r', 0) as jt:
        job_template = jt.readlines()
    job_template = tuple(job_template)
    for geo in lib:
        temp = list(job_template)
        job_file = library + '/' + geo
        job_dir = cwd + '/jobpool/short/' + geo.split('.')[0]
        chk_mkdir(job_dir)
        shutil.copy(job_file, job_dir + '/' + geo)
        temp.append('* xyzfile 0 1 ' + geo + '\n')
        with open(job_dir + '/' + geo.split('.')[0] + '.inp', 'w') as tmp:
            tmp.writelines(temp)


def prioritize_pool():
    """(prioritize_pool):
        This function prioritizes the jobs pool.
    """

# TODO: write this function creating all the stuff in the dummy_project, including the config file containing the project name
