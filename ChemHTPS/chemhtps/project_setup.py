#!/usr/bin/env python

_MODULE_NAME = "project_setup"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2016-02-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module sets up the project, including the file and directory structure."

# Version history timeline:
# v0.0.1 (2015-06-24): basic implementation
# v0.1.0 (2016-02-24): alpha version

###################################################################################################
# TASKS OF THIS MODULE:
# -set up the project, including the file and directory structure
###################################################################################################

###################################################################################################
# TODO:
# 
###################################################################################################

import os
import fnmatch


###################################################################################################

def setup_project(project_name):
    """(setup_project):
        This function sets up the project, including the file and directory structure.
    """
    dir_list = ['/archive', '/db', '/jobpool/short', '/jobpool/priority', '/jobpool/long', '/lost+found',
                '/screeninglib/geometrylib', '/screeninglib/structurelib']
    cwd = os.getcwd()  # just in case we need this
    user = os.getlogin()

    for dir in dir_list:
        os.makedirs(project_name + dir, 0755)

    # Copy template files into the project directory
    current_path = os.path.realpath(__file__).rsplit('/', 2)[0]
    job_templates = current_path + '/job_templates'
    tmp_str = 'cp -r ' + job_templates + ' ' + cwd + '/' + project_name
    os.system(tmp_str)
    building_blocks = current_path + '/libgen/building_blocks.dat'
    gener_rules = current_path + '/libgen/generation_rules.dat'
    tmp_str = 'cp ' + building_blocks + ' ' + gener_rules + ' ' + cwd + '/' + project_name + '/screeninglib'
    os.system(tmp_str)
    for root, directories, filenames in os.walk(cwd + '/' + project_name + '/job_templates'):
        for filename in fnmatch.filter(filenames, '*.sh'):
            current_name = os.path.join(root, filename)
            tmp = project_name + filename  # prepend project name to slurm scripts
            new_name = os.path.join(root, tmp)
            tmp = 'mv ' + current_name + ' ' + new_name
            os.system(tmp)

    with open(project_name + '/' + project_name + '.config', 'w') as config:
        config.write('project_name = ' + project_name + '\n')
        config.write('user_name = ' + user + '\n')

    with open(project_name + '/queue_list.dat', 'w') as queue_list:
        lines = ['ub-hpc general-compute 0 long\n', 'ub-hpc debug 3 short\n', 'ub-hpc gpu 0 long\n', 'ub-hpc largemem 0 long\n',
                 'chemistry beta 0 long\n']
        queue_list.writelines(lines)
        queue_list.close()
