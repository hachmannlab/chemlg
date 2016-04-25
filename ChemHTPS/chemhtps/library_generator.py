#!/usr/bin/env python

_MODULE_NAME = "library_generator"
_MODULE_VERSION = "v0.0.1"
_REVISION_DATE = "2015-06-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and Mohammad Atif Faiz Afzal (m27@buffalo.edu)"
_DESCRIPTION = "This module generates the high-throughput screening libraries."

# Version history timeline:
# v0.0.1 (2015-06-24): basic implementation
# v0.1.0 (2016-02-24): fixed menu to use curses

###################################################################################################
# TASKS OF THIS MODULE:
# -generate screening libraries
###################################################################################################

###################################################################################################
# TODO:
# 
###################################################################################################

import os
import curses
import fnmatch
from template_generator import runmenu, showresult
from misc import menu_input


def generate_structurelib():
    """
    .. function:: generate_structurelib()
        This function generates a screening library of structures (SMILES).
    """


def generate_geometries():
    """
    .. function:: generate_geommetries()
        This function generates the guess geometries for a screening library.
    """
    logfile = open('lib_gen.log', 'a', 0)
    error_file = open('lib_gen.err', 'a', 0)


    cwd = os.getcwd()

    menu = curses.initscr()
    curses.noecho()
    curses.cbreak()
    curses.curs_set(0)
    curses.start_color()
    menu.keypad(1)

    # Menu for the library generator options
    inp = ['Building Block File', 'Please enter the name of the building block file: ', 'Manual input', 'Exit Program']
    rule = ['Generator Rules File', 'Please enter the name of the rule file: ', 'Manual input', 'Previous Menu']
    moltyp = ['Molecule Type', 'Please enter the input molecule type: ', 'smiles', 'inchi', 'Previous Menu']
    combi = ['Combination Type', 'Please enter the combination type for the library: ', 'link', 'fusion',
             'Previous Menu']
    gens = ['Generation Level', 'Please enter the number of generations: ', 'Generation: ', 'Previous Menu']
    outtyp = ['Output Type', 'Please enter the output format', 'smi', 'xyz', 'Previous Menu']
    lib_name = ['Library Name', 'Press Enter: ', 'Press here to enter a name', 'Previous Menu']

    # dictionary of menus for navigation purposes
    menus = {0: inp, 1: rule, 2: moltyp, 3: combi, 4: gens, 5: outtyp, 6: lib_name}
    menu_names = {0: 'input_file', 1: 'rule_file', 2: 'molecule_type', 3: 'combination_type',
                  4: 'generation_levels', 5: 'output_type', 6: 'lib_name'}
    options = {0: '', 1: '', 2: '', 3: '', 4: '', 5: '', 6: ''}
    # Build out the inp and rule menu lists
    for root, directories, filenames in os.walk(cwd + '/screeninglib'):
        for filename in fnmatch.filter(filenames, '*building*'):
            menus[0].insert(2,os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, '*rules*'):
            menus[1].insert(2,os.path.join(root, filename))
    i = 0
    while i <= len(menus):
        if i ==len(menus):
            pos = showresult(menu, menu_names, options)
            if pos == 1:
                i = 0
            else:
                break
        pos = runmenu(menu, menus[i])
        options[i] = menus[i][pos]
        if i == 0 and pos == menus[0].index('Exit Program'):
            break
        elif i != 0 and pos == menus[i].index('Previous Menu'):
            i += -1
        elif i == 0 and pos == menus[0].index('Manual input'):
            options[0] = menu_input(menu, pos, len(menus[0][pos]))
            i += 1
        elif i == 1 and pos == menus[1].index('Manual input'):
            options[1] = menu_input(menu, pos, len(menus[1][pos]))
            i += 1
        elif i == 4 and pos == menus[4].index('Generation: '):
            options[4] = menu_input(menu, pos, len(menus[4][pos]))
            i += 1
        elif i == 6 and pos == menus[6].index('Press here to enter a name'):
            options[6] = menu_input(menu, pos, len(menus[6][pos])) + '_'
            i += 1
        else:
            i += 1

    curses.nocbreak()
    curses.echo()
    curses.curs_set(1)
    menu.keypad(0)
    curses.endwin()


    ops = ''
    for i in xrange(len(menus)):
        ops += '--' + menu_names[i] + ' ' + options[i] + ' '
    libgen_path = os.path.realpath(__file__).rsplit('/', 2)[0]
    libgen_path += '/libgen/libgen_play.py '
    tmp_str = 'python ' + libgen_path + ops + '--ChemHTPS'
    os.system(tmp_str)

# TODO: write this function creating all the stuff in the dummy_project, including the config file containing the project name
