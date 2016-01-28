#!/usr/bin/env python

_MODULE_NAME = "library_generator"
_MODULE_VERSION = "v0.0.1"
_REVISION_DATE = "2015-06-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and Mohammad Atif Faiz Afzal (m27@buffalo.edu)"
_DESCRIPTION = "This module generates the high-throughput screening libraries."

# Version history timeline:
# v0.0.1 (2015-06-24): basic implementation

###################################################################################################
# TASKS OF THIS MODULE:
# -generate screening libraries
###################################################################################################

###################################################################################################
# TODO:
# 
###################################################################################################
import sys
import os


def generate_structurelib():
    """(generate_structurelib):
        This function generates a screening library of structures (SMILES).
    """


def generate_geometries():
    """(generate_geometries):
        This function generates the guess geometries for a screening library.
    """
    logfile = open('lib_gen.log', 'a', 0)
    error_file = open('lib_gen.err', 'a', 0)

    defaults = {'input_file': 'screeninglib/building_blocks.dat', 'rule_file': 'screeninglib/generation_rules.dat',
                'molecule_type': 'smiles', 'combination_type': 'link', 'generation_levels': '1', 'output_type': 'smi'}

    inp = raw_input('Please enter the name of the building block file: ')
    if inp != '':
        defaults['input_file'] = 'screeninglib/' + inp
    rule = raw_input('Please enter the name of the rule file: ')
    if rule != '':
        defaults['rule_file'] = 'screeninglib/' + rule
    moltyp = raw_input('Please enter the input molecule type eg. \'smiles\': ')
    if moltyp != '':
        defaults['molecule_type'] = moltyp
    combi = raw_input('Please enter link or fusion for combination type: ')
    if combi != '':
        defaults['combination_type'] = combi
    gens = raw_input('Please enter the number of generations: ')
    if gens != '':
        defaults['generation_levels'] = gens
    outtyp = raw_input('Please enter the output format eg. \'smi\' or \'xyz\': ')
    if outtyp != '':
        defaults['output_type'] = outtyp

    options = ''
    for key, value in defaults.iteritems():
        options += '--' + key + ' ' + value + ' '
    tmp_str = 'python ~/chemhtps/ChemHTPS/libgen/libgen_play.py ' + options + '--ChemHTPS'
    os.system(tmp_str)

# TODO: write this function creating all the stuff in the dummy_project, including the config file containing the project name
