#!/usr/bin/env python

PROGRAM_NAME = "ChemHTPS"
PROGRAM_VERSION = "v0.0.1"
REVISION_DATE = "2015-06-24"
AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
CONTRIBUTORS = """Mohammad Atif Faiz Afzal (library generation)"""
DESCRIPTION = "ChemHTPS is a virtual high-throughput screening program suite for the chemical and materials sciences."

# Version history timeline (move to CHANGES periodically):
# v0.0.1 (2015-06-24): complete refactoring of original ChemHTPS code in new package format


###################################################################################################
# TASKS OF THIS MODULE:
# -main function
###################################################################################################


###################################################################################################
# Desired workflow (corresponding to modules, that can be called from this main module):
# 1) set up file/dir structure for a new project
# 2) create a library of structures
# 3) read in the library of structure into database
# 4) create geometries
# 5) connect geometries to the library
# 6) create jobs
# 7) put jobs into pool
# 8) run jobs
# 9) parse jobs
# 10) ...
###################################################################################################


###################################################################################################
#TODO:
# -restructure more general functions into modules
###################################################################################################

import sys
import os
import time
# TODO: this should at some point replaced with argparser
from optparse import OptionParser


###################################################################################################

def main(opts,commline_list):
    """(main):
        Driver of ChemHTPS.
    """
    time_start = time.time()

# TODO: add banner
# TODO: add parser function
    
    return 0    #successful termination of program
    
##################################################################################################

if __name__=="__main__":
    usage_str = "usage: %prog [options] arg"
    version_str = "%prog " + PROGRAM_VERSION
# TODO: replace with argparser
    parser = OptionParser(usage=usage_str, version=version_str)    

    # it is better to sort options by relevance instead of a rigid structure
    parser.add_option('--job', 
                      dest='input_file', 
                      type='string', 
                      default='input.dat', 
                      help='input/job file [default: %default]')


    opts, args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 2:
        sys.exit("You tried to run ChemHTPS without options.")
    main(opts,sys.argv)   #numbering of sys.argv is only meaningful if it is launched as main
    
else:
    sys.exit("Sorry, must run as driver...")
    

if __name__ == '__main__':
    pass