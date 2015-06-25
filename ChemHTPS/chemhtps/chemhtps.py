#!/usr/bin/env python

PROGRAM_NAME = "ChemHTPS"
PROGRAM_VERSION = "v0.0.1"
REVISION_DATE = "2015-06-24"
AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
CONTRIBUTORS = """   Mohammad Atif Faiz Afzal (library generation)"""
DESCRIPTION = "ChemHTPS is a virtual high-throughput screening program suite for the chemical and materials sciences."

# Version history timeline (move to CHANGES periodically):
# v0.0.1 (2015-06-24): complete refactoring of original ChemHTPS code in new package format


###################################################################################################
# TASKS OF THIS MODULE:
# -main function
###################################################################################################


###################################################################################################
# Desired workflow (corresponding to modules, that can be called from this main module):
# 1) set up file/dir structure for a new project (setup_project)
# 2) create a library of structures (generate_structurelib)
# 3) read in the library of structure into database (populate_db)
# 4) create geometries (generate_geometries)
# 5) connect geometries to the library (populate_db)
# 6) create jobs, put jobs into pool (generate_jobs)
# 7) prioritize pool
# 8) run jobs
# 9) parse jobs
# 10) ...
###################################################################################################


###################################################################################################
#TODO:
# -restructure more general functions into modules
# -put in a printlevel restriction for each print statement 
###################################################################################################


import sys
import os
import time
# TODO: this should at some point replaced with argparser
from optparse import OptionParser

from misc import (banner,
                  format_invoked_opts,
                  tot_exec_time_str,
                  intermed_exec_timing,
                  intermed_process_timing,
                  std_datetime_str,
                  chk_rmfile,
                  chk_mkdir)

from project_setup import setup_project
from library_generator import (generate_structurelib,
                               generate_geometries)
from db_feeder import populate_db
from job_generator import (generate_jobs,
                           prioritize_pool)

###################################################################################################

def main(opts,commline_list):
    """(main):
        Driver of ChemHTPS.
    """
    time_start = time.time()
    logfile = open(opts.logfile,'a',0)
    error_file = open(opts.error_file,'a',0)

    banner_list = banner(PROGRAM_NAME, PROGRAM_VERSION, REVISION_DATE, AUTHORS, CONTRIBUTORS, DESCRIPTION)
    for line in banner_list:
        print line
        logfile.write(line + '\n')

# TODO: implement read function for local config file which populates additional opts

    fopts_list = format_invoked_opts(opts,commline_list)
    for line in fopts_list:
        print line
        logfile.write(line + '\n')

    tmp_str = "------------------------------------------------------------------------------ "
    print tmp_str
    logfile.write(tmp_str + '\n')


    if opts.setup_project:
# TODO: write test that opts.project_name exists
        setup_project(opts.project_name)

    if opts.generatelib:
        generate_structurelib()
        populate_db("moleculegraph")
        generate_geometries()
        populate_db("moleculegeom")

# TODO: we may want to put a db-based bookkeeping step in here                

    if opts.generatejobs:
        generate_jobs()

    if opts.prioritizepool:
        prioritize_pool()

# TODO: add parser function

    tmp_str = "------------------------------------------------------------------------------ "
    print tmp_str
    logfile.write(tmp_str + '\n')

    tmp_str = tot_exec_time_str(time_start) + "\n" + std_datetime_str()
    print tmp_str  + '\n\n\n'
    logfile.write(tmp_str + '\n\n\n\n')
    logfile.close()    
    error_file.close()
        
    # check whether error_file contains content
    chk_rmfile(opts.error_file)
    
    return 0    #successful termination of program
    
##################################################################################################

if __name__=="__main__":
    usage_str = "usage: %prog [options] arg"
    version_str = "%prog " + PROGRAM_VERSION
# TODO: replace with argparser; make sure that this doesn't break things in the code that rely on the opts info (e.g., the banner)
    parser = OptionParser(usage=usage_str, version=version_str)    

# TODO: we should implement a way of automatically recognizing that chemhtps is launched from within a project folder (e.g., by hte presence of a config file), so that all this doesn't have to be specified every time
# TODO: see additional comments about use of a chemhtps.config file above
    parser.add_option('--project_name', 
                      dest='project_name', 
                      type='string', 
                      help='name of the current project [default: %default]')

    parser.add_option('--setup_project', 
                      dest='setup_project', 
                      action='store_true', 
                      default=False,
                      help='sets up the infrastructure for a new project [default: %default]')

    parser.add_option('--generatelib', 
                      dest='generatelib', 
                      action='store_true', 
                      default=False,
                      help='generates a new screening library [default: %default]')

    parser.add_option('--generatejobs', 
                      dest='generatejobs', 
                      action='store_true', 
                      default=False,
                      help='generates computational chemistry jobs [default: %default]')

    parser.add_option('--prioritizepool', 
                      dest='prioritizepool', 
                      action='store_true', 
                      default=False,
                      help='prioritizes the jobs pool [default: %default]')


    # specify log files 
    parser.add_option('--logfile', 
                      dest='logfile', 
                      type='string', 
                      default='ChemHTPS.log',  
                      help='specifies the name of the log-file [default: %default]')

    parser.add_option('--errorfile', 
                      dest='error_file', 
                      type='string', 
                      default='ChemHTPS.err',  
                      help='specifies the name of the error-file [default: %default]')

    parser.add_option('--print_level', 
                      dest='print_level', 
                      type='int', 
                      default=2, 
                      help='specifies the print level for on screen and the logfile [default: %default]')


    opts, args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 2:
        sys.exit("You tried to run ChemHTPS without options.")
    main(opts,sys.argv)   #numbering of sys.argv is only meaningful if it is launched as main
    
else:
    sys.exit("Sorry, must run as driver...")
    

if __name__ == '__main__':
    pass