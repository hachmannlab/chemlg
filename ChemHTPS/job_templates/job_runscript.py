#!/usr/bin/env python

_MODULE_NAME = "job_runscript"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2016-02-25"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module runs the jobs on the cluster."

# Version history timeline:
# v0.0.1 (2016-02-26): basic implementation
# v0.1.0 (2016-02-26): alpha version

###################################################################################################
# TASKS OF THIS MODULE:
# -Run the job after being sent to the cluster
###################################################################################################

###################################################################################################
# TODO:
#
###################################################################################################

import os
import sys
import fnmatch
import signal
import argparse
import subprocess

###################################################################################################

def main(args, commline_list):
    """
    .. function:: main(args, commline_list)
        Function to run execute jobs on the cluster

        :param object args: arguments from the parser
        :param list commline_list: the entered commandline to start the program
    """
    scratch = args.scratch_dir
    submit = args.submit_dir

    def handler(signum, frame):
        """
        .. function:: handler(signum, frame)
            Function to move files back when termination signal is sent to a job

            :param signum:
            :param frame:
        """
        os.system('mv ' + scratch + '/* ' + submit + '/.')

    signal.signal(signal.SIGTERM, handler)

    # move files from submit directory to the local scratch of the node
    tmp_str = 'mv ' + submit + '/!(slurm_orca.out) ' + scratch
    subprocess.Popen("bash -O extglob -c '" + tmp_str + "'", shell=True).communicate()
    os.chdir(scratch) # change directory to the local scratch

    # find the name of the input file
    input_file=''
    for inp in fnmatch.filter(os.listdir(scratch), '*.inp'):
        input_file = inp
    job_name = input_file.rsplit('.')[0]

    # start the job
    with open(job_name + '.out', 'w', 0) as output:
        orca = subprocess.Popen(['orca',job_name+'.inp','>',job_name+'.out','&'], stdout=output)

    # wait for program to end
    orca.communicate()
    os.system('mv ' + scratch + '/* ' + submit + '/.')

if __name__ == "__main__":
    usage_str = "usage: %(prog)s [options] arg"
    version_str = "%(prog)s " + _MODULE_VERSION
    parser = argparse.ArgumentParser(usage=usage_str)

    parser.add_argument('--version',
                        action='version',
                        version=version_str)

    parser.add_argument('--scratch_dir',
                        dest='scratch_dir',
                        default=None,
                        help='path of the scratch directory [default: %(default)s]')

    parser.add_argument('--submit_dir',
                        dest='submit_dir',
                        default=None,
                        help='path of the submit directory [default: %(default)s]')

    args = parser.parse_args(sys.argv[1:])

    main(args, sys.argv)   #numbering of sys.argv is only meaningful if it is launched as main

else:
    sys.exit("Sorry, must run as driver...")
