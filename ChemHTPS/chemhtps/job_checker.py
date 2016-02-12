#!/usr/bin/env python

_MODULE_NAME = "job_checker"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2016-02-07"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module checks for completed jobs and attempts error handling on crashed jobs."

# Version history timeline:
# v0.0.1 (2016-02-07): basic implementation
# v0.0.3 (2016-02-10): basic functionality only finds last lines for now
# v0.1.0 (2016-02-11): alpha version with full functionality for several cases

####################################################################################################
# TASKS OF THIS MODULE:
# -check for completed or crashed jobs
####################################################################################################

####################################################################################################
# TODO:
# Could probably come up with something a little neater looking maybe take some pieces out as functions
# Maybe make a class for a job that could have things like that jobs important file paths and last lines of output files
####################################################################################################

import os
import sys
import fnmatch
import datetime


def check_jobs(scratch, archive, lost):
    """(check_jobs):
        This function checks for completed or crashed job and processes them.
        :param scratch: The path of the scratch directory where jobs are run
        :param archive: The path of the archive to place finished jobs
        :param lost: The path to the lost and found folder for user attention
    """

    cwd = os.getcwd()
    logfile = open('job_checker.log', 'a', 0)
    error_file = open('job_checker.err', 'a', 0)

    if not os.path.isdir(scratch):
        error_file.write("There is no scratch folder with jobs\n")
        sys.exit('Missing scratch folder')

    for root, directories, filenames in os.walk(scratch):
        if not filenames:
            continue
        job_id = root.rsplit('/')[-1]
        job_path = root
        slurm_last = ""
        orca_last = ""
        # Find last line of both orca and slurm output
        for filename in fnmatch.filter(filenames, '*.out'):
            if filename == job_id + '.out':
                # tmp = 'ORCA output ' + os.path.join(root, filename)
                orca_out = os.path.join(root, filename)
                with open(orca_out, 'rb') as f:
                    f.seek(-2, 2)
                    while f.read(1) != b"\n":
                        f.seek(-2, 1)
                    orca_last = f.readline().split('\n')[0]
            elif 'slurm_orca' in filename:
                # tmp = 'slurm output ' + os.path.join(root, filename)
                slurm_out = os.path.join(root, filename)
                with open(slurm_out, 'rb') as f:
                    f.seek(-2, 2)
                    while f.read(1) != b"\n":
                        f.seek(-2, 1)
                    slurm_last = f.readline().split('\n')[0]
        # This is the case where the job has completed succesfully
        if slurm_last == "All Done!" and "TOTAL RUN TIME:" in orca_last:
            tmp = "tar -cjf " + job_path + ".tbz " + job_path
            os.system(tmp)
            tmp = "mv " + job_path + ".tbz " + archive
            os.system(tmp)
            tmp = "rm -rf " + job_path
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write('Job ' + job_id + ' has finished and been moved to the archive: ' + str(now) + '\n')
        # This is the case where one of the coordinates is way off
        # This also probably pops up for a number of different errors
        elif slurm_last == "All Done!" and orca_last == "ABORTING THE RUN":
            tmp = "tar -cjf " + job_path + ".coordoff.tbz " + job_path
            os.system(tmp)
            tmp = "mv " + job_path + ".coordofftbz " + lost
            os.system(tmp)
            tmp = "rm -rf " + job_path
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write(
                'Job ' + job_id + ' has not finished due to the geomery being very off moved to lost+found: ' + str(
                    now) + '\n')
        # This is the case where a coordinate is missing from the geometry file
        elif slurm_last == "All Done!" and orca_last == "No atoms to convert in Cartesian2Internal":
            tmp = "tar -cjf " + job_path + ".missingcoord.tbz " + job_path
            os.system(tmp)
            tmp = "mv " + job_path + ".missingcoord.tbz " + lost
            os.system(tmp)
            tmp = "rm -rf " + job_path
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write(
                'Job ' + job_id + ' has not finished due to a missing coordinate in the geometry: ' + str(now) + '\n')
        # This is the case where we have run out of memory
        # The orca_last here might be ORCA finished by error termination in ORCA_GTOInt but this could depend
        # on when the memory ran out
        elif slurm_last == "slurmstepd: Exceeded step memory limit at some point.":
            for root2, directories2, filenames2 in os.walk(job_path):
                for filename in fnmatch.filter(filenames2, '*.sh'):
                    slurm_script = os.path.join(root, filename)
            with open(slurm_script, 'r') as slurm:
                lines = slurm.readlines()
            with open(slurm_script, 'w') as nslurm:
                for i, line in enumerate(lines):
                    if '--mem' in line:
                        value = int(line.rsplit('=')[-1])
                        new_value = str(3 * value)
                        lines[i] = '#SBATCH --mem=' + new_value + '\n'
                nslurm.writelines(lines)
            tmp = 'sbatch ' + slurm_script
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write(
                'Job ' + job_id + " crashed due to a memory limit issue, and has been restarted: " + str(now) + '\n')
        # This is the case where the job ran out of time
        elif "DUE TO TIME LIMIT" in slurm_last:
            for root2, directories2, filenames2 in os.walk(job_path):
                for filename in fnmatch.filter(filenames2, '*.inp'):
                    input_file = os.path.join(root, filename)
            for root2, directories2, filenames2 in os.walk(job_path):
                for filename in fnmatch.filter(filenames2, '*.sh'):
                    slurm_script = os.path.join(root, filename)
            for root2, directories2, filenames2 in os.walk(job_path):
                for filename in fnmatch.filter(filenames2, '*.gbw'):
                    gbw = os.path.join(root, filename)
                    tmp = "mv " + gbw + " " + os.path.join(root, "old.gbw")
                    os.system(tmp)
            with open(input_file, 'r') as inp:
                lines = inp.readlines()
                print lines
            with open(input_file, 'w') as ninp:
                lines.insert(1, "!MORead\n")
                lines.insert(2, '%moinp "old.gbw"\n')
                ninp.writelines(lines)
            tmp = 'sbatch ' + slurm_script
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write('Job ' + job_id + ' ran out of time and has been restarted: ' + str(now) + '\n')
        else:
            tmp = "tar -cjf " + job_path + ".bad.tbz " + job_path
            os.system(tmp)
            tmp = "mv " + job_path + ".bad.tbz " + lost
            os.system(tmp)
            tmp = "rm -rf " + job_path
            os.system(tmp)
            now = datetime.datetime.now()
            logfile.write('Job ' + job_id + ' has not finished for unknown reason: ' + str(now) + '\n')
            error_file.write(
                'Job ' + job_id + ' has not finished due to a previously unknown issue: ' + str(now) + '\n')

    return 0
