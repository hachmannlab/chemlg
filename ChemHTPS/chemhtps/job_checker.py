#1/usr/bin/env python

_MODULE_NAME = "job_checker"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2016-03-21"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module checks for completed jobs and attempts error handling on crashed jobs."

# Version history timeline:
# v0.0.1 (2016-02-07): basic implementation
# v0.0.3 (2016-02-10): basic functionality only finds last lines for now
# v0.1.0 (2016-02-11): alpha version with full functionality for several cases
# v0.2.0 (2016-03-21): replaced with a class based system

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
import subprocess
import fnmatch
import datetime
from misc import chk_mkdir


class Job(object):
    """
    .. class:: Job(name, cluster, sbatch)
        A class to handle all aspects of a job unit

        :param str name: The name of the job unit
        :param str cluster: Which cluster the job was submitted to
        :param str sbatch: The slurm submission string
        :param str path: The path to the location the job is run from
        :param str slurm_id: The slurm job id of the job
        :param str slurm_last_line: The last line of the slurm output file
        :param str out_last_line: The last line of the quantum software output file
        :param bool is_running: A boolean that is true if the job is running and False if its not
        :param str rm_path: a path for removal, if populated can be deleted because the job has been tarred
    """

    def __init__(self, name, cluster, sbatch):
        """
        .. method:: __init__(self, name, cluster, sbatch)
            Initialize a Job object.
        """
        self.name = name
        self.cluster = cluster
        self.path = os.getcwd()
        self.slurm_id = subprocess.check_output(sbatch, shell=True).split()[3]
        self.is_running = True
        self.rm_path = ''
        self.slurm_last_line = ''
        self.out_last_line = ''


    def check_status(self, user_name):
        """
        .. method:: check_status(self, user_name)
            Checks if the slurm job has finished or not and updates the is_running variable

            :return self.is_running: The state of the job
            :rtype: bool
        """
        tmp = "squeue -M " + self.cluster + " -u " + user_name + " | grep " + self.slurm_id + " | wc -l"
        number = int(subprocess.check_output(tmp, shell=True))
        if number == 0:
            self.is_running = False
        else:
            pass
        return self.is_running

    def time_limit_restart(self):
        """
        .. method:: time_limit_restart(self)
            Makes the necessary changes to the input file and gbw file to restart a job that ran out of time
        """
        input_file = self.path + '/' + self.name + '.inp'
        gbw = self.path + '/' + self.name + '.gbw'
        tmp = "mv " + gbw + " " + self.path + '/' + 'old.gbw'
        os.system(tmp)
        with open(input_file, 'r') as inp:
            lines = inp.readlines()
        lines.insert(1, "!MORead\n")
        lines.insert(2, '%moinp "old.gbw"\n')
        with open(input_file, 'w') as ninp:
            ninp.writelines(lines)

    def nth_line(self, n, file_name):
        """
        .. method:: nth_line(self, n, file_name)
            Returns the nth from the end line of the specified file

            :param str file_name: The name of the file
            :param int n: Number of lines from the end
            :return line: The nth line
            :rtype: str
        """
        with open(self.path + '/' + file_name) as out:
            out.seek(-1, 2)
            i = 0
            while i <= n:
                if out.read(1) == '\n':
                    i += 1
                out.seek(-2,1)
            out.seek(2,1)
            tmp = out.readlines()
            #tmp = filter(lambda a: a != '\n', tmp)
            line = tmp[0].strip('\n')
            return line

    def slurm_last(self):
        """
        .. method:: slurm_last(self)
            Get the last line of the slurm output
        """
        self.slurm_last_line = self.nth_line(1, 'slurm_orca.out')
            
    def out_last(self):
        """
        .. method:: orca_last(self)
            Get the last line of the quantum software output
        """
        self.out_last_line = self.nth_line(1, self.name + ".out")

    def tar_job_unit(self, tbz='.tbz'):
        """
        .. method:: tar_job_unit(self, tbz='.tbz')
            Tars the jobunit preparing it for transport

            :param str tbz: The tar file extension
        """
        cwd = os.getcwd()
        os.chdir(self.path.rsplit('/', 1)[0])
        tmp = "tar -cjf " + self.name + tbz + " " + self.name
        os.system(tmp)
        self.rm_path = self.path
        self.path = self.path + tbz
        os.chdir(cwd)

    def move_job(self, dest_path):
        """
        .. method:: move_job(self, dest_path)
            moves the job unit to the specified location

            :param str dest_path: The destination path where the jobunit is going
        """
        tmp = "mv " + self.path + " " + dest_path
        print tmp
        os.system(tmp)

    def rm_job(self):
        """
        .. method:: rm_job(self)
            Deletes the job folder after the job has been tarred
        """
        tmp = "rm -r " + self.rm_path
        os.system(tmp)


def check_jobs(user_name, scratch, archive, lost, job_list):
    """
    .. function:: check_jovs(user_name, scratch, archive, lost, job_list)
        This function checks for completed or crashed jobs and processes them.

        :param str scratch: The path of the scratch directory where jobs are run
        :param str archive: The path of the archive to place finished jobs
        :param str lost: The path of the lost and found folder for user attention
        :param str job_list: A list of current job units to check
    """

    cwd = os.getcwd()
    logfile = open('job_checker.log', 'a', 0)
    error_file = open('job_checker.err', 'a', 0)

    if not os.path.isdir(scratch):
        error_file.write("There is no scratch folder with jobs\n")
        sys.exit('Missing scratch folder')

    rem_list = []
    for job in job_list:
        if not job.check_status(user_name):
            job.slurm_last()
            job.out_last()
            if job.slurm_last_line == "All Done!" and "TOTAL RUN TIME:" in job.out_last_line:
                job.tar_job_unit()
                job.move_job(archive)
                job.rm_job()
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' has finished and been moved to the archive: ' + str(now) + '\n')
            elif job.slurm_last_line == "All Done!" and job.out_last_line == "ABORTING THE RUN":
                job.tar_job_unit('.coordoff.tbz')
                job.move_job(lost)
                job.rm_job()
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' has not finished due to the geometry being very off, moved to lost+found: ' + str(now) + '\n')
            elif job.slurm_last_line == "All Done!" and job.out_last_line == "No atoms to convert in Cartesian2Internal":
                job.tar_job_unit('.missingcoord.tbz')
                job.move_job(lost)
                job.rm_job()
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' has not finished due to a missing coordinate in the geometry, moved to lost+found: ' + str(now) + '\n')
            elif job.slurm_last_line == "slurmstepd: Exceeded step memory limit at some point.":
                restarts = cwd + '/jobpool/priority/restarts'
                chk_mkdir(restarts)
                job.move_job(restarts)
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' crashed due to a memory issue, and has been restarted: ' + str(now) + '\n')
            elif "DUE TO TIME LIMIT" in job.slurm_last_line:
                job.time_limit_restart()
                restarts = cwd + '/jobpool/priority/restarts'
                chk_mkdir(restarts)
                job.move_job(restarts)
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' ran out of time and has been restarted: ' + str(now) + '\n')
            else:
                job.tar_job_unit('.bad.tbz')
                job.move_job(lost)
                job.rm_job()
                now = datetime.datetime.now()
                logfile.write('Job ' + job.name + ' has not finished for unknown reason: ' + str(now) + '\n')
                error_file.write('Job ' + job.name + ' has not finished due to a previously unknown issue: ' + str(now) + '\n')
        else:
            rem_list.append(job)
            pass
    

    return rem_list
