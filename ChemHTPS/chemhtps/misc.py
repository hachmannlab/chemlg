#!/usr/bin/env python

_MODULE_NAME = "misc"
_MODULE_VERSION = "v1.0.0"
_REVISION_DATE = "2015-06-24"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This is a module for miscellaneous, general purpose functions."

# Version history timeline:
# v1.0.0 (2015-06-24): adaptation of old lib_jcode
# v1.0.1 (2016-02-24): added a menu input function

###################################################################################################
# TASKS OF THIS MODULE:
# -provides miscellaneous, general purpose functions
###################################################################################################

###################################################################################################
# TODO:
# 
###################################################################################################

import sys
import os
import struct
import time
import datetime
import subprocess
import hashlib
import mmap
import curses
from numpy import fromstring

###################################################################################################

# TODO: check out similar codes and see if we should add anything
def banner(SCRIPT_NAME, SCRIPT_VERSION, REVISION_DATE, AUTHORS, CONTRIBUTORS, DESCRIPTION,):
    """
    .. function:: banner(SCRIPT_NAME, SCRIPT_VERSION, REVISION_DATE, AUTHORS, CONTRIBUTORS, DESCRIPTION)
        Banner for python scripts.

        :param SCRIPT_NAME: The name of the script
        :param SCRIPT_VERSION: The version number
        :param REVISION_DATE: The latest revision date
        :param AUTHORS: The main authors
        :param CONTRIBUTORS: Any contributing authors
        :param DESCRIPTION: A brief description
        :return banner_list: A nicely formatted banner
        :rtype: list
    """
    banner_list = []
    banner_list.append("============================================================================== ")
    banner_list.append(SCRIPT_NAME + " " + SCRIPT_VERSION + " (" + REVISION_DATE + ")")
    banner_list.append(AUTHORS)
    banner_list.append("============================================================================== ")
    banner_list.append(time.ctime())
    banner_list.append("")    
    banner_list.append(DESCRIPTION)
    banner_list.append("With contributions from:")
    banner_list.append(CONTRIBUTORS)
    banner_list.append("")

    return banner_list

##################################################################################################

def format_invoked_opts(opts,commline_list=[]):
    """
    .. function:: format_invoked_opts(opts, commline_list=[])
        Prints the invoked options to stdout and the logfile.

        :param object opts: An object storing the relevant options
        :param list commline_list: The entered command line arguments
        :return fopts_list: A formatted list of options
        :rtype: list
    """    
    fopts_list = []
    if len(commline_list) != 0:
        fopts_list.append("Invoked command line: ")
        fopts_list.append(" ".join(commline_list))
        fopts_list.append("")
        
    fopts_list.append("Invoked options: ")    
    for key, value in opts.__dict__.items():
        fopts_list.append("   " + key + ": " + str(value))
    fopts_list.append("")
    return fopts_list
           

##################################################################################################

def wc_dir(dir):
    """
    .. function:: wc_dir(dir)
        Returns the number of dirs in a given dir via ls -1d | wc -l. 
        Note that this becomes a rather expensive function call when dir contains many subdirs.
    """    
#TODO: take care of error for empty dirs
    tmp_str = "ls -1d " + dir + "/*/ | wc -l"
    # this is a quite new python feature and may is only available in 2.6 or so 
    # n = subprocess.getoutput(tmp_str)
    # ... and for older python versions     
    return int(subprocess.Popen(tmp_str,shell=True,stdout=subprocess.PIPE).stdout.read())

##################################################################################################

def wc_all(dir):
    """
    .. function:: wc_all(dir)
        Returns the number of files and dirs in a given dir via ls -1 | wc -l. 
        Not that this becomes a rather expensive function call when dir contains many entries.
    """    
#TODO: take care of error for empty dirs
    tmp_str = "ls -1 " + dir + " | wc -l"
    # this is a quite new python feature and may is only available in 2.6 or so 
    # n = subprocess.getoutput(tmp_str)
    # ... and for older python versions     
    return int(subprocess.Popen(tmp_str,shell=True,stdout=subprocess.PIPE).stdout.read())

##################################################################################################

def line_count(file_namestr):
    """
    .. line_count(file_namestr)
        Returns the number of lines in a file.
    """    
    if os.path.getsize(file_namestr) == 0:
        return 0
    with open(file_namestr) as file:
        for i, l in enumerate(file):
            pass
    return i + 1

##################################################################################################

def mksubdir_struct(dir,max_n_entries=10000,run_always=False):
    """
    .. function:: mksubdir_struct(dir, max_n_entries=10000, run_always=False)
        This function takes the content of a dir and makes numbered substructure dirs with each n_entries of the original dir.
        The motivation was to have a function with limits the number of entries in a directory to a certain threshold
        (e.g., 10,000 or 30,000) in order to avoid performance issues with the OS/filesystem.

         :param str dir: Path of the relevant directory
         :param int max_n_entries: Max entries in a subdir
         :param bool run_always: whether to run even if there are less than :param:'max_n_entries' in :param:'dir'
    """
    entry_list = []
    for entry in os.listdir(dir):
        entry_list.append(entry)
    entry_list.sort()
    
    n_entries = len(entry_list)
    
    if n_entries >= max_n_entries or run_always:
        subdir_counter = 0
        subdir_entry_counter = 0
        subdir_pathstr = dir + "/%05d"  %(subdir_counter)
        
        if chk_mkdir(subdir_pathstr,True) == False:
            sys.exit("Naming conflict!")
        
        for entry in entry_list:
            tmp_str = "mv " + entry + " " + subdir_pathstr + "/." 
            os.system(tmp_str)
            subdir_entry_counter +=1
            if subdir_entry_counter >= max_n_entries:
                subdir_counter += 1
                subdir_entry_counter = 0
                subdir_pathstr = dir + "/%05d"  %(subdir_counter)
                if chk_mkdir(subdir_pathstr,True) == False:
                    sys.exit("Naming conflict!")
                
##################################################################################################

def chk_mkdir(dir,warning=False):
    """
    .. function:: chk_mkdir(dir, warning=False)
        This function checks whether a directory exists and if not creates it.
    """
    if not os.path.isdir(dir):
        tmp_str = "mkdir -p " + dir
        os.system(tmp_str)
    elif warning:
        return False

##################################################################################################

def chk_rmdir(dir,check='any'):
    """
    .. function:: chk_rmdir(dir,check='any')
        This function checks whether a directory exists and removes it, if it is empty.
    """
    if os.path.isdir(dir):
        n_dirs = 0
        n_files = 0
        for i in os.listdir(dir):
            if os.path.isdir(dir + '/' + i):
                n_dirs += 1
            elif os.path.isfile(dir + '/' + i):
                n_files += 1
        if n_dirs == 0 and n_files == 0:
            tmp_str = "rm -rf " + dir
        elif n_dirs == 0 and check=='dirs':
            tmp_str = "rm -rf " + dir
        elif n_files == 0 and check=='files':
            tmp_str = "rm -rf " + dir
        else:
            tmp_str = " "
        os.system(tmp_str)

##################################################################################################

def chk_rmfile(file_namestr):
    """
    .. function:: chk_rmfile(file_namestr)
        This function checks whether a file is empty and if yes deletes it.
    """
    file = open(file_namestr,'r')
    test_str = file.read()
    file.close()
    if len(test_str) == 0:
        os.remove(file_namestr)
    
##################################################################################################

def target_dir_struct(target_dir_path, maxitems = 10000, digits=5):
    """
    .. function:: target_dir_struct(target_dir_path, maxitems = 10000, digits=5)
        This function checks whether a target dir exists and establishes/checks the subdir structure.

        :param str target_dir_path: The path of the directory
        :param int maxitems: The max number of items in a directory
        :param int digits: The number of digits in the folder names
        :return int target_subdir: The number/name of the subdir
        :return int target_subdir_n: The number of items in the subdir
        :return str target_subdir_pathstr: The path of the subdir
    """
    # check if target_dir exists and if not create it
    chk_mkdir(target_dir_path)
    # establish target_dir structure
    # 1) get all the present subdirs
    target_subdir_list = [] # fill with all present subfolders
    for i in os.listdir(target_dir_path):
        if os.path.isdir(target_dir_path + '/' + i) and i not in target_subdir_list:
            target_subdir_list.append(i)
    # 2a) if there are no subfolders present
    if len(target_subdir_list)==0:
        target_subdir = 0   # this is the highest folder
        target_subdir_n = 0 # this is the number of items in it
    # 2b) if there are subfolders present    
    else:
        target_subdir_list.sort()
        target_subdir = int(target_subdir_list[-1]) # pick the highest folder
        target_subdir_n = wc_all(target_dir_path + '/' + target_subdir_list[-1])
        if target_subdir_n >= maxitems:     # this limit is more important for folders rather than files (in this case tarballs); but we do it anyways
            target_subdir += 1
            target_subdir_n = 0

    target_subdir_pathstr = target_dir_path + '/' + '{num:{fill}{width}}'.format(num=target_subdir, fill='0', width=digits)
#    target_subdir_pathstr = target_dir_path + "/%05d"  %(target_subdir)
    chk_mkdir(target_subdir_pathstr)
    return target_subdir, target_subdir_n, target_subdir_pathstr

##################################################################################################

def mv2subdir_struct(source_dir_pathstr, target_subdir, target_subdir_n, target_subdir_pathstr, maxitems = 10000):
    """
    .. function:: mv2subdir_struct(source_dir_pathstr, target_subdir, target_subdir_n, target_subdir_pathstr, maxitems = 10000)
        This function moves a source folder into a target subdir structure and updates it.

        :param str source_dir_pathstr: The path of the dir to be moved
        :param int target_subdir: The number of the target dir
        :param int target_subdir_n: The number of items in the target dir
        :param str target_subdir_pathstr: The path of the target dir
        :param int maxitems: The max number if items in a directory
        :return int target_subdir: The number of the subdir
        :return int target_subdir_n: The number of items in target subdir
        :return str target_subdir_pathstr: The path of the target subdir
    """
    # move
    tmp_str = 'mv ' + source_dir_pathstr + ' ' + target_subdir_pathstr + '/. ' 
    os.system(tmp_str)
    target_subdir_n += 1

    # check if limit is reached
    if target_subdir_n >= maxitems:     # this limit is more important for folders rather than files (in this case tarballs); but we do it anyways
        target_subdir += 1
        target_subdir_n = 0
        
        # make new target subdir
        tmp_str = target_subdir_pathstr.split('/')[-1]
        digits = len(tmp_str)
        target_subdir_pathstr = target_subdir_pathstr[:-digits] + '{num:{fill}{width}}'.format(num=target_subdir, fill='0', width=digits)
        chk_mkdir(target_subdir_pathstr)
    return target_subdir, target_subdir_n, target_subdir_pathstr
    
##################################################################################################

def std_datetime_str(mode='datetime'):
    """
    .. function:: std_datetime_str(mode='datetime')
        This function gives out the formatted time as a standard string, i.e., YYYY-MM-DD hh:mm:ss.
    """
    if mode == 'datetime':
        return str(datetime.datetime.now())[:19]
    elif mode == 'date':
        return str(datetime.datetime.now())[:10]
    elif mode == 'time':
        return str(datetime.datetime.now())[11:19]
    elif mode == 'datetime_ms':
        return str(datetime.datetime.now())
    elif mode == 'time_ms':
        return str(datetime.datetime.now())[11:]
    else:
        sys.exit("Invalid mode!")

##################################################################################################
def tot_exec_time_str(time_start):
    """
    .. function:: tot_exec_time_str(time_start)
        This function gives out the formatted time string.
    """
    time_end = time.time()
    exec_time = time_end-time_start
    tmp_str = "Total execution time: %0.2fs (%dh %dm %0.2fs)" %(exec_time, exec_time/3600, (exec_time%3600)/60,(exec_time%3600)%60)
    return tmp_str

##################################################################################################

def intermed_exec_timing(time_start,intermed_n,total_n,n_str="n"):
    """
    .. function:: intermed_exec_timing(time_start, intermed_n, total_n, n_str='n')
        This function gives out the intermediate timing, speed, pace, projected remaining and end time.
    """
    tmp_time = time.time()
    tmp_exec_time = tmp_time-time_start
    sec_per_n = 1.0*tmp_exec_time/intermed_n
    n_per_hour = 3600.0/sec_per_n
    proj_rest_sec = sec_per_n*(total_n-intermed_n)
    proj_end_time = int(round(tmp_time + proj_rest_sec))
    tmp_str = "   Current speed: %0.2f " %(n_per_hour)
    tmp_str += n_str + "'s/hour; current pace: %0.3f " %(sec_per_n)
    tmp_str += "sec/" + n_str + "\n" 
#    tmp_str +="   Projected remaining time: %0.2fs (%dh %dm %0.2fs) " %(proj_rest_sec, proj_rest_sec/3600, (proj_rest_sec%3600)/60,(proj_rest_sec%3600)%60)
    tmp_str +="   Projected remaining time: %0.2fs (%dh %dm %0.2fs) \n" %(proj_rest_sec, proj_rest_sec/3600, (proj_rest_sec%3600)/60,(proj_rest_sec%3600)%60)
    tmp_str +="   Projected end time: " + time.ctime(proj_end_time) 
    return tmp_str

##################################################################################################

def intermed_process_timing(time_start,process_n,intermed_n,total_n,n_str="n"):
    """
    .. function:: intermed_process_timing(time_start, process_n, intermed_n, total_n, n_str='n')
        This function gives out the intermediate timing, speed, pace, projected remaining and end time of a particular process with restarted time.
    """
    tmp_time = time.time()
    tmp_exec_time = tmp_time-time_start
    if process_n == 0:
        return ''
    
    sec_per_n = 1.0*tmp_exec_time/process_n
    n_per_hour = 3600.0/sec_per_n
    proj_rest_sec = sec_per_n*(total_n-intermed_n)
    proj_end_time = int(round(tmp_time + proj_rest_sec))
    tmp_str = "   Current speed: %0.2f " %(n_per_hour)
    tmp_str += n_str + "'s/hour; current pace: %0.3f " %(sec_per_n)
    tmp_str += "sec/" + n_str + "\n" 
#    tmp_str +="   Projected remaining time: %0.2fs (%dh %dm %0.2fs) " %(proj_rest_sec, proj_rest_sec/3600, (proj_rest_sec%3600)/60,(proj_rest_sec%3600)%60)
    tmp_str +="   Projected remaining time: %0.2fs (%dh %dm %0.2fs) \n" %(proj_rest_sec, proj_rest_sec/3600, (proj_rest_sec%3600)/60,(proj_rest_sec%3600)%60)
    tmp_str +="   Projected end time: " + time.ctime(proj_end_time) 
    return tmp_str

##################################################################################################

def timeit(func):
    """
    .. function:: timeit(func)
        Annotate a function with its elapsed execution time.
    """
    def timed_func(*args, **kwargs):
        t1 = time.time()
        
        try:
            func(*args, **kwargs)
        finally:
            t2 = time.time()

        timed_func.func_time = ((t2 - t1) / 60.0, t2 - t1)

        if __debug__:
            sys.stdout.write("%s took %0.3fm %0.3fs %0.3fms\n" % (
                func.func_name,
                timed_func.func_time[0],
                timed_func.func_time[1],
            ))

    return timed_func

###################################################################################################

def dsu_sort(list, index, reverse=False):
    """
    .. function:: dsu_sort(list, index, reverse=False)
    """
# TODO: infoline
    for i, e in enumerate(list):
        list[i] = (e[index], e)
    if reverse:
        list.sort(reverse=True)
    else:
        list.sort()
    for i, e in enumerate(list):
        list[i] = e[1]
    return list

###################################################################################################

def dsu_sort2(list, index, reverse=False):
    """
    .. function:: dsu_sort2(list, index, reverse=False)
        This function sorts only based on the primary element, not on secondary elements in case of equality.
    """
    for i, e in enumerate(list):
        list[i] = e[index]
    if reverse:
        list.sort(reverse=True)
    else:
        list.sort()
    for i, e in enumerate(list):
        list[i] = e[1]
    return list

###################################################################################################

def isFloat(x):
    if x in ('NAN','NaN','Nan','nan'):
        return 0
    elif '#IND' in x:
        return 0
    try:
        float(x)
        return 1
    except:
        return 0

###################################################################################################

def md5checksum(file_path, blocksize=8192): # 4kB blocks
    """
    .. function:: md5checksum(file_path, blocksize=8192)
        Compute md5 hash of the specified file.
    """
    file = open(file_path, 'rb')
    md5sum = hashlib.md5()
    while True:
        data = file.read(blocksize)
        if not data:
            break
        md5sum.update(data)
    file.close()
    return md5sum.hexdigest()

###################################################################################################

def filelinecount(filename):
    """
    .. function:: filelinecount(filename)
        Counts the number of lines in a file.
    """
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

###################################################################################################

def revdict_lookup(dict,lookup_val):
    """
    .. function:: revdict_lookup(dict, lookup_val)
        Performs a reverse dictionary lookup. Careful: only returns first match, but there may be others.
    """    
    key = (key for key,value in dict.items() if value==lookup_val).next()
    return key


###################################################################################################

def queryset_iterator(queryset, chunksize=1000, reverse=False, id_only=False, values= False):
    """
    .. function:: queryset_iterator(queryset, chunksize=1000, reverse=False, id_only=False, values-False)
        Django incremental queryset iterator.
        Found on: http://www.poeschko.com/2012/02/memory-efficient-django-queries/
    """    
    ordering = '-' if reverse else ''
    queryset = queryset.order_by(ordering + 'pk')
    last_pk = None
    new_items = True
    while new_items:
        new_items = False
        chunk = queryset
        if last_pk is not None:
            func = 'lt' if reverse else 'gt'
            chunk = chunk.filter(**{'pk__' + func: last_pk})
        chunk = chunk[:chunksize]
        if id_only:
            chunk = chunk.values('pk')            
        row = None
        for row in chunk:
            yield row
        if row is not None:
            if id_only or values:
                last_pk = row['pk']                
            else:
                last_pk = row.pk
            new_items = True

###################################################################################################

def list_chunks(l, n):
    """
    .. function:: list_chunks(l, n)
        Yield successive n-sized chunks from l.
        Found on: http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

###################################################################################################

def menu_input(menu, row, col):
    """
    .. function:: menu_input(menu, row, col)
        Function to allow the user to enter text in a curses menu

        :param menu: The curses window
        :type menu: window object
        :param int row: The row where the cursor shows up
        :param int col: The col where the cursor shows up
        :return user_input: The input from the user
         :rtype: str
    """
    curses.echo()
    curses.curs_set(1)
    user_input = menu.getstr(row+5, col+8, 100)
    curses.noecho()
    curses.curs_set(0)
    return user_input
