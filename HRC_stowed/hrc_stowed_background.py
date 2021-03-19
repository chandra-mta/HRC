#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#############################################################################################
#                                                                                           #
#      hrc_stowed_background.py: main script to update hrc stowed background data/site      # 
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Mar 19, 2021                                                   #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import random
import time
import math
import unittest
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project1/Scripts/house_keeping/dir_list'

with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var   = atemp[1].strip()
    line  = atemp[0].strip()
    exec("%s = %s" %(var, line))

sys.path.append(mta_dir)
sys.path.append(bin_dir)
#
#--- import several functions
#
import mta_common_functions       as mcf    #---- contains other functions commonly used in MTA scripts
import hrc_stowed_common_function as hcf    #---- contains hrc stored related functions
import run_hrc_process            as rhp    #---- contains run_hrc_process related function
import run_next_in_line           as rnl    #---- contains functions related extraction of next_in_line data 
import hrc_nil_table              as hnt    #---- contains funcions to create stat tables
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail) + 'hrc'
#
#--- set  a few things
#
NULL      = 'null'
comb_dirs = ('hrc_i_90','hrc_i_115','hrc_s_90','hrc_s_125','hrc_s_90_hi','hrc_s_125_hi')

#-----------------------------------------------------------------------------------------------
#-- hrc_stowed_background: controlling script to set up directory and run all script          --
#-----------------------------------------------------------------------------------------------

def hrc_stowed_background(syear, smonth, eyear, emonth):
    """
    controlling script to set up directory and run all script
    input:  syear   --- starting year
            smonth  --- starting month
            eyear   --- ending year
            emonth  --- ending month
    output: evt0, evt1, ss0, and hk00 fits files corresponding to next_in_line condition
    """
#
#--- update hrc gain list
#
    rhp.update_gain_selection_file()
#
#--- run over the given time period
#
    for year in range(syear, eyear+1):
        for month in range(1, 13):
            if year == syear and month < smonth:
                continue
            elif year == eyear and month > emonth:
                break
#
#--- start and stop time in mm/dd/yy,hh:mm:ss format
#
            begin = hcf.conv_time_format(year, month)         
            end   = hcf.conv_time_format(year, month, next=1)
#
#--- start and stop time in seconds from 1998.1.1
#
            start = hcf.convertto1998sec(begin)
            stop  = hcf.convertto1998sec(end)
#
#--- make saving directory
#
            lmon   = hcf.find_month(month)
            outdir = data_dir  + str(year) + lmon + '/'

            if not os.path.isdir(outdir):
                cmd = 'mkdir ' + outdir
                os.system(cmd)

            for ent in comb_dirs:
                udir = outdir + ent
                if not os.path.isdir(udir):
                    cmd = 'mkdir ' +  udir
                    os.system(cmd)
#
#--- now run the main script
#
            rnl.extract_next_in_line_data(begin, end, start, stop, outdir)
#
#--- move other files to appropriated directories
#
            cmd = 'mv tscpos_positive ' + outdir
            os.system(cmd)

            cmd = 'ls * > ' + zspace
            os.system(cmd)

            data = mcf.read_data_file(zspace, remove=1)

            for ent in data:
                mc = re.search('period', ent)
                if mc is None:
                    continue 

                atemp = re.split('hrc_', ent)
                btemp = re.split('_', atemp[1])
                dname = 'hrc_' + btemp[0] + '_' + btemp[1]
                mc    = re.search('_hi', ent)
                if mc is not None:
                    dname = dname + '_hi' 
                
                cmd = 'mv ' + ent + ' ' + outdir + dname + '/'
                os.system(cmd)

            cmd = 'rm -f *fits  ./Temp_dir/* 2> /dev/null'
            os.system(cmd)

            cmd  = 'gzip -f ' + outdir + '*/*.fits  2> /dev/null' 
            os.system(cmd)
#
#---- update stat tables
#
            hnt.hrc_nil_table(year, month)
            hnt.status_bit_report(year, month)
#
#---- clean up stat files
#
    clean_up_stat_lists()

#-----------------------------------------------------------------------------------------------
#-- clean_up_stat_lists: sort and clean up the statistics related files                       --
#-----------------------------------------------------------------------------------------------

def clean_up_stat_lists():
    """
    sort and clean up the statistics related files
    input:  none, but read from <data_dir>/Stats/*
    output: clean uped data files
    """

    for part in ('*dead*', '*stat_results*', '*status*'):
#
#--- make a list of related file names
#
        cmd  = 'ls ' + data_dir + 'Stats/' + part + '  > ' + zspace
        os.system(cmd)

        slist = mcf.read_data_file(zspace, remove=1)
#
#--- update each file
#
        for lname in slist:
            data = mcf.read_data_file(lname)
#
#--- dead time file
#
            if part == '*dead*':
                update_ddate_list(data, lname)
#
#--- two others
#
            else:
                update_ldate_list(data, lname)

#-----------------------------------------------------------------------------------------------
#-- update_ddate_list: sort and clean up a file with digit time stamp                         --
#-----------------------------------------------------------------------------------------------

def update_ddate_list(data, lname):
    """
    sort and clean up a file with digit time stamp
    input:  data    --- data list with the first entry is digit time stamp
            lname   --- output file name
    output: lname   --- updated file
    """
#
#--- clean up the list
#
    out = hcf.remove_duplicate_by_column(data, 0, sep = '\s+')

    sline = ''
    for line in out:
        sline = sline + line + '\n'

    with open(lname, 'w') as fo:
        fo.write(sline)

#-----------------------------------------------------------------------------------------------
#-- update_ldate_list: sort and clean up a file with literal date time stamp                  --
#-----------------------------------------------------------------------------------------------

def update_ldate_list(data, lname):
    """
    sort and clean up a file with literal date time stamp
    input:  data    --- data list with the first entry is literal time stamp
            lname   --- output file name
    output: lname   --- updated file
    """
#
#--- make a dictionary with the time element as index
#
    save = {}
    tind = []
    for ent in data:
        atemp = re.split('\s+', ent)
        save[atemp[0]] = ent
        tind.append(atemp[0])
#
#--- remove duplicated index
#
    clist = hcf.remove_duplicate_from_list(tind)
#
#--- print out the result
#
    line = ''
    for ent in clist:
        line = line + save[ent] + '\n'

    with open(lname, 'w') as fo:
        fo.write(line)

#-----------------------------------------------------------------------------------------------
 
if __name__ == "__main__":
#
#--- if you like to specify the date period, give
#---  a year and starting yday and ending yday
#
    test = 0
    if len(sys.argv) == 2:
        if sys.argv[1] == 'test':
            test = 1
            del sys.argv[1:]
        else:
            exit(1)

    elif len(sys.argv) == 3:
        syear  = int(float(sys.argv[1]))
        smonth = int(float(sys.argv[2]))
        eyear  = syear
        emonth = smonth
    elif len(sys.argv) == 5:
        syear  = int(float(sys.argv[1]))
        smonth = int(float(sys.argv[2]))
        eyear  = int(float(sys.argv[3]))
        emonth = int(float(sys.argv[4]))
#
#--- if the date period is not specified,
#--- the period is set from the last entry date to
#--- the last day of the last month
#
    else:
        syear  = ''
        smonth = ''
        eyear  = ''
        emonth = ''

    if test == 0:
        hrc_stowed_background(syear, smonth, eyear, emonth)

        cmd = 'rm /tm/zspace*hrc 2>/dev/null'
        os.system(cmd)
    else:
        unittest.main()

