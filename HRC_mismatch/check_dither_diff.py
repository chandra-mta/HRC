#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#########################################################################################
#                                                                                       #
#       check_dither_diff.py: compare backstop dither commands and hrc hk record        #
#                             to find mismatch                                          #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Aug 14, 2019                                                   #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import time
import Chandra.Time
import numpy
import astropy.io.fits  as pyfits
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project2/Script3.6/house_keeping/dir_list'
with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
#
#--- append path to a private folders
#
sys.path.append(mta_dir)
sys.path.append(bin_dir)

import mta_common_functions     as mcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#-------------------------------------------------------------------------------------
#-- check_dither_diff: check wether dither command mismatches                       --
#-------------------------------------------------------------------------------------

def check_dither_diff(start, stop):
    """
    check wether dither command mismatches
    input:  start   --- start time
            stop    --- stop time
                        if not given, the period of the las two days is set
    output: warning email, if there are mismatch
    """
#
#--- remove a warning email file, if it still exists
#
    ifile = exc_dir + 'dither_cont'
    mcf.rm_files(ifile)
#
#--- set data collection period
#
    if start == "":
        ltime = time.strftime('%Y:%j:00:00:00', time.gmtime())    
        stop  = Chandra.Time.DateTime(ltime).secs
        start = stop - 2.0 * 86400.0
#
#--- read backstop data
#
    ifile   = data_dir + 'hrc_backstop_extracted2007'
    bs_data = mcf.read_data_file(ifile)
#
#--- find dither disable command time
#
    dis_list = []
    for ent in bs_data:
        atemp = re.split('\s+', ent)
        stime = float(atemp[1])
        if stime < start:
            continue
        elif stime > stop:
            break

        if atemp[3].upper() == 'AODSDITH':
            dis_list.append(stime)
#
#--- remove duplicates and sort
#
    dis_list = sorted(list(set(dis_list)))

    if len(dis_list) == 0:
        exit(1)
#
#--- extract dither informaiton from pcad data; widen the data collection period
#
    wstart = start - 7200 
    wstop  = stop  + 3600
    [di_time, di_data] = extract_dither_info_from_pacd(wstart, wstop)
    comp_list = []

    for k in range(0, len(di_data)):
        stime = di_time[k]
        cdata = di_data[k]
        if stime < start:
            continue
        elif stime > stop:
            break
        if cdata == 'DISA':
            comp_list.append(stime)
#
#---- remove duplicates and sort
#
    comp_list = sorted(list(set(comp_list)))

    mismatch = compare_data(dis_list, comp_list)
#
#--- if there are mismatches, create email
#
    if len(mismatch) > 0:
        line = 'HRC Dither Disable Command Mis-match: \n'
        for ent in mismatch:
            line = line + Chandra.Time.DateTime(ent).date + ' (' + str(ent) + ')\n'

        ifile = exc_dir + 'dither_cont'
        with open(ifile, 'w') as fo:
            fo.write(line)

#-------------------------------------------------------------------------------------
#-- compare_data: create lists of matched cases and mismatched cases                --
#-------------------------------------------------------------------------------------

def compare_data(data1, data2):
    """
    create lists of matched cases and mismatched cases
    input:  data1   --- data set 1
            data2   --- data set 2; data1 is compared to this
    output: save    --- a list of mismatched cases
    """
    save = []
    for ent in data1:
        range1 = ent - 180
        range2 = ent + 180
        chk = 0
        for comp in data2:
            if comp > range1 and comp < range2:
                chk = 1
                break

        if chk == 0:
            save.append(ent)

    return  save

#-------------------------------------------------------------------------------------
#-- extract_dither_info_from_pacd: extract dither information from archieved pacad data
#-------------------------------------------------------------------------------------

def extract_dither_info_from_pacd(start, stop):
    """
    extract dither information from archieved pacad data
    input:  start   --- start time
            stop    --- stop time
    output: t_list  --- a list of time
            d_list  --- a list of data (DISA/ENAB)
    """
#
#--- set arc5gl command
#
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector=pcad\n'
    line = line + 'subdetector=eng\n'
    line = line + 'level=0\n'
    line = line + 'filetype=pcad8eng\n'
    line = line + 'tstart=' + str(start) + '\n'
    line = line + 'tstop='  + str(stop)  + '\n'
    line = line + 'go\n'
#
#--- run arc5gl
#
    f_list = mcf.run_arc5gl_process(line)
    if len(f_list) < 1:
        return []
#
#--- extract the data
#
    t_list = []
    d_list = []
    for ent in f_list:
        hout   = pyfits.open(ent)
        data   = hout[1].data
        t_list = t_list + list(data['time'])
        d_list = d_list + list(data['AODITHEN'])
        hout.close()

        mcf.rm_files(ent)
    t_list = [float(ent) for ent in t_list]

    return [t_list, d_list]

#-------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) > 2:
        start = sys.argv[1].strip()
        stop  = sys.argv[2].strip()
        mc  = re.search(':', start)
        if mc is not None:
            start = Chandra.Time.DateTime(start).secs
            stop  = Chandra.Time.DateTime(stop).secs
    else:
        start = ''
        stop  = ''

    check_dither_diff(start, stop)
