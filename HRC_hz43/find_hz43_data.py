#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################
#                                                                               #
#       find_hz43_data.py: find unprocessed hz43 data                           #
#                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                           #
#                                                                               #
#           Last Update: Jan 14, 2021                                           #
#                                                                               #
#################################################################################

import sys
import os
import string
import re
import math
import time
import random
import numpy
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

with  open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
html_top = html_top.replace('#', ':')
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

#-----------------------------------------------------------------------------------------
#-- extract_calb_hz43: create a list of calibration observations of HZ43 from the database 
#-----------------------------------------------------------------------------------------

def extract_calb_hz43():
    """
    create a list of calibration observations of HZ43 from the database
    input:  none
    output: <house_keeping>/hrc_i   --- a list of hz43 observation on hrc i
            <house_keeping>/hrc_s   --- a list of hz43 observation on hrc s
            chk                     --- if htere are new entries, return 1, otherwise 0
    """
#
#--- read the past data
#
    hrc_i_p   = make_past_obsid_list('i')
    hrc_s_p   = make_past_obsid_list('s')
#
#--- read database
#
    data = mcf.read_data_file('/data/mta4/obs_ss/sot_ocat.out')
#
#--- extract HZ43 calibration data
#
    hrc_i = []
    hrc_s = []
    for ent in data:
        mc1  = re.search('HZ43',      ent)
        mc1a = re.search('HZ 43',     ent)
        mc1b = re.search('hz43',      ent)
        mc2  = re.search('CAL',       ent)
        mc3  = re.search('archived',  ent)
        mc4  = re.search('HRC-I',     ent)
        mc5  = re.search('HRC-S',     ent)

        if (mc1 is not None) or (mc1a is not None) or (mc1b is not None):
            if mc2 is not None:
                if mc3 is not None:
                    atemp = re.split('\^', ent)
                    obsid = atemp[1].strip()
                    obsid = int(float(obsid))
                    if mc4 is not None:
                        hrc_i.append(obsid)
                    elif mc5 is not None:
                        hrc_s.append(obsid)
#
#--- check whether there are any new hz43 calibration data
#
    hrc_i_p = list(set(hrc_i_p))
    hrc_s_p = list(set(hrc_s_p))
    hrc_i   = list(set(hrc_i))
    hrc_s   = list(set(hrc_s))
    idiff   = []
    sdiff   = []
    new     = [idiff, sdiff]
    chk     = 0
    idiff   = numpy.setdiff1d(hrc_i, hrc_i_p)
    sdiff   = numpy.setdiff1d(hrc_s, hrc_s_p)

    if (len(idiff) > 0) or (len(sdiff) > 0):
        chk = 1
#
#--- if so, update the list
#
    if chk > 0:
        ifile_i = house_keeping + 'hrc_i_list'
        with open(ifile_i, 'w') as fo:
            for ent in list(hrc_i):
                fo.write(str(ent) + '\n')
    
        ifile_s = house_keeping + 'hrc_s_list'
        with  open(ifile_s, 'w') as fo:
            for ent in list(hrc_s):
                fo.write(str(ent) + '\n')

        new = [idiff, sdiff]

    chk_save = house_keeping + 'chk_save'
    with open(chk_save, 'w') as fo:
        fo.write(str(chk) + '\n')

    return new

#-----------------------------------------------------------------------------------------
#-- make_past_obsid_list: create a list of obsids of the past observations              --
#-----------------------------------------------------------------------------------------

def make_past_obsid_list(inst):
    """
    create a list of obsids of the past observations
    input: inst     --- instrument "i" or "s"
           <data_dir>/*_list_<inst>
    output: olist   --- a list of obsid
    """
    cmd   = 'ls ' + data_dir + '*_list_' + inst + ' > ' + zspace
    os.system(cmd)
    data  = mcf.read_data_file(zspace, remove=1)
    olist = []
    for ifile  in data:
        out   = mcf.read_data_file(ifile)
        for ent in out:
            atemp = re.split('\s+', ent)
            olist.append(atemp[1])

    olist = list(set(olist))
    olist = [int(x)  for x in olist]
#
#--- read past failed obsid lists
#
    ffile = house_keeping + 'skip_obsids'
    flist = mcf.read_data_file(ffile)
    flist = [int(x)  for x in flist]

    olist = olist + flist
    olist = list(set(olist))

    return olist

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    out = extract_calb_hz43()
    print("I AM HERE: " + str(out))
