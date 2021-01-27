#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#########################################################################################
#                                                                                       #
#           extract_vega_data.py: extract vega calibration data                         #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Jan 27, 2021                                                   #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import math
import unittest
import time
import random
import astropy.io.fits  as pyfits
from datetime import datetime
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/Vega/Scripts/house_keeping/dir_list'

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

#-----------------------------------------------------------------------------------------
#-- extract_calb_vega: create a list of calibration observations of Vega from the database 
#-----------------------------------------------------------------------------------------

def extract_calb_vega():
    """
    create a list of calibration observations of Vega from the database
    input:  none
    output: <house_keeping>/hrc_i   --- a list of vega observation on hrc i
            <house_keeping>/hrc_s   --- a list of vega observation on hrc s
            chk                     --- if htere are new entries, return 1, otherwise 0
    """
#
#--- read the past data
#
    ifile_i = house_keeping + 'hrc_i_list'
    hrc_i_p = mcf.read_data_file(ifile_i)
    hrc_i_p = [int(x) for x in hrc_i_p]

    ifile_s = house_keeping + 'hrc_s_list'
    hrc_s_p = mcf.read_data_file(ifile_s)
    hrc_s_p = [int(x) for x in hrc_s_p]
#
#--- read database
#
    data = mcf.read_data_file('/data/mta4/obs_ss/sot_ocat.out')
#
#--- extract Vega calibration data
#
    hrc_i = []
    hrc_s = []
    for ent in data:
        mc1  = re.search('VEGA',     ent)
        mc1a = re.search('Vega',     ent)
        mc2  = re.search('CAL',      ent)
        mc3  = re.search('archived', ent)
        mc4  = re.search('HRC-I',    ent)
        mc5  = re.search('HRC-S',    ent)
        if (mc1 is not None) or (mc1a is not None):
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
#--- check whether there are any new vega calibration data
#
    hrc_i = set(hrc_i)
    hrc_s = set(hrc_s)
    chk = 0
    if hrc_i != set(hrc_i_p):
        chk = 1
    if hrc_s != set(hrc_s_p):
        chk = 1
#
#--- if so, update the list
#
    if chk > 0:
        with open(ifile_i, 'w') as fo:
            for ent in list(hrc_i):
                fo.write(str(ent) + '\n')
    
        with open(ifile_s, 'w') as fo:
            for ent in list(hrc_s):
                fo.write(str(ent) + '\n')

    chk_save = house_keeping + 'chk_save'
    with open(chk_save, 'w') as fo:
        fo.write(str(chk) + '\n')

    return chk

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    extract_calb_vega()
