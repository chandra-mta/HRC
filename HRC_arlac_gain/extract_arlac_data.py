#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#################################################################################
#                                                                               #
#   extract_arlac_data.py: extract arlac calibration data                       #
#                                                                               #
#   author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                               #
#   Last Update: Jul 22, 2019                                                   #
#                                                                               #
#################################################################################

import sys
import os
import string
import re
import math
import unittest
import time
import random
import numpy
import astropy.io.fits  as pyfits
from datetime import datetime
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts2/house_keeping/dir_list_py'

with  open(path, 'r') as f:
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

import mta_common_functions as mcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- extract_calb_arlac: create a list of calibration observations of ArLac from the database 
#-----------------------------------------------------------------------------------------

def extract_calb_arlac():
    """
    create a list of calibration observations of ArLac from the database
    input:  none
    output: <house_keeping>/hrc_i   --- a list of arlac observation on hrc i
    <house_keeping>/hrc_s   --- a list of arlac observation on hrc s
    chk --- if htere are new entries, return 1, otherwise 0
    """
#
#--- read the past data
#
    hrc_past = make_past_obsid_list()
#
#--- read database
#
    data = mcf.read_data_file('/data/mta4/obs_ss/sot_ocat.out')
#
#--- extract ArLac calibration data
#
    hrc_i = []
    hrc_s = []
    for ent in data:
        mc1  = re.search('ARLAC', ent)
        mc1a = re.search('ArLac', ent)
        mc1b = re.search('AR LAC',ent)
        mc2  = re.search('CAL',   ent)
        mc3  = re.search('archived',  ent)
        if (mc1 is not None) or (mc1a is not None) or (mc1b is not None):
            if mc2 is not None:
                if mc3 is not None:
                    atemp = re.split('\^', ent)
                    obsid = atemp[1].strip()
                    obsid = int(float(obsid))
                    inst  = atemp[12].strip()

                    if inst.lower() == 'hrc-i':
                        hrc_i.append(obsid)
                    elif inst.lower() == 'hrc-s':
                        hrc_s.append(obsid)
#
#--- check whether there are any new arlac calibration data
#
    hrc_past = set(hrc_past)
    hrc_i    = set(hrc_i)
    hrc_s    = set(hrc_s)
    idiff = []
    sdiff = []
    new   = [idiff, sdiff]
    chk   = 0
    idiff = list(hrc_i - hrc_past)
    sdiff = list(hrc_s - hrc_past)
    
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

def make_past_obsid_list():
    """
    create a list of obsids of the past observations
    input: <data_dir>/hrcf*data.gz
    output: olist   --- a list of obsid
    """
    cmd   = 'ls ' + data_dir + 'hrcf*_pha.dat* > ' + zspace
    os.system(cmd)
    flist = mcf.read_data_file(zspace, remove=1)

    olist = []
    for ent in flist:
        atemp = re.split('hrcf', ent)
        btemp = re.split('_pha', atemp[1])
        try:
            obsid = int(btemp[0])
        except:
            continue

        olist.append(obsid)

    olist = list(set(olist))
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

    out = extract_calb_arlac()
    print('DATA: ' + str(out))
