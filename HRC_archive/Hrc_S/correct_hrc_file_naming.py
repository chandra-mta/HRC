#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#           correct_hrc_file_naming.py: correct wrong naming of data files                      #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Apr 01, 2021                                                           #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import math
import time
import Chandra.Time
import astropy.io.fits  as pyfits
#
path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/house_keeping/dir_list'
#path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_I/house_keeping/dir_list'

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

import mta_common_functions as mcf
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#------------------------------------------------------------------------------------------------
#-- correct_hrc_file_naming: correct wrong naming of data files                               ---
#------------------------------------------------------------------------------------------------

def correct_hrc_file_naming():
    """
    correct wrong naming of data files
        example: hrcf59_evt2.fits ---> hrcf00059_evt2.fits
    input:  none, but read from /data/mta/<inst>/
    output: corrected file names
    """
    for inst in ['s', 'i']:
        cmd = 'ls -d /data/hrc/'+ inst + '/* > ' + zspace
        os.system(cmd)
    
        data = mcf.read_data_file(zspace, remove=1)
        obsid_list = []
        for ent in data:
            atemp = re.split('\/', ent)
            obsid = atemp[-1]
            if mcf.is_neumeric(obsid):
                correct_naming(obsid, inst)

#------------------------------------------------------------------------------------------------
#-- correct_naming: check secondary and analysis directories and correct wrongly named fits and par file
#------------------------------------------------------------------------------------------------

def correct_naming(obsid, inst):
    """
    check secondary and analysis directories and correct wrongly named fits and par file
    input:  obsid   --- obsid   
            inst    --- instrument. either "i" or "s"
    """
    cobsid = str(int(float(obsid)))
    if len(cobsid) == 5:
        return 

    lobsid = mcf.add_leading_zero(obsid, 5)
    
    for sdir in ['secondary', 'analysis', 'repro']:

        cmd = 'ls /data/hrc/' + inst  + '/' + lobsid + '/' + sdir + '/hrcf* >' + zspace
        os.system(cmd)

        data = mcf.read_data_file(zspace, remove=1)
        for ent in data:
            atemp = re.split('\/', ent)
            fname = atemp[-1]
            mc = re.search(lobsid, fname)
            if mc is not None:
                continue
            else:
                atemp = re.split('hrcf', fname)
                btemp = re.split('_',   atemp[1])
                sobs  = btemp[0]
                new   = fname.replace(sobs, lobsid)
                full  = '/data/hrc/' + inst + '/' + lobsid + '/' + sdir + '/' + new

                cmd = 'mv ' + ent + ' ' + full
                print("I AM HERE: " + full)
                #os.system(cmd)


#------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    correct_hrc_file_naming()
