#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#################################################################################################
#                                                                                               #
#   analyze_hz43_data.py: extract data, analyze and update data tables                          #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Jul 18, 2019                                                           #
#                                                                                               #
#################################################################################################

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
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

with open(path, 'r') as f:
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
sys.path.append(hrc_common)

import mta_common_functions as mcf
import hrc_common_functions as hcf
import adjust_hz43_position as ahp
import find_hz43_data       as fhd
import extract_hz43_stat    as ehs
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

olist  = ['hrc_i_coords', 'hrc_s_coords']

#-----------------------------------------------------------------------------------------
#-- analyze_hz43_data: extract data, analyze and update data tables                    ---
#-----------------------------------------------------------------------------------------

def analyze_hz43_data():
    """
    extract data, analyze and update data tables
    input:  none but download data with arc5gl
    output: <house_keeping>/hrc_i_coords, <hosue_keeping>/hrc_s_coords
            both has <obsid> <ra> <dec> <skyx> <skyy> 
            <data_dir>/pi_* and samp_*
    """
#
#--- find new HZ43 observations data_set: [hrc-i, hrc-s]
#
    data_set = fhd.extract_calb_hz43()

    if (len(data_set[0]) == 0) and (len(data_set[1]) == 0):
        print("No new HZ43 data found")
        exit(1)

    for k in range(0, 2):
        data_list = data_set[k]
        if len(data_list) == 0:
            continue
        if k == 0:
            print("Analyzing HRC I")
        else:
            print("Analyzing HRC S")

        zline = ''
        for obsid in data_list:
            print("OBSID: " + str(obsid))
#
#--- extract evt1 data
#
            evt1 = hcf.run_arc5gl(0, 0,  obsid = obsid, operation='retrieve',\
                                  level ='1',   filetype='evt1')
            if evt1 == "":
                print("No evt1 file extracted")
                continue
#
#--- get the ra/dec of HZ43 adjusted to proper motion
#
            fout      = pyfits.open(evt1)
            fhead     = fout[1].header
            otime     = fhead['DATE-OBS']
            [ra, dec] = ahp.adjust_hz43_position(otime)
#
#--- convert coordinates from cel to sky
#
            cmd = 'dmcoords ' + evt1 + ' opt=cel ra=' + str(ra) + ' dec=' + str(dec) 
            cmd = cmd + ' verbose=1 > ' + zspace
            hcf.run_ciao(cmd)
#
#--- extract sky coordindats
#
            info = mcf.read_data_file(zspace, remove=1)
            out  = get_skyp_coord(info)
            if out == False:
                continue

            skyx  = out[0]
            skyy  = out[1]
            
            zout  = str(obsid) + '\t' + str(ra) + '\t' + str(dec) + '\t' 
            zout  = zout       + str(skyx) + '\t' + str(skyy) + '\n'
            zline = zline + zout
#
#--- upate the data table
#
            ehs.extract_hz43_stat(obsid,  evt1)

        if zline != "":
            ofile = house_keeping + olist[k]
            with open(ofile, 'a') as fo:
                fo.write(zline)
#
#--- reorder data file in time order
#
        order_by_time()
#
#--- clean up coord lists
#
        coords_list(olist[0])
        coords_list(olist[1])

#-----------------------------------------------------------------------------------------
#-- order_by_time: re-order data file by time                                           --
#-----------------------------------------------------------------------------------------

def order_by_time():
    """
    re-order data file by time
    input: none but read from <data_dir>/*_list*
    output: re-rodered data file
    """
    cmd = 'ls ' +  data_dir + '*_list_* > ' + zspace
    os.system(cmd)

    data = mcf.read_data_file(zspace, remove=1)

    for ifile in data:

        out   = mcf.read_data_file(ifile)
        edict = {}
        elist = []
        for  ent in out:
            atemp = re.split('\s+', ent)
            atime = int(float(atemp[0]))
            edict[atime] = ent
            elist.append(atime)
#
#--- remove duplicate entry and reorder
#
        elist = list(set(elist))
        elist = sorted(elist)

        line  = ''
        for otime in elist:
            line = line + edict[otime] + '\n'

        with open(ifile, 'w') as fo:
            fo.write(line)

#-----------------------------------------------------------------------------------------
#-- coords_list: remove duplicates from hrc_<*>_coords list                             --
#-----------------------------------------------------------------------------------------

def coords_list(alist):
    """
    remove duplicates from hrc_<*>_coords list
    input:  alist                   --- a name of the list
    output: <hosue_keeping>/alist   --- a cleaned up list
    """
    ifile = house_keeping + alist
    data  = mcf.read_data_file(ifile)
    dlist = []
    ddict = {}
    for ent in data:
        atemp = re.split('\s+', ent)
        obsid = int(atemp[0])
        dlist.append(obsid)
        ddict[obsid] = ent

    dlist = list(set(dlist))

    line  = ''
    for ent in dlist:
        line = line + ddict[ent] + '\n'

    with open(ifile, 'w') as fo:
        fo.write(line)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def get_skyp_coord(info):

    for ent in info:
        mc = re.search('SKY', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            skyx  = float(atemp[1])
            skyy  = float(atemp[2])

            return [skyx, skyy]

    return False

#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    analyze_hz43_data()

