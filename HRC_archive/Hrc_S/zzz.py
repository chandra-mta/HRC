#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Aug 06, 2019                                                           #
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
#--- from ska
#
from Ska.Shell import getenv, bash
#
#--- set ciao environment 
#
ciaoenv  = getenv('source /soft/ciao/bin/ciao.csh; \
                   source /home/mta/bin/reset_param; setenv PFILES "${PDIRS}";\
                   set path=(/soft/ciao/bin/ $path);', shell='tcsh')
#
#--- reading directory list
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
import fits_operation       as mfo
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#------------------------------------------------------------------------------------------------
#-- find_hrc_calib_obsid: find obsids of hrc calibraiton data starting 61* or 62*              --
#------------------------------------------------------------------------------------------------

def find_hrc_calib_obsid(inst):
    """
    find obsids of hrc calibraiton data starting 61* or 62*
    input:  inst    --- s: hrc-s or i: hrc-i
    output: a list of obsids
    """
#
#--- create a list of already processed data
#
    cmd = 'ls -d /data/hrc/' + str(inst) + '/61* >' + zspace
    os.system(cmd)
    cmd = 'ls -d /data/hrc/' + str(inst) + '/62* >' + zspace
    os.system(cmd)

    data = mcf.read_data_file(zspace, remove=1)
    prev_list = []
#    for ent in data:
#        atemp = re.split('\/', ent)
#        prev_list.append(int(float(atemp[-1])))

#
#--- find today's date and set checking range for the last 30 days
#
    today = time.strftime('%Y:%j:%H:%M:%S', time.gmtime())
    today = int(Chandra.Time.DateTime(today).secs)
    start = today - 10 * 86400
#
#--- extract hrc obsid information
#
    line = 'operation=browse\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'filetype=evt1\n'
    line = line + 'tstart=' + str(start) + '\n'
    line = line + 'tstop='  + str(today) + '\n'
    line = line + 'go\n'

    with open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > ' + zspace
    os.system(cmd)

    mcf.rm_files('./zline')

    data = mcf.read_data_file(zspace, remove=1)
#
#--- select obsids with 61* and 62* starting
#
    h_list = []
    for ent in data:
        mc = re.search('hrcf', ent)
        if mc is not None:
            atemp = re.split('hrcf', ent)
            btemp = re.split('_', atemp[1])
            obsid = int(float(btemp[0]))
            if obsid > 61000 and obsid < 63000:
#
#--- if it is already observed skip it
#
                if obsid in prev_list:
                    continue
#
#--- check which instrument
#
                chk = check_inst(obsid)
                if chk == inst:
                    h_list.append(obsid)

    return h_list

#------------------------------------------------------------------------------------------------
#-- check_inst: find which instremnt this obsid based on                                       --
#------------------------------------------------------------------------------------------------

def check_inst(obsid):
    """
    find which instremnt this obsid based on
    input:  obsid   --- obsid
    output: inst    --- instrument i or s
    """
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'filetype=evt1\n'
    line = line + 'obsid='  + str(obsid) + '\n'
    line = line + 'go\n'

    with open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline >' + zspace 
    os.system(cmd)

    data = mcf.read_data_file(zspace, remove=1)
    dfile = ''
    for ent in data:
        mc = re.search('hrcf', ent)
        if mc is not None:
            dfile = ent

    if dfile == '':
        return 'na'

    flist = pyfits.open(dfile)
    try:
        inst  = flist[0].header['DETNAM']
    except:
        inst  = flist[1].header['DETNAM']
    flist.close()

    mc   = re.search('HRC-I', inst)
    if mc is not None:
        inst = 'i'
    else:
        inst = 's'

    mcf.rm_file(dfile)

    return inst



#------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    inst = 's'
    out = find_hrc_calib_obsid(inst)
    print("I AM HERE: " + str(out))
