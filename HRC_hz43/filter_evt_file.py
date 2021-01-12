#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#################################################################################
#                                                                               #
#           filter_evt_file.py: filter the events with dtf and status           #
#                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                           #
#                                                                               #
#           Last Update: Jul 18, 2019                                           #
#                                                                               #
#################################################################################

import sys
import os
import string
import re
import copy
import math
import numpy
import astropy.io.fits  as pyfits
import time
import random
#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; \
                   source /home/mta/bin/reset_param ', shell='tcsh')
ciaoenv  = getenv('source /soft/ciao/bin/ciao.sh')
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
#
#--- append path to a private folders
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
sys.path.append(hrc_common)

import mta_common_functions as mcf
import hrc_common_functions as hcf

rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

cos45  =  0.707106781        #---- cos(45 deg)

#-------------------------------------------------------------------------------------------
#-- filter_evt_file: fiter the events in the fits file by dtf and status                  --
#-------------------------------------------------------------------------------------------

def filter_evt_file(fits, outfile='./cleaned_file.fits', obsid = ''):
    """
    fiter the events in the fits file by dtf and status
    input:  fits    --- fits file
    output: fits    --- filtered fits file
    """
#
#--- get obsid
#
    if obsid == '':
        atemp = re.split('hrcf', fits)
        btemp = re.split('_',    atemp[1])
        obsid = btemp[0]
#
#--- copy the fits file and make sure that it is not zipped
#
    mc = re.search('gz', fits)
    if mc is not None:
        cmd = 'cp ' + fits + ' ztemp1.fits.gz'
        os.system(cmd)
        cmd = 'gzip -d ztemp1.fits.gz'
        os.system(cmd)
    else:
        cmd = 'cp ' + fits + ' ztemp1.fits'
        os.system(cmd)
#
#--- filter out dtf <= 0.98
#
    [pstart, pstop] = get_dead_period(obsid)
    
    t      = pyfits.open('ztemp1.fits')
    tdata  = t[1].data
    t_list = tdata['time']
    t.close()
#
#--- if there are clipped parts based of dtf, create mask for that
#
    clen = len(pstart)
    if clen > 0:
        m    = 0
        mask = []
        cend = 0
        for t in t_list:
            if cend == 0:
                for k in range(m, clen):
                    if t < pstart[k]:
                        mask.append(True)
                        break
                    elif (t >= pstart[k]) and (t <= pstop[k]):
                        mask.append(False)
                        break
                    elif t > pstop[k]:
                        if k == clen -1:
                            cend = 1
                            break
                        m = k  
                        continue
            else:
                mask.append(True)
        mask.append(False)
    
        mask = numpy.array(mask)
    
        newtdata = tdata[mask]
#
#--- remove the masked part and create a new fits file
#
        hdu = pyfits.BinTableHDU(data=newtdata)
        hdu.writeto('ztemp2.fits')
    
        cmd = 'mv -f ztemp2.fits ztemp1.fits'
        os.system(cmd)
#
#--- filter by status
#
    filter_by_status('ztemp1.fits')

    cmd = 'mv -f  ztemp1.fits ' + outfile
    os.system(cmd)

#-------------------------------------------------------------------------------------------
#-- get_dead_period: extract time periods which is dtf <= 0.98                            --
#-------------------------------------------------------------------------------------------

def get_dead_period(obsid):
    """
    extract time periods which is dtf <= 0.98
    input: obsid    --- obsid
    output: pstart  --- starting time
            pstop   --- stoppping time
                          these are +/- 2 of the center value listed on dtf fits file
    """
#
#--- get hrc dtf fits file
#
    dtf   = hcf.run_arc5gl(start='', stop='', obsid=obsid, operation='retrieve',\
                           level ='1', filetype='dtf')
#
#--- read out time and dtf column from the fits file
#
    tdata = hcf.read_fits_file(dtf)
    tout  = hcf.extract_col_data(tdata[1], ['time', 'dtf'])

    mcf.rm_files(dtf)
#
#--- find the period which match the condition (dtf<=0.98)
#
    pstart = []
    pstop  = []
    for k in range(0, len(tout)):
        try:
            cdtf = float(tout[k][1])
        except:
            continue
        if cdtf <= 0.98:
            time = int(float(tout[k][0]))
#
#--- take the period to be +/- 2 of the center value listed
#
            pstart.append(time - 2)
            pstop.append(time  + 2)

    return [pstart, pstop]

#-------------------------------------------------------------------------------------------
#-- filter_by_status: filter fits file by the status column                               --
#-------------------------------------------------------------------------------------------

def filter_by_status(fits):
    """
    filter fits file by the status column
    input:  fits    --- fits file
            filter = '0000xxxx000xxxxx'
    output: fits    --- status filtered fits file
    """
    cmd  = ' dmcopy "' + fits + '[status=0000xxxx000xxxxx]" outfile=zxc.fits clobber=yes'
#
#--- dmcopy may not work with ascds; if so, use ciao
#
    hcf.run_ascds(cmd)
    if not os.path.isfile('zxc.fits'):
        hcf.run_ciao(cmd)

    cmd  = 'mv zxc.fits ' + fits
    os.system(cmd)

#----------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) == 2:
        fits = sys.argv[1]
        fits.strip()

        filter_evt_file(fits)
    else:
        print("Need fits file name!")
