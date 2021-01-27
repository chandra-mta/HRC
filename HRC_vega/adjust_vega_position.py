#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################################
#                                                                                               #
#           adjust_vega_position.py: adjust Vega position with proper motion                    #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Jan 27, 2021                                                           #
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

import mta_common_functions as mcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- adjust_vega_position: adjust Vega position with proper motion                       --
#-----------------------------------------------------------------------------------------

def adjust_vega_position(cyear=''):
    """
    adjust Vega position with proper motion
    input:  cyear   --- the year of the observed. it can be in year date. 
                        if it is not given, no corrections for the propoer motion
    output: coordinate  [ra, dec] in decimal format.
    """
    if not mcf.is_neumeric(cyear):
        cyear = 2000.0
#
#--- set several initial values
#
    pfile = house_keeping + 'vega_pos'
    with open(pfile, 'r') as f:
        out = f.read()
    atemp = re.split('::', out)
    tra   = float(atemp[0])
    tdec  = float(atemp[1])
    pra   = float(atemp[2]) /3600.0 /1000.0
    pdec  = float(atemp[3]) /3600.0 /1000.0

    dyear = float(cyear) - 2000.0

    tra  += dyear * pra
    tdec += dyear * pdec
    
    return [tra, tdec]

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        year   = float(sys.argv[1].strip())

        [ra, dec]  = adjust_vega_position(year)
        print('RA: ' + str(ra) + '   DEC: ' + str(dec))

    else:
        print("provide year --- fractional year is fine")

