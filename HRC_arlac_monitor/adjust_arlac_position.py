#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################################
#                                                                                               #
#           adjust_arlac_position.py: adjust ArLac position with proper motion                  #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Jan 21, 2021                                                           #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import math
import time
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts/house_keeping/dir_list'

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
sys.path.append(bin_dir)
sys.path.append(hrc_common)

import hrc_common_functions as hcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- adjust_arlac_position: adjust ArLac position with proper motion                       --
#-----------------------------------------------------------------------------------------

def adjust_arlac_position(date):
    """
    adjust ArLac position with proper motion
            22 08 40.8180372773::+45 44 32.104045703::-52.190::47.190
    input:  date    in <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
                        if it is not given, no corrections for the propoer motion
    output: coordinate  [ra, dec] in decimal format.
    """
#
#--- evt1 DATE-OBS to fractional year
#
    cyear = hcf.convert_time_to_fyear(date)
#
#--- set several initial values
#
    pfile = house_keeping + 'arlac_pos'
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
        date = sys.argv[1].strip()

        [ra, dec]  = adjust_arlac_position(date)
        print('RA: ' + str(ra) + '   DEC: ' + str(dec))

    else:
        print("provide year --- fractional year is fine")

