#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#############################################################################################
#                                                                                           #
#       get_yearly_evt_count.py: counts the numbers of events in the yearly evt1 files      #
#                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                           #
#           last update: Mar 19, 2021                                                       #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import random
import time
import datetime
import math
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project1/Scripts/house_keeping/dir_list'

with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
#
#--- append a path to a private folder to python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- converTimeFormat contains MTA time conversion routines
#
import mta_common_functions         as mcf
import hrc_stowed_common_function   as hscf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

inst_list  = ('Hrc_i_115', 'Hrc_s_125', 'Hrc_s_125_hi')

#---------------------------------------------------------------------
#-- get_yearly_evt_count: counts the numbers of events in the yearly evt1 files
#---------------------------------------------------------------------

def get_yearly_evt_count():
    """
    counts the numbers of events in the yearly evt1 files
    input: none, but read from yearly evt1.fits files
    output: <data_dir>/<inst>/<inst>_yearly_evt_counts
    """
#
#--- find today's date
#
    out   = time.localtime()
    lyear = out.tm_year
    mon   = out.tm_mon
#
#--- if the date is the first half of the year, compute the last year
#
    if mon > 6:
        lyear += 1
#
#--- check each instrument
#
    for inst in inst_list:

        line = ''
        for year in range(2000, lyear):
            infile = data_dir + inst + '/'+ inst.lower() + '_' + str(year) + '_evt1.fits.gz'
#
#--- read the numbers of events
#
            try:
                out  = hscf.fitsTableStat(infile, 'time')
                cnt  = out[-1]
            except:
                cnt  = 0


            line = line +  str(year) + '\t' + str(cnt) + '\n'

        ofile = data_dir + inst + '/' + inst.lower() + '_yearly_evt_counts'
        with open(ofile, 'w') as fo:
            fo.write(line)

#---------------------------------------------------------------------

if __name__ == "__main__":

        get_yearly_evt_count()

