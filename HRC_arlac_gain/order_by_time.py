#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################################
#                                                                                               #
#           order_by_time.py: order the data by date                                            #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Jan 25, 2021                                                           #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import time
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts2/house_keeping/dir_list_py'

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

import mta_common_functions as mcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- order_by_time: order the data by date                                               --
#-----------------------------------------------------------------------------------------

def order_by_time():
    """
    order the data by date
    input: none but read <data_dir>/fitting_results
    output: re-ordered <data_dir>/fitting_results
    """
    ifile = data_dir + 'fitting_results'

    data  = mcf.read_data_file(ifile)
    ddict = {}
    dlist = []
    for ent in data:
        atemp = re.split('\s+', ent)
        try:
            stime = int(float(atemp[2]))
        except:
            continue
        ddict[stime] = ent
        dlist.append(stime)
#
#--- sort the time and update the data file
#
    dlist = sorted(dlist)
    line  = ''
    for ent in dlist:
        line = line + ddict[ent] + '\n'

    with open(ifile, 'w') as fo:
        fo.write(line)


#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    order_by_time()

