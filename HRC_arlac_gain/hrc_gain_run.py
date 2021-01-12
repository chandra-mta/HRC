#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#####################################################################################################
#                                                                                                   #
#       hrc_gain_fit_voigt.py: extract hrc evt2 files and fit a normal distribution on pha values   #
#                                                                                                   #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                                #
#                                                                                                   #
#           Last Update:    Sep 05, 2019                                                            #
#                                                                                                   #
#####################################################################################################

import os
import sys
import re
import string
import random
import math
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
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import mta_common_functions     as mcf
import hrc_gain_fit_voigt       as hgfv  
import hrc_gain_trend_plot      as hgtp
import extract_arlac_data       as ead
import order_by_time            as obt
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

#---------------------------------------------------------------------------------------------------
#-- hrc_gain_run: extract hrc evt2 file, find the brightest object and create pha distribution     -
#---------------------------------------------------------------------------------------------------

def  hrc_gain_run(c_input):
    """
    extract hrc evt2 file, find the brightest object and create pha distribution
    this is a control script

    Input:  c_inut      --- if it is obsid, use it as a input
                            otherwise, a list of new candidates will be create based database
    Output: <header>_pha.dat    --- pha distribution data
            <header>_gfit.png   --- a plot of pha data
            fitting_results     --- a table of fitted results
    """
#
#--- if an obsid is provided, analyize that, else get new obsids from databases
#
    if mcf.is_neumeric(c_input):
        candidate_list = [c_input]

    elif isinstance(c_input, list):
        candidate_list = c_input

    else:
        out = ead.extract_calb_arlac()
        candidate_list = out[0] + out[1]
#
#--- run the main script only when candidates exist
#
    if len(candidate_list) > 0:
#
#--- analyze data
#
        hgfv.hrc_gain_fit_voigt(candidate_list)
#
#--- create plots
#
        obt.order_by_time()
        hgtp.hrc_gain_trend_plot()

        send_notification()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def send_notification():

    text = '\nNew AR Lac Observations are added to:\n\n'
    text = text + 'https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/arlac_energy_trend.html'
    text = text + '\n'
    
    with open(zspace, 'w') as fo:
        fo.write(text)
    
    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" vkashyap@cfa.harvard.edu'
    os.system(cmd)
    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" ' + admin
    os.system(cmd)
    
    mcf.rm_files(zspace)


#--------------------------------------------------------------------

if __name__ == '__main__':
#
#--- you can give an obsid of ar lac as an argument
#
    if len(sys.argv) >= 2:
        c_input = sys.argv[1]
        c_input.strip()
    else:
        c_input = ""

    hrc_gain_run(c_input)
