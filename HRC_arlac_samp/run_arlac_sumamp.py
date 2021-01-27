#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#####################################################################################
#                                                                                   #
#       run_arlac_sumamp.py: find new ar lac observations and update the page       #
#                                                                                   #
#           author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                   #
#           Last Update: Jan 26, 2021                                               #
#                                                                                   #
#####################################################################################

import sys
import os
import string
import re
import time
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts3/house_keeping/dir_list'

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

import mta_common_functions     as mcf
import extract_arlac_stat       as eas
import create_arlac_trend_plots as catp
import create_count_rate_plots  as ccrp
import create_profile_plot      as cpp
import create_html_page         as chp

import random
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

admin  = 'tisobe@cfa.harvard.edu'

#-----------------------------------------------------------------------------------------
#-- run_arlac_sumamp: find new ar lac observations and update the page                  --
#-----------------------------------------------------------------------------------------

def run_arlac_sumamp():
    """
    find new ar lac observations and update the page
    input: none
    output: updated <web_dir>/arlac_trends.html
    """
#
#--- check whether there are new Al Lac observation
#
    dout = find_new_arlac_observations(ret=1)
#
#--- if there is no data, stop
#
    if (len(dout[0]) == 0) and (len(dout[1]) == 0):
        exit(1)
    else:
#
#---- run the processes
#
        eas.run_process(dout)               #--- update data tables
        order_data_file()                   #--- order the data tables in time order
        catp.create_plots()                 #--- create binned trend plots
        ccrp.create_count_rate_plots('')    #--- create count rate and background rate plots
        cpp.run_for_all_profile_plot()      #--- fit gamma profile on the data and creat a plot
        chp.create_html_page('')            #--- update Ar Lac SumAmp page

        cmd = 'chgrp -R mtagroup /proj/web-cxc/htdocs/contrib/cxchrc/HRC_trendings/ArLac/*'
        os.system(cmd)
        cmd = 'chgrp -R mtagroup /data/aschrc6/wilton/isobe/Project8/ArLac/Data3/*'
        os.system(cmd)
#
#--- send out notificaiton
#
        send_notification()

#-----------------------------------------------------------------------------------------
#-- find_new_arlac_observations: create a list of calibration observations of ArLac from the database 
#-----------------------------------------------------------------------------------------

def find_new_arlac_observations(ret='0', nk=1):
    """
    create a list of calibration observations of ArLac from the database
    input:  ret --- if 0, print out the result. if >0, print out and return the results
            nk  --- if 0, don't read skip list, otherwise, read skip list
                        skip list contains obsid which we could not analyze in past
    output: <house_keeping>/hrc_i   --- a list of arlac observation on hrc i
            <house_keeping>/hrc_s   --- a list of arlac observation on hrc s
            [<hrc_i list>, <hrc_s list>] --- a list of lists of  new hrc i and hrc s obsids
    """
#
#--- read obsids which we could not process in past
#
    if nk == 0:
        s_list = []
    else:
        ifile   = house_keeping + 'skip_list'
        s_list  = mcf.read_data_file(ifile)
        s_list  = [int(x) for x in s_list]
#
#--- find the past hrc i observations
#
    ifile   = data_dir + 'pi_list_i'
    hrc_i_p = set(get_past_hrc_list(ifile) + s_list)
#
#--- find the past hrc s observations
#
    ifile   = data_dir + 'pi_list_s'
    hrc_s_p = set(get_past_hrc_list(ifile) + s_list)
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
        mc1  = re.search('ARLAC',     ent, re.IGNORECASE)
        mc1a = re.search('ArLac',     ent, re.IGNORECASE)
        mc1b = re.search('AR LAC',    ent, re.IGNORECASE)
        mc2  = re.search('CAL',       ent, re.IGNORECASE)
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
    hrc_i = set(hrc_i)
    hrc_s = set(hrc_s)

    chk = 0
    if hrc_i != set(hrc_i_p):
        chk = 1
    if hrc_s != set(hrc_s_p):
        chk = 1

    if chk == 0:
        if ret == 1:
            return [[], []]
    else:
        diff_i = hrc_i - hrc_i_p
        diff_s = hrc_s - hrc_s_p

        if len(diff_i) > 0:
            ofile = house_keeping + 'hrci_new'
            with open(ofile, 'w') as fo:
                for ent in diff_i:
                    fo.write(str(ent) + '\n')

        if len(diff_s) > 0:
            ofile = house_keeping + 'hrcs_new'
            with open(ofile, 'w') as fo:
                for ent in diff_s:
                    fo.write(str(ent) + '\n')

        if ret == 1:
            return [list(diff_i), list(diff_s)]
        
#-----------------------------------------------------------------------------------------
#-- get_past_hrc_list: create a list of obsids from the processed data list file        --
#-----------------------------------------------------------------------------------------

def get_past_hrc_list(ifile):
    """
    create a list of obsids from the processed data list file
    input:  ifile   --- input file name
    output: hrc_p   --- a set of obsids
    """
    h_out   = mcf.read_data_file(ifile)

    hrc_p = []
    for ent in h_out:
        atemp = re.split('\s+', ent)
        hrc_p.append(int(float(atemp[1])))

    return hrc_p

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def send_notification():

    text = '\nNew AR Lac Observations are added to:\n\n'
    text = text + 'https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/arlac_trends.html'
    text = text + '\n'
    
    with open(zspace, 'w') as fo:
        fo.write(text)
    
    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" vkashyap@cfa.harvard.edu'
    ###os.system(cmd)
    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" ' + admin
    os.system(cmd)
    
    mcf.rm_files(zspace)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def order_data_file():

    for ifile in ['pi_list_i', 'pi_list_s', 'samp_list_i', 'samp_list_s']:
        dfile = data_dir + ifile
        data  = mcf.read_data_file(dfile)
        dlist = []
        ddict = {}
        for ent in data:
            atemp = re.split('\s+', ent)
            atime = int(float(atemp[0]))
            dlist.append(atime)
            ddict[atime] = ent

        dlist = sorted(set(dlist))

        line  = ''
        for atime in dlist:
            line = line + ddict[atime] + '\n'

        with open(dfile, 'w') as fo:
            fo.write(line)

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    run_arlac_sumamp()
