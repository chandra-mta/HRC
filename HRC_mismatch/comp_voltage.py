#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#       comp_voltage.py: compare commanded HV setting and actual reading HV output and  #
#                        if there are any discrepancies, send out a warning email       #
#                                                                                       #
#               Note: run check_cmd_diff.py to extract needed hrc hk0 fits files        #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Mar 19, 2021                                                   #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import time
import Chandra.Time
import numpy
import astropy.io.fits  as pyfits
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project2/Script/house_keeping/dir_list'
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
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)
#
#--- col list
#
hk_cols = ['time', 'hvpsstat', 'imtpast', 'imhblv', 'imbpast', 'imhvlv', 'sptpast',\
           'sphblv', 'spbpast', 'sphvlv', 's1hvst', 's1hvlv', 's2hvst', 's2hvlv']
ref_val = [0, 0, 0, 127, 0, 127, 0, 127, 0, 127, 0, 127, 0, 127]
hv_list = ['sptpast', 'spbpast', 's1hvst', 'imtpast', 'imbpast', 's2hvst']
in_list = ['HRC-S',   'HRC-S',  'ANTI-CO', 'HRC-I',   'HRC-I',  'ANTI-CO']

#----------------------------------------------------------------------------------------
#-- comp_voltage: compare commanded HV setting and actual reading HV output           ---
#----------------------------------------------------------------------------------------

def comp_voltage():
    """
    compare commanded HV setting and actual reading HV output and
    if there are any discrepancies, send out a warning email
        Note: run check_cmd_diff.py to extract needed hrc hk0 fits files
    input:  none but check hrc hk0 fits files and read them in
    output: warning email sentout
    """
#
#--- remove a warning email file if it still exists
#
    ifile = exc_dir + 'voltage_cont'
    mcf.rm_files(ifile)
#
#--- read value correspoinding tables
#
    [hrc_i_cm, hrc_i_hv]   = corresponding_table('/hrc_i_list')
    [hrc_s_cm, hrc_s_hv]   = corresponding_table('/hrc_s_list')
    [antico_cm, antico_hv] = corresponding_table('/antico_list')
#
#--- find available hrc hk0 fits files
#
    cmd     = 'ls ./**_hk0.fits* > ' + zspace
    os.system(cmd)
    hk_list = mcf.read_data_file(zspace, remove=1)

    if len(hk_list) < 1:
        exit(1)
#
#--- read the data into a dictionary form
#--- ref_dict contains initial values <msid> <--> <value>, 
#--- but hk_dict is <msid> <---> <a list of data>
#
    hk_dict  = {}
    ref_dict = {}
#
#--- initialize
#
    for k in range(0, len(hk_cols)):
        col = hk_cols[k]
        hk_dict[col]  = numpy.array([])
        ref_dict[col] = ref_val[k]
#
#--- read hrc hk0 fits files and fill hk_dict
#
    for fits in hk_list:
        hout = pyfits.open(fits)
        data = hout[1].data

        for col in hk_cols:
#
#--- drop the first few entries since they are not stable
#
            out   = data[col][4:]
            psave = hk_dict[col]
            hk_dict[col] = numpy.append(psave, out)
        hout.close()

        mcf.rm_files(fits)
#
#--- start checking the data
#
    w_save   = [[] for i in range(0, len(hv_list))]
    t_list   = hk_dict['time']
    d_len    = len(t_list)
#
#--- check whether high voltages are on (from HVPSSTAT)
#
    hvpsstat = data['hvpsstat']
    for k in range(0, d_len):
        try:
            test = list(hvpsstat[k])
        except:
            continue

        if test == [1,1,1,1,1,1,1,1]:
            continue
#
#--- chechking hrc s hv
#
        if test[0] == 1:
            line = check_values(hk_dict, ref_dict, k, hrc_s_cm, hrc_s_hv, 'sptpast', 'sphblv')
            if line != '': w_save[0].append(line)

            line = check_values(hk_dict, ref_dict, k, hrc_s_cm, hrc_s_hv, 'spbpast', 'sphvlv')
            if line != '': w_save[1].append(line)
#
#--- checking antico hv: s1
        if test[2] == 1:
            line = check_values(hk_dict, ref_dict, k, antico_cm, antico_hv, 's1hvst', 's1hvlv')
            if line != '': w_save[2].append(line)
#
#--- checking hrc i hv
#
        if test[3] == 1:
            line = check_values(hk_dict, ref_dict, k, hrc_i_cm, hrc_i_hv, 'imtpast', 'imhblv')
            if line != '': w_save[3].append(line)

            line = check_values(hk_dict, ref_dict, k, hrc_i_cm, hrc_i_hv, 'imbpast', 'imhvlv')
            if line != '': w_save[4].append(line)
#
#--- checking antico hv: s2
#
        if test[5] == 1:
            line = check_values(hk_dict, ref_dict, k, antico_cm, antico_hv, 's2hvst', 's2hvlv')
            if line != '': w_save[5].append(line)
#
#--- update ref_dict
#
        for msid in hk_cols:
            ref_dict[msid] = hk_dict[msid][k]
#
#--- if there are any mismatches, create a warning email
#
    sline = ''
    for k in range(0, len(hv_list)):
        line = check_warning(w_save, hv_list, in_list, k)
        sline = sline + line

    if sline != '':
        ifile = exc_dir + 'voltage_cont'
        with open(ifile, 'w') as fo:
            fo.write(sline)

#----------------------------------------------------------------------------------------
#-- check_warning: check warning and create contents for email                         --
#----------------------------------------------------------------------------------------

def check_warning(w_save, hv_list, in_list, k):
    """
    check warning and create contents for email
    input:  w_save  --- a list of lists of warnings
            hv_list --- a list of msids
            in_list --- a list of instrument correspond to hv_list
            k       --- index of which msid etc are used
    output: line    --- waring setence for email (if there are any warnings)
    """
    line = ''
    if len(w_save[k]) > 0:
        line = line + 'Mismatch in ' + in_list[k] + ' ' + hv_list[k].upper() + ':\n\n'
        line = line + 'Time\t\t  Command (Expected)\tRecorded\n'
        for ent in w_save[k]:
            line = line + ent + '\n'

    return line

#----------------------------------------------------------------------------------------
#-- corresponding_table: read cm <--> hv correspoinding table                          --
#----------------------------------------------------------------------------------------

def corresponding_table(ifile):
    """
    read cm <--> hv correspoinding table
    input:  ifile   --- a name of file
    output: hrc_mc  --- cm list
            hrc_hv  --- hv list
    """
    ifile = house_keeping + ifile
    data  = mcf.read_data_file(ifile)
    hrc_cm = []
    hrc_hv = []
    for ent in data:
        atemp = re.split('\s+', ent)
        hrc_cm.append(int(float(atemp[0])))
        hrc_hv.append(int(float(atemp[1])))

    return [hrc_cm, hrc_hv]

#----------------------------------------------------------------------------------------
#-- check_values: check whether the values are in the range                            --
#----------------------------------------------------------------------------------------

def check_values(hk_dict, ref_dict, k, hrc_cm, hrc_hv, msid1, msid2):
    """
    check whether the values are in the range
    input:  hk_dict ---  a dictionary of data (msid <---> <a list of  all data>)
            ref_dict---  a dictionary of data (msid <---> <a current value>)
            k       ---  a position of the data to be checked in hk_dict[<msid>][k]
            hrc_cm  ---  a list of data value
            hrc_hv  ---  a list of data value correspoinding to hrc_cm
            msid1   ---  msid of value 1
            msid2   ---  msid of value 2
    output: line    ---  if the value is out of range, return information about that value
    """
    line = ''
    if(hk_dict[msid1][k] != ref_dict[msid1]) \
                    or (hk_dict[msid2][k] != ref_dict[msid2]):

        expected = find_value(hk_dict[msid1][k], hrc_cm, hrc_hv)
        range1   = expected - 1
        range2   = expected + 1
        if (hk_dict[msid2][k] < range1) or (hk_dict[msid2][k] > range2):
#
#--- convert time format 
#
            ltime = Chandra.Time.DateTime(hk_dict['time'][k]).date
            atemp = re.split('\.', ltime)       #--- remove fraction of seconds part
            ltime = atemp[0]

            line = line + ltime  + '\t' +  str(int(hk_dict[msid1][k])) 
            line = line + ' (' + str(expected)+ ')\t\t'+ str(hk_dict[msid2][k]) 

    return line

#----------------------------------------------------------------------------------------
#-- find_value: find a correspoinding value                                            --
#----------------------------------------------------------------------------------------

def find_value(dval, o_list, e_list):
    """
    find a correspoinding value
    input:  dval    --- the value to be converted
            o_list  --- a list of input values
            e_list  --- a list of output values
    outpu:  a selected corresponding value
    """
    chk = 0
    for k in range(0, len(o_list)):
        if dval == o_list[k]:
            pos = k
            chk = 1
            break
    if chk == 0:
        return 127

    else:
        return e_list[k]

#----------------------------------------------------------------------------------------

if __name__ == '__main__':

    comp_voltage()
