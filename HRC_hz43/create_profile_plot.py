#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#############################################################################################
#                                                                                           #
#           create_profile_plot.py: fit gamma profile on the data and creat a plot          #
#                                                                                           #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                        #
#                                                                                           #
#           Last Update:    Jan 21, 2021                                                    #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import random
import math
import numpy
import time
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

with  open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var   = atemp[1].strip()
    line  = atemp[0].strip()
    exec("%s = %s" %(var, line))
html_top = html_top.replace('#', ':')
#
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import mta_common_functions       as mcf  #---- contains other functions commonly used in MTA scripts
import gamma_function             as gamma
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)
#
head_list = ['samp_center_list', 'pi_center_list', 'samp_p_list_',\
            'samp_n_list_', 'pi_p_list_', 'pi_n_list_']

#---------------------------------------------------------------------------------------------------
#-- run_for_all_profile_plot: read which obsids are available for profile plotting and plot them   -
#---------------------------------------------------------------------------------------------------

def run_for_all_profile_plot():
    """
    read which obsids are available for profile plotting and plot them
    input: none
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/pi_p_list_<#>_vifits.png etc
    """
#
#--- check whether there are new data
#
    chk_save = 0
    try:
        c_file = house_keeping + 'chk_save'
        with open(c_file, 'r') as f:
            chk_save = int(f.read())
    except:
        pass

    if chk_save == 0:
        exit(1)

    sdir  = data_dir + 'Fittings/'
    cmd   = 'ls -d ' + sdir + '/*  > ' + zspace
    os.system(cmd)

    olist = mcf.read_data_file(zspace, remove=1)

    for ent in olist:
        atemp = re.split('\/', ent)
        obsid = atemp[-1]
#
#--- check the data are already processed
#
        if check_plot_exist(obsid) == 1:
            continue

        print('Creating fitting plot of ObsID: ' + str(obsid))
        create_profile_plot(obsid)

#---------------------------------------------------------------------------------------------------
#-- create_profile_plot: fit gamma profile on the data and creat a plot                          ---
#---------------------------------------------------------------------------------------------------

def  create_profile_plot(obsid):
    """
    fit gamma profile on the data and creat a plot
    input:  obsid   --- obsid. the data is read from <data_dir>
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/<head>_<col #>_vfits.png
    """
    sdir = data_dir + 'Fittings/' + str(obsid) + '/'

    for m in range(0, len(head_list)):
        sline = ''
        if m < 2:
            fname = sdir + head_list[m]
            sline  = create_plot_and_table_entry(head_list[m], '', obsid,  fname)
            if (sline == False) or (sline == ''):
                continue
        else:
            for k in range(0, 20):
                fname = sdir + head_list[m] + str(k)
                out = create_plot_and_table_entry(head_list[m], k, obsid,  fname)
                if (out == False) or (out == ''):
                    continue
                else:
                    sline = sline + out

        if sline != '':
            out  = sdir + head_list[m] +  'fit_results'
            with  open(out, 'w') as fo:
                fo.write(sline)
            
#---------------------------------------------------------------------------------------------------
#-- create_plot_and_table_entry: fit a gamma profile on the data, 
#-- create a plot, and return the fitting results 
#---------------------------------------------------------------------------------------------------

def create_plot_and_table_entry(head, pos, obsid, fname):
    """
    fit a gamma profile on the data, create a plot, and return the fitting results
    input:  head    --- head part of the file
            pos     --- position of the data. either center of cell <#>
            obsid   --- obsid
            fname   --- the data file name
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/<head>_<col #>_vfits.png
            line    --- fitted results in a table column form
    """

    if not  os.path.isfile(fname):
        return False
    try:
        out = binning_data(fname)
        if out == False:
            return False
        else:
            [abin, data, avg, std] = out
    except:
        return False
#
#--- fit a gamma distribution on the data
#
    title = 'ObsID: ' + str(obsid)
    [a, b, h]  = gamma.fit_gamma_profile(abin, data, avg, std, plot_title=title)
#
#--- rename a plotting file
#            
    outdir  = web_dir + 'Plots/Indivisual_Plots/' + str(obsid) + '/' 

    if not os.path.isdir(outdir):
        cmd = 'mkdir ' + outdir
        os.system(cmd)

    if not os.path.isdir(outdir):
        cmd = 'mkdir ' + outdir
        os.system(cmd)

    outfile = outdir + head + str(pos) + '_vfit.png'
    cmd = 'mv out.png ' + outfile
    os.system(cmd)
#
#--- keep the fitting results
#
    line = str(pos) + '\t'
    line = line + str(a) + '\t'
    line = line + str(b) + '\t'
    line = line + str(h) + '\n'

    return line

#---------------------------------------------------------------------------------------------------
#-- binning_data: binning the data into 256 bins                                                 ---
#---------------------------------------------------------------------------------------------------

def binning_data(fname):
    """
    binning the data into 256 bins
    input:  fname   --- data file name
    output: abin    --- a list of bin # 0 - 255
            out     --- a list of the counts of each bin
    """
    data  = mcf.read_data_file(fname)

    out = [0] * 256

    fdata = []
    for ent in data:
        fdata.append(int(float(ent)))

    avg = numpy.mean(fdata)
    std = numpy.std(fdata)

    for val in fdata:
        if val >= 256:
            continue
        out[val]  += 1

    abin = []
    for k in range(0, 256):
        abin.append(k)

    if str(avg) == 'nan':
        return False

    else:
        return [abin, out, avg, std]

#--------------------------------------------------------------------
#-- check_plot_exist: check data exist and if so whether the plots already exist
#--------------------------------------------------------------------

def check_plot_exist(obsid):
    """
    check data exist and if so whether the plots already exist
    input: obsid    --- obsid
    output: chk     --- 1 if plots needs to be created; 0: skip the plotting
    """
#
#--- first check whether there are data
#
    cfile  = data_dir + 'Fittings/' + str(obsid) + '/pi_center_list'
    if os.path.isfile(cfile):
        sout = os.stat(cfile)
        if int(sout.st_size) == 0:
            return  1
        else:
            chk = 0
    else:
        return 1
#
#--- then check whether the plots are already created
#
    odir = web_dir + 'Plots/Indivisual_Plots/' + str(obsid)

    if os.path.isdir(odir):
        cmd = 'ls ' + odir + '/*.png > ' + zspace
        os.system(cmd)
        out = mcf.read_data_file(zspace, remove=1)
        if len(out) > 0:
            chk = 1

        else:
            chk = 0
    else: 
        chk = 0

    return chk

#--------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) > 1:
        obsid = sys.argv[1]
        obsid = int(float(obsid))

        create_profile_plot(obsid)
    else:
        run_for_all_profile_plot()


