#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#####################################################################################
#                                                                                   #
#   create_profile_plot.py: fit gamma profile on the data and creat a plot          #
#                                                                                   #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                #
#                                                                                   #
#           Last Update:    Jul 30, 2019                                            #
#                                                                                   #
#####################################################################################

import os
import sys
import re
import string
import math
import numpy
import time
import random
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
html_top = html_top.replace('#', ':')
#
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import mta_common_functions as mcf    #---- contains other functions commonly used in MTA scripts
import fit_gamma_profile    as gamma
#
#--- temp writing file name
#
rtail  = int(time.time() *random.random())
zspace = '/tmp/zspace' + str(rtail)

#---------------------------------------------------------------------------------------------------
#-- run_for_all_profile_plot: read which obsids are available for profile plotting and plot them   -
#---------------------------------------------------------------------------------------------------

def run_for_all_profile_plot():
    """
    read which obsids are available for profile plotting and plot them
    input: none
    output: <web_dir>/Plots3/Indivisual_Plots/<obsid>/pi_p_list_<#>_vifits.png etc
    """
#
#--- check whether new data came in
#
    ifile = house_keeping + 'hrci_new'
    sfile = house_keeping + 'hrcs_new'
    if (not os.path.isfile(ifile)) and (not os.path.isfile(sfile)):
        exit(1) 

    data1 = mcf.read_data_file(ifile)
    data2 = mcf.read_data_file(sfile)
    olist = data1 + data2
#
#--- get information dictionary
#
    odict = get_obsid_info()

    for ent in olist:
        atemp = re.split('\/', ent)
        obsid = atemp[-1]

        print('ObsID: ' + str(obsid))
        try:
            obs_info = odict[str(obsid)]
        except:
            obs_info = ['', '', '']

        create_profile_plot(obsid, obs_info)

#---------------------------------------------------------------------------------------------------
#-- create_profile_plot: fit gamma profile on the data and creat a plot                          ---
#---------------------------------------------------------------------------------------------------

def  create_profile_plot(obsid, obs_info=''):
    """
    fit gamma profile on the data and creat a plot
    input:  obsid   --- obsid. the data is read from <data_dir>
            obs_info    --- information about the observations
    output: <web_dir>/Plots3/Indivisual_Plots/<obsid>/<head>_<col #>_vfits.png
    """
    if obs_info == '':
        obs_info = get_obsid_info()

    head_list = ['samp_list', 'pi_list' ]

    sdir = data_dir + 'Fittings/' + str(obsid) + '/'

    for m in range(0, len(head_list)):
        cmd = 'mkdir -p ' + sdir
        os.system(cmd)
        fname = sdir + head_list[m]
        line  = create_plot_and_table_entry(head_list[m], obsid,  fname, obs_info)
        if line:
            out  = sdir + head_list[m] +  'fit_results'
            with open(out, 'w') as fo:
                fo.write(line)
        else:
            continue

#---------------------------------------------------------------------------------------------------
#-- create_plot_and_table_entry: fit a gamma profile on the data, create a plot, 
#--                              and return the fitting results 
#---------------------------------------------------------------------------------------------------

def create_plot_and_table_entry(head, obsid, fname, obs_info):
    """
    fit a gamma profile on the data, create a plot, and return the fitting results
    input:  head    --- head part of the file
            obsid   --- obsid
            fname   --- the data file name
    output: <web_dir>/Plots3/Indivisual_Plots/<obsid>_<head>_vfits.png
            line    --- voight fitted results in a table column form
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
    title = 'ObsID: ' + str(obsid) + '(' + obs_info[0] + ': Loc [' 
    title = title + str(obs_info[1]) + ',' + str(obs_info[2]) + '])'

    [a, b, h]  = gamma.fit_gamma_profile(abin, data, avg, std, plot_title=title)
#
#--- rename a plotting file
#            
    outfile = web_dir + 'Plots3/Indivisual_Plots/' + str(obsid) + '_' +  head + '_vfit.png'
    cmd = 'mv out.png ' + outfile
    os.system(cmd)
#
#--- keep the fitting results
#
    line = str(a) + '\t'
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
    data = mcf.read_data_file(fname)

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

#---------------------------------------------------------------------------------------------------
#-- get_obsid_info: create information dictionary  with obsid as a key                            --
#---------------------------------------------------------------------------------------------------

def get_obsid_info():
    """
    create information dictionary  with obsid as a key
    input: none but read from <data_dir>/pi_list_<i/s>
    output: odict   --- dictionary of [inst, yoffset, zoffset] with obsid as a key
    """
    infile = data_dir + 'pi_list_i'
    data1  = mcf.read_data_file(infile)

    infile = data_dir + 'pi_list_s'
    data2  = mcf.read_data_file(infile)

    odict = {}
    for ent in data1:
        atemp = re.split('\s+', ent)
        obsid   = atemp[1]
        yoffset = int(round(float(atemp[10])))
        zoffset = int(round(float(atemp[11])))
        odict[obsid] = ['Hrc I', yoffset, zoffset]

    for ent in data2:
        atemp = re.split('\s+', ent)
        obsid   = atemp[1]
        yoffset = int(round(float(atemp[10])))
        zoffset = int(round(float(atemp[11])))
        odict[obsid] = ['Hrc S', yoffset, zoffset]

    return odict

#--------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) > 1:
        obsid = sys.argv[1]
        obsid.strip()
        if mcf.is_neumeric(obsid):
            obsid = int(float(obsid))
            create_profile_plot(obsid)
        else:
            run_for_all_profile_plot(obsid)
    else:
        run_for_all_profile_plot()


