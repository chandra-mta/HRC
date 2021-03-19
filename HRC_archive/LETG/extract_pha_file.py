#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#    extract_pha_file.py: create a pha2 file and a tg directory for LETG observation    #
#                                                                                       #
#               author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                       #
#               last update: Feb 25, 2021                                               #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import math
import time
import astropy.io.fits  as pyfits
#
#--- from ska
#
from Ska.Shell import getenv, bash
#
#--- set ciao environment 
#
ciaoenv = getenv('source /soft/ciao/bin/ciao.csh -o',  shell='tcsh')
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project9/Scripts/LETG/house_keeping/dir_list'

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
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)
#
#--- a couple of other settings
#
acmd = '/usr/bin/env PERL5LIB=""   ;'
ilist = ['i', 's']

#------------------------------------------------------------------------------------
#-- extract_pha_file: create a pha2 file and a tg directory for LETG observation   --
#------------------------------------------------------------------------------------

def extract_pha_file():
    """
    create a pha2 file and a tg directory for LETG observation
    input:  none, but read from /data/hrc/<inst>
    output: /data/hrc/<inst>/<obsid>/analysis/*pha2.fits
            /data/hrc/<inst>/<obsid>/analysis/tg/*
    """
#
#--- get obsid list; d_list = [hrci_list, hrcs_list]
#
    d_list = make_process_list()

    for k in range(0, 2):
        dlen = len(d_list[k])
        for j in range(0, dlen):
            obsid = d_list[k][j]
            if j > 0:
                cmd = 'rm -rf ./' + d_list[k][j-1]
                os.system(cmd)
            print("OBSID: " + str(obsid))
#
#--- original data directory name
#
            d_dir = '/data/hrc/' + ilist[k] + '/' + str(obsid) 
#
#--- copy the data to the local directory
#
            cmd    = 'cp -r ' + d_dir + ' ./' + str(obsid)
            os.system(cmd)
#
#--- remove the files with "new" in the file names
#
            cmd    = 'rm -f ' + str(obsid) + '/secondary/*new*'
            os.system(cmd)
#
#--- extract_pha_file chandra_repro
#
            try:
                pcmd = 'chandra_repro indir=' + str(obsid) + ' outdir=' 
                pcmd = pcmd + str(obsid) + '/new cleanup=no'
                cmd  = acmd + pcmd
                bash(cmd,  env=ciaoenv)
            except:
                continue
#
#--- move pha2 file and tg directroy to the analysis directory
#
            outdir = d_dir + '/analysis/'
            cmd  = 'mkdir -p ' + outdir     #--- just in a case analysis dir does not exist
            os.system(cmd)
            cmd  = 'mv -f '  + str(obsid) + '/new/*_pha2.fits* ' + outdir + '/.'
            os.system(cmd)
    
            cmd  = 'mv  -f ' + str(obsid) + '/new/tg '           + outdir + '/.'
            os.system(cmd)
#
#--- remove the copied data
#
            cmd  = 'rm -rf ./' + str(obsid)
            os.system(cmd)
#
#--- send email
#
#    line = 'HRC pha process finished\n'
#    with open(zspace, 'w') as fo:
#        fo.write(line)
#
#    cmd  = 'cat ' + zspace + '|mailx -s "Subject: HRC PHA finished" tisobe@cfa.harvard.edu'
#    os.system(cmd)
#    cmd  = 'rm -rf ' + zspace
#    os.system(cmd)

#------------------------------------------------------------------------------------
#-- make_process_list: create a list of unprocessed obsid lists                    --
#------------------------------------------------------------------------------------

def make_process_list():
    """
    create a list of unprocessed obsid lists
    input:  none
    output: a list of lists of [<hrc_i obsids>, <hrc_s obsids>]
    """
#
#--- create a dict: obsid <---> grating
#
    [obs_list, dict_inst, dict_grat] = make_inst_dict()
    

    save = []
    for inst in ['i', 's']:
        hdir  = '/data/hrc/' + inst + '/'
        olist = []
#
#--- choose data with evt1 exists in the directory
#
        cmd   = 'ls -d ' + hdir + '*/secondary/*evt1.fits* > ' + zspace + ' 2>/dev/null'
        os.system(cmd)
        out   = mcf.read_data_file(zspace, remove=1)
        for ent in out:
            atemp = re.split('\/', ent)
            obsid = atemp[-3]
#
#--- check whether the pha2 file already exists
#
            cmd   = 'ls ' + hdir + obsid + '/analysis/* > ' + zspace + ' 2>/dev/null'
            os.system(cmd)
            with open(zspace, 'r') as f:
                ochk = f.read()
            cmd   = 'rm -rf ' + zspace
            os.system(cmd)
            mc    = re.search('pha2', ochk)
#
#--- check whether it is an grating observation
#
            if mc is None:
                try:
                    iobsid = str(int(float(obsid)))
                    grat   = dict_grat[iobsid]
                except:
                    grat = 'NONE'
#
#--- special treatment for 6**** level calib observations
#
                if obsid[0] == '6':
                    try:
                        grat = check_grating_from_header(inst, obsid)
                    except: 
                        grat = 'NONE'

                if grat == 'LETG':
                    olist.append(obsid)

        save.append(olist)

    return save

#------------------------------------------------------------------------------------
#-- make_inst_dict: create obsid <---> inst, obsid <---> grating dictionaries       -
#------------------------------------------------------------------------------------

def make_inst_dict():
    """
    create obsid <---> inst, obsid <---> grating dictionaries
    input:  none, but read from /data/mta4/obs_ss/sot_ocat.out
    output: a list of <obsid list>, <dict of instruments>, <dict of grating>
            note: only letg is taken as grating. hetg is ignored
    """
    ifile = '/data/mta4/obs_ss/sot_ocat.out'
    data  = mcf.read_data_file(ifile)

    obs_list  = []
    dict_inst = {}
    dict_grat = {}
    for ent in data:
        mc1 = re.search('HRC', ent)
        if mc1 is None:
            continue
        mc2 = re.search('archived', ent)
        mc3 = re.search('observed', ent)
        if (mc2 is None) and (mc3 is None):
            continue

        mc4  = re.search('LETG', ent)
        if mc4 is not None:
            grat = 'LETG'
        else:
            grat = 'NONE'

        atemp = re.split('\^', ent)
        obsid = atemp[1].strip()
        obsid = str(int(float(obsid)))
        inst  = atemp[12].strip()
        if inst in ['HRC-I', 'HRC-S']:
            obs_list.append(obsid)
            dict_inst[obsid] = inst
            dict_grat[obsid] = grat

    return [obs_list, dict_inst, dict_grat]

#------------------------------------------------------------------------------------
#-- check_grating_from_header: checking grating from a header of the evt1 file of obsid
#------------------------------------------------------------------------------------

def check_grating_from_header(hrc, obsid):
    """
    checking grating from a header of the evt1 file of obsid
    input:  hrc     --- either i or s
            obsid   --- obsid
    output: grat    --- gating, such as LETG, HETG,  or NONE
    """
    cmd  = ' ls /data/hrc/' + hrc + '/' + obsid + '/secondary/*evt1.fits* > ' + zspace + '  2>/dev/null'
    os.system(cmd)
    data = mcf.read_data_file(zspace, remove=1)
    try:
        fits = data[0].strip()
    except:
        return 'NONE'

    flist = pyfits.open(fits)
    try:
        grat  = flist[1].header['GRATING']
    except:
        grat  = 'NONE'

    flist.close()

    return grat

#------------------------------------------------------------------------------------

if __name__ == "__main__":

    extract_pha_file()

