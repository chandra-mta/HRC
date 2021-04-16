#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#           re_process_hrc_data.py: a control script to run reprocess csh scripts               #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Apr 14, 2021                                                           #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import math
import time
import Chandra.Time
import astropy.io.fits  as pyfits
#
#--- from ska
#
from Ska.Shell import getenv, bash
#
#--- set ciao environment 
#
ciaoenv  = getenv('source /soft/ciao/bin/ciao.csh; \
                   source /home/mta/bin/reset_param; setenv PFILES "${PDIRS}";\
                   set path=(/soft/ciao/bin/ $path);', shell='tcsh')
#
#--- reading directory list
#
#path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/house_keeping/dir_list'
path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_I/house_keeping/dir_list'

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

#-----------------------------------------------------------------------------------------
#-- run_process: a control script to run reprocess csh scripts                         ---
#-----------------------------------------------------------------------------------------

def run_process(hrc):
    """
    a control script to run reprocess csh scripts 
    input:  hrc: either "hrc_i" or "hrc_s"
    output: hrc_i_list  --- a list of hrc i obsids which need to be re-processed
            hrc_s_list  --- a list of hrc s obsids which need to be re-processed
            <data_dir>/<obsid>    --- re-processed data direcotry
    """
#
#--- set conditions for either hrc-i or hrc s
#
    if hrc == 'hrc_i':
        out_list = 'hrc_i_list'
        data_dir = '/data/hrc/i/'
        inst     = 'i'
    else:
        out_list = 'hrc_s_list'
        data_dir = '/data/hrc/s/'
        inst     = 's'
#
#--- find un process data 
#
    hlist = find_un_processed_data(inst)

    if hrc == 'hrc_i':
        print("HRC I : " + str(hlist))
    else:
        print("HRC S : " + str(hlist))
    
    for obsid in hlist:
        with open(out_list, 'w') as fo:
            fo.write(str(obsid) + '\n')
#
#--- extract fits data needed for analysis
#
        chk = extract_hrc_data(obsid, data_dir)
        if chk == False:
            print("Not all data are available")
            continue

        if hrc == 'hrc_i':
            cmd = 'csh -f ' + bin_dir + 'repro_all_new.csh hrc_i_list'
        else:
            cmd = 'csh -f ' + bin_dir + 'repro_all_S_new.csh hrc_s_list'

        try:
            run_ciao(cmd)
        except:
            pass

        cdir = data_dir + '/' + str(obsid)
        if os.path.isdir(cdir):
            cmd = 'chgrp -R hat ' + cdir 
            os.system(cmd)
            cmd = 'chmod -R 775 ' + cdir 
            os.system(cmd)
#
#--- directory name should be 5 digit
#
        test = int(float(obsid))
        if test < 10000:
            chk  = mcf.add_leading_zero(obsid, 5)
            odir = data_dir + '/' + str(chk)
            if os.path.isdir(odir):
                cmd = 'mv ' + odir + ' ' + data_dir + '/Duplicate/.'
                os.system(cmd)

            cmd = 'mv ' + cdir + ' ' + odir
            os.system(cmd)

        mcf.rm_files(out_list)
#
#--- correct data file name format
#
        correct_naming(obsid, inst)
#
#--- correct pcad.lst
#
        correct_pacd_path(inst)
#
#--- compress fits files
#
        lobsid = mcf.add_leading_zero(obsid, 5)
        cmd = 'gzip /data/hrc/' + inst + '/' + lobsid + '/*/*fits'
        os.system(cmd)

    chk_proccess_status(inst, hlist)

#-----------------------------------------------------------------------------------------
#-- find_un_processed_data: find hrc obsids which need to be reprocessed                --
#-----------------------------------------------------------------------------------------

def find_un_processed_data(inst):
    """
    find hrc obsids which need to be reprocessed
    input: inst     --- insturment indicator "i" or "s"
    output: uhrc    --- a list of obsids of either hrc i or hrc s
    """
#
#--- extract all hrc obsid listed in database
#
    infile = '/data/mta4/obs_ss/sot_ocat.out'
    data   = mcf.read_data_file(infile)

    h_list = []
    h_dict = {}
    for ent in data:
        atemp = re.split('\^', ent)
        if inst == 'i':
            mc = re.search('HRC-I', atemp[12])
        else:
            mc = re.search('HRC-S', atemp[12])

        if mc is not None:
            atemp = re.split('\^\s+', ent)
            atemp[0].strip()
            try:
                val   = int(float(atemp[1]))
            except:
                continue
            h_list.append(val)
            h_dict[val] = check_status(ent)

        else:
            continue

    uhrc  = clean_the_list(h_list, h_dict, inst)
#
#--- calibration data starting 61* or 62* obsids
#
    clist = find_hrc_calib_obsid(inst)
#
#--- combine them before return 
#
    uhrc  = uhrc + clist
#
#--- check whether the data is already processed: check 5 digit directory
#
    clean = []
    cdir  = '/data/hrc/' + str(inst) + '/'
    for obsid in uhrc:
        chk =  mcf.add_leading_zero(obsid, 5)
        chk = cdir + chk
        if os.path.isdir(chk):
            pass
        else:
            clean.append(obsid)

    return clean

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def check_status(line):

    mc1 = re.search('archived',   line)
    mc2 = re.search('observed',   line)
    mc3 = re.search('scheduled',  line)
    mc4 = re.search('unobserved', line)
    mc5 = re.search('canceled',   line)
    mc6 = re.search('discarded',  line)
    if mc1 is not None:
        return 'archived'
    elif (mc2 is not None) and (mc4 is None):
        return 'observed'
    elif mc3 is not None:
        return 'scheduled'
    elif mc4 is not None:
        return 'unobserved'
    elif mc5 is not None:
        return 'canceled'
    elif mc6 is not None:
        return 'discarded'
    else:
        return 'tbd'

#-----------------------------------------------------------------------------------------
#-- find_processed_data: find the hrc obsids which are already re-processed             --
#-----------------------------------------------------------------------------------------

def find_processed_data(inst):
    """
    find the hrc obsids which are already re-processed
    input:  inst    --- instrument designation: "i" or "s"
    output: out     --- a list of obsids
    """
    if inst == 'i':
        data_dir = '/data/hrc/i/'
    else:
        data_dir = '/data/hrc/s/'

    cmd  = 'ls -d ' + data_dir + '/* > ' + zspace
    os.system(cmd)
    data = mcf.read_data_file(zspace, remove=1)

    out = []
    for ent in data:
        atemp = re.split('\/', ent)
        try:
            val   = int(float(atemp[-1]))
        except:
            continue
        if mcf.is_neumeric(val):
            out.append(val)
#
#--- remove duplicate
#
    oset = set(out)
    out  = list(oset)

    return out

#-----------------------------------------------------------------------------------------
#-- clean_the_list: select out obsids which need re-process                            ---
#-----------------------------------------------------------------------------------------

def clean_the_list(current, cdict,  inst):
    """
    select out obsids which need re-process
    input:  current --- a list of all hrc obsid found in database
            inst    --- instrument designation; "i" or "s"
    output: good    --- a list of obsids to be re-processed
            <house_keeping>/cancelled_list  --- a list of observations canceleld or discarded
    """
#
#--- read the past cancelled obsids
#
    rfile  = house_keeping + 'cancelled_list'
    out    = mcf.read_data_file(rfile)
    cout   = [int(float(x)) for x in out]
    remove = set(cout)
#
#--- find obsids already re-processed
#
    phrc = find_processed_data(inst)
    uhrc = set(current) - set(phrc)
    uhrc = list(uhrc - remove)
#
#--- select out obsids which need to be reprocessed
#
    good = []
    bad  = []
    for obsid in uhrc:
        try:
            status = cdict[obsid]
        except:
            status = 'tbd'

        if status in ['archived', 'observed']:
            good.append(obsid)

        elif status in ['unobserved', 'scheduled']:
            continue

        elif status in ['canceled', 'discarded']:
            bad.append(obsid)

        else:
            continue
#
#--- update cancelled_list if there are new cancelled observations
#
    sbad = set(bad)
    nbad = (set(bad) - remove)

    if len(bad) > 0:
        ncancel = list(remove) + list(nbad)
        fo = open(rfile, 'a')
        for ent in  ncancel:
            fo.write(str(ent))
            fo.write('\n')

        fo.close()

    return good

#-----------------------------------------------------------------------------------------
#-- run_ciao: running ciao comannds                                                    ---
#-----------------------------------------------------------------------------------------

def run_ciao(cmd, clean =0):
    """
    run the command in ciao environment
    input:  cmd --- command line
    clean   --- if 1, it also resets parameters default: 0
    output: command results
    """
    if clean == 1:
        acmd = '/usr/bin/env PERL5LIB=""  source /home/mta/bin/reset_param ;' + cmd
    else:
        acmd = '/usr/bin/env PERL5LIB="" LD_LIBRARY_PATH=""   ' + cmd
    
    try:
        bash(acmd, env=ciaoenv)
    except:
        try:
            bash(acmd, env=ciaoenv)
        except:
            pass

#------------------------------------------------------------------------------------------------
#-- extract_hrc_data: extract fits data needed for analysis                                   ---
#------------------------------------------------------------------------------------------------

def extract_hrc_data(obsid, data_dir):
    """
    extract fits data needed for analysis
    input:  obsid       --- obsid
            data_dir    --- the directory where the data saved
    output: <obsid>/primary/*fits etc
    """
#
#--- extract fits data
#
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'obsid=' + str(obsid) + '\n'
    line = line + 'go\n'

    with open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > zout'
    os.system(cmd)
#
#--- create directories and move the data into them
#
    cmd  = 'mkdir primary secondary'
    os.system(cmd)
 
    cmd  = 'mv *dtf1*fits* *fov*fits* ./primary/.'
    os.system(cmd)

    cmd  = 'mv *bpix1*fits* *evt1*fits* *msk1*fits* *mtl1*fits* \
               *std_dtfstat1.fits* *std_flt1.fits* ./secondary/.'
    os.system(cmd)

    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=pcad\n'
    line = line + 'subdetector=aca\n'
    line = line + 'obsid=' + str(obsid) + '\n'
    line = line + 'go\n'

    with  open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > zout'
    os.system(cmd)
    cmd  = 'mv *asol*fits* ./primary/.'
    os.system(cmd)

    cmd  = 'rm -rf *fits* zline zout'
    os.system(cmd)

    hdir = data_dir + '/' + str(obsid)
    if os.path.isdir(hdir):
        cmd = 'rm -rf ' + hdir + '/*'
        os.system(cmd)
    else:
        cmd = 'mkdir ' + hdir 
        os.system(cmd)

    cmd = 'chmod 774 primary/* secondary/*'
    os.system(cmd)

#
#--- check whether there are duplicated fits files extracted; if so, remove older ones
#
    h_list = ['dtf1', 'fov1', 'asol1']
    sdir   = 'primary'
    remove_duplicate(h_list, sdir)

    h_list = ['bpix1', 'evt1', 'msk1', 'mtl1', 'std_dtfstat1', 'std_flt1']
    sdir   = 'secondary'
    remove_duplicate(h_list, sdir)

    cmd = 'mv primary secondary ' + hdir + '/.'
    os.system(cmd)

    return check_data_exist(hdir)

#------------------------------------------------------------------------------------------------
#-- remove_duplicate: remove duplicated fits files                                             --
#------------------------------------------------------------------------------------------------

def remove_duplicate(h_list, sdir):
    """
    remove duplicated fits files
    input:  h_list  --- a list of head part of the fits files
    sdir--- 'primary' or 'secondary'
    output: cleaned up 'primary' or 'secondary' directory
    """
    for head in h_list:
        cmd = 'ls ./' + sdir + '/*' + head + '* > ' + zspace
        os.system(cmd)
        out  = mcf.read_data_file(zspace, remove=1)
        if len(out) > 1:
            for ent in out[:-1]:
                mcf.rm_files(ent)

#------------------------------------------------------------------------------------------------
#-- check_data_exist: check all needed fits data are extracted                                 --
#------------------------------------------------------------------------------------------------

def check_data_exist(hdir):
    """
    check all needed fits data are extracted
    input:  hdir    --- directory where the data will be kept
    output: True/False. if False, the mail notification is also sent out
    """
    cmd = 'ls ' + hdir + '*/* > ' + zspace
    os.system(cmd)
    with open(zspace, 'r') as f:
        out = f.read()
    mcf.rm_files(zspace)

    for name in ['dtf1', 'fov', 'bpix1', 'evt1', 'msk1', 'mtl1', 'std_dtfstat1', 'std_flt1', 'asol']:
        mc = re.search(name, out)
        if mc is None:
            line = 'Some data files are missing and the re-process is terminated. Check: ' + hdir + '\n'
            with open(zspace, 'w') as fo:
                fo.write(line)

            cmd = 'cat ' + zspace + ' | mailx -s "Subject: hrc re-process failed" tisobe@cfa.harvard.edu'
            os.system(cmd)

            return False

    return True


#------------------------------------------------------------------------------------------------
#-- chk_proccess_status: check whether new processed data actually exist, send out a notification
#------------------------------------------------------------------------------------------------

def chk_proccess_status(inst, hlist):
    """
    check whether new processed data actually exist, and if so send out a notification
    input:  inst    --- instrument: i or s
            hlist   --- a list of obsid
    output: email sent out 
    """

    if inst == 'i':
        data_dir = '/data/hrc/i/'
    else:
        data_dir = '/data/hrc/s/'

    cmd = 'ls ' + data_dir + '* > ' + zspace
    data = mcf.read_data_file(zspace, remove=1)
    d_list = []
    for ent in data:
        if mcf.is_neumeric(ent):
            d_list.append(int(float(ent)))

    done = []
    for obsid in hlist:
        if obsid in d_list:
            done.append(obsid)

    if len(done) > 0:
        line = 'Following obsids are processed for hrc-' + str(inst) + ':\n'
        for obsid in done:
            line = line + '\t' + str(obsid) + '\n'
#
#--- change the status of processed data
#
            cmd = 'chgrp -R hat /data/hrc/i/' + str(obsid)
            os.system(cmd)
            cmd = 'find /data/hrc/i/ -type d -user isobe -exec chmod a+rx,ug+w,o-w {}'
            os.system(cmd)
            cmd = 'chmod -fR a+r,g+w,o-w /data/hrc/i/' + str(obsid)
            os.system(cmd)


        with opne(zspace, 'w') as fo:
            fo.write(line)

        cmd = 'cat ' + zspace + ' |mailx -s "Subject: HRC Obs Re-processed" vkashyap@cfa.harvard.edu'
        os.system(cmd)
        cmd = 'cat ' + zspace + ' |mailx -s "Subject: HRC Obs Re-processed" tisobe@cfa.harvard.edu'
        os.system(cmd)

        mcf.rm_files(zspace)

#------------------------------------------------------------------------------------------------
#-- find_hrc_calib_obsid: find obsids of hrc calibraiton data starting 61* or 62*              --
#------------------------------------------------------------------------------------------------

def find_hrc_calib_obsid(inst):
    """
    find obsids of hrc calibraiton data starting 61* or 62*
    input:  inst    --- s: hrc-s or i: hrc-i
    output: a list of obsids
    """
#
#--- create a list of already processed data
#
    cmd = 'ls -d /data/hrc/' + str(inst) + '/6*  > '+ zspace
    os.system(cmd)
    with open(zspace, 'r') as f:
        ftest = f.read()
    wrd = str(inst) + '/61'
    mc  = re.search(wrd, ftest)
    if mc is not None:
        cmd = 'ls -d /data/hrc/' + str(inst) + '/61* >' + zspace
        os.system(cmd)

    cmd = 'ls -d /data/hrc/' + str(inst) + '/62* >' + zspace
    os.system(cmd)

    data = mcf.read_data_file(zspace, remove=1)
    prev_list = []
    for ent in data:
        atemp = re.split('\/', ent)
        prev_list.append(int(float(atemp[-1])))

#
#--- find today's date and set checking range for the last 30 days
#
    today = time.strftime('%Y:%j:%H:%M:%S', time.gmtime())
    today = int(Chandra.Time.DateTime(today).secs)
    start = today - 10 * 86400
#
#--- extract hrc obsid information
#
    line = 'operation=browse\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'filetype=evt1\n'
    line = line + 'tstart=' + str(start) + '\n'
    line = line + 'tstop='  + str(today) + '\n'
    line = line + 'go\n'

    with open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > ' + zspace
    os.system(cmd)

    mcf.rm_files('./zline')

    data = mcf.read_data_file(zspace, remove=1)
#
#--- select obsids with 61* and 62* starting
#
    h_list = []
    for ent in data:
        mc = re.search('hrcf', ent)
        if mc is not None:
            atemp = re.split('hrcf', ent)
            btemp = re.split('_', atemp[1])
            obsid = int(float(btemp[0]))
            if obsid > 61000 and obsid < 63000:
#
#--- if it is already observed skip it
#
                if obsid in prev_list:
                    continue
#
#--- check which instrument
#
                chk = check_inst(obsid)
                if chk == inst:
                    h_list.append(obsid)

    return h_list

#------------------------------------------------------------------------------------------------
#-- check_inst: find which instremnt this obsid based on                                       --
#------------------------------------------------------------------------------------------------

def check_inst(obsid):
    """
    find which instremnt this obsid based on
    input:  obsid   --- obsid
    output: inst    --- instrument i or s
    """
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'filetype=evt1\n'
    line = line + 'obsid='  + str(obsid) + '\n'
    line = line + 'go\n'

    with open('zline', 'w') as fo:
        fo.write(line)

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline >' + zspace 
    os.system(cmd)

    data = mcf.read_data_file(zspace, remove=1)
    dfile = ''
    for ent in data:
        mc = re.search('hrcf', ent)
        if mc is not None:
            dfile = ent

    if dfile == '':
        return 'na'

    flist = pyfits.open(dfile)
    try:
        inst  = flist[0].header['DETNAM']
    except:
        inst  = flist[1].header['DETNAM']
    flist.close()

    mc   = re.search('HRC-I', inst)
    if mc is not None:
        inst = 'i'
    else:
        inst = 's'

    mcf.rm_file(dfile)

    return inst


#------------------------------------------------------------------------------------------------
#-- correct_naming: check secondary and analysis directories and correct wrongly named fits and par file
#------------------------------------------------------------------------------------------------

def correct_naming(obsid, inst):
    """
    check secondary and analysis directories and correct wrongly named fits and par file
    input:  obsid   --- obsid   
            inst    --- instrument. either "i" or "s"
    """
    cobsid = str(int(float(obsid)))
    if len(cobsid) == 5:
        return 

    lobsid = mcf.add_leading_zero(obsid, 5)
    
    for sdir in ['secondary', 'analysis']:

        cmd = 'ls /data/hrc/' + inst  + '/' + lobsid + '/' + sdir + '/hrcf* >' + zspace
        os.system(cmd)

        data = mcf.read_data_file(zspace, remove=1)
        for ent in data:
            atemp = re.split('\/', ent)
            fname = atemp[-1]
            mc = re.search(lobsid, fname)
            if mc is not None:
                continue
            else:
                atemp = re.split('hrcf', fname)
                btemp = re.split('_',   atemp[1])
                sobs  = btemp[0]
                new   = fname.replace(sobs, lobsid)
                full  = '/data/hrc/' + inst + '/' + lobsid + '/' + sdir + '/' + new

                cmd = 'mv ' + ent + ' ' + full
                os.system(cmd)

#-----------------------------------------------------------------------------------------
#-- correct_pacd_path: correcting pcad.lst data path to a full path                    ---
#-----------------------------------------------------------------------------------------

def correct_pacd_path(obsid, inst):
    """
    correcting pcad.lst data path to a full path
    input:  obsid   --- obsid
            inst    --- instrument either 'i' or 's'
    output: corrected pcad.lst file
    """
    lobsid = mcf.add_leading_zero(obsid, 5)

    cmd = 'ls /data/hrc/' + inst + '/' + lobsid + '/primary/pcad.lst > ' +  zspace
    os.system(cmd)
    data = mcf.read_data_file(zspace, remove=1)

    for ifile in data:
        out = mcf.read_data_file(ifile)
        aline = ''
        for line in out:
            atemp = re.split('\/', line)
            if len(atemp) > 2:
                obsid = atemp[4]
                nform = mcf.add_leading_zero(obsid, 5)
                if obsid != nform:
                    test = line
                    line = line.replace(obsid, nform)
     
                    aline = aline + line + '\n'
            else:
                atemp = re.split('\/', ifile)
                obsid = atemp[4]
                ipath = '/data/hrc/' + inst + '/' + mcf.add_leading_zero(obsid, 5) + '/primary/'
                line  = ipath + line + '\n'
                aline = aline + line


        if aline == '':
            continue

        with open(ifile, 'w') as fo:
            fo.write(aline)

#------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) > 1:
        hrc = sys.argv[1].strip()
    else:
        hrc = 'hrc_i'

    run_process(hrc)
