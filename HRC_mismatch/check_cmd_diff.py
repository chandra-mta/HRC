#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#########################################################################################
#                                                                                       #
#   check_cmd_diff.py: compare backstop commands and hrc hk record to find mismatched   #
#                      cases and create info file, if there is the mismatch.            #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Aug 30, 2019                                                   #
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
#-- col name lists
#
hrchk_list = ['time', 'p15cast', 'sptpast', 'n15cast', 'spbpast', 'p24cast', 'imtpast',\
              'imbpast', 'mtrselct', 'mtrcmndr', 'scthast', 'mtritmp', 'mlswenbl', 'fcpuast',\
              'fcpvast', 'cbhuast', 'cbluast', 'cbhvast', 'cblvast', 'wdthast', 'scidpren',\
              'hvpsstat', 's1hvst', 's2hvst', 'imhvlv', 'imhblv', 'sphvlv', 'sphblv', 's1hvlv',\
              's2hvlv', 'uldialv', 'lldialv', 'grdvalv', 'rsrfalv']

cast_list  = ['p15cast', 'p24cast', 'n15cast']
mtrselect  = ['nymtast', 'pymtast', 'clmtast', 'drmtast', 'almtast']
mtrcmndr   = ['mcmdars', 'mcnbamd', 'mcnaamd', 'mclbamd', 'mclaamd', 'mcpramd', 'mdrvast']
mtritmp    = ['drotast', 'droiast']
mlswenbl   = ['sflgast', 'oslsast', 'oplsast', 'cslsast', 'cplsast']
scidpren   = ['clmdast', 'fifoavr', 'obnlasl', 'spmdasl', 'eblkavr', 'cblkavr',\
              'uldiavr', 'wdthavr', 'shldavr']
hvpsstat   = ['sponst', 'spclst', 's1onst', 'imonst', 'imclst', 's2onst']

ifile = house_keeping + 'description'
out = mcf.read_data_file(ifile)
hkall_list = ['time']
for ent in out:
    atemp = re.split('#', ent)
    hkall_list.append(atemp[0].strip())

#----------------------------------------------------------------------------------------
#-- check_cmd_diff: compare backstop commands and hrc hk record to find mismatched cases 
#----------------------------------------------------------------------------------------

def check_cmd_diff(start='', stop=''):
    """
    compare backstop commands and hrc hk record to find mismatched cases
    and create info file containing the mismatch info, if there is the mismatch
    """
#
#--- remove email notificaiton file if it still exists
#
    ifile = exc_dir + 'cmd_diff_cont'
    mcf.rm_files(ifile)
#
#--- if start and stop dates are not given, find today's date and set start and stop time
#
    if start == '':
        today = time.strftime('%Y:%j:00:00:00', time.gmtime())
        stop  = Chandra.Time.DateTime(today).secs
        start = stop - 86400.0 * 2.0                #--- 2 days ago 
#
#--- extract needed information from hrc hk files (extracted from archive)
#
    hk_dict  = read_hk_data(start, stop)
#
#--- read backstop data file for a given period
#
    bs_dict  = read_bs_data(start, stop)
#
#--- compare bk stop data against hrc hk data and send out email 
#
    mismatch1 = compare_two_dicts(bs_dict, hk_dict, offset1= 0, offset2=180)
    out1 = create_notification(mismatch1, pchk=1)

    mismatch2 = compare_two_dicts(hk_dict, bs_dict, offset1= 180, offset2=180)
    out2 = create_notification(mismatch2, pchk=2)

    send_out_notification(out1, out2)

#----------------------------------------------------------------------------------------
#-- read_hk_data: read hrc hk data from fits files and create a data dictionary        --
#----------------------------------------------------------------------------------------

def read_hk_data(start, stop):
    """
    read hrc hk data from fits files and create a data dictionary
    input:  start   --- starting time in seconds from 1998.1.1
            stop    --- stopping time in seconds from 1998.1.1
    output: d_dict  --- data dictionary
    """
    flist = extract_data(start, stop)
    if len(flist) == 0:
        print("No new HK data extracted --- exit")
        exit(1)
#
#--- extracted each data are kept in a dictionary form
#
    d_dict = {}
    for fits in flist:
#
#--- read data from the fits file
#
        hout     = pyfits.open(fits)
        head     = hout[1].header
        fdata    = hout[1].data
        hout.close()
#
#--- select data between start and stop time period
#
        mask     = fdata['TIME'] >= start
        fdata    = fdata[mask]
        mask     = fdata['TIME'] < stop
        fdata    = fdata[mask]
        t_list   = list(fdata['TIME'])
        tlen     = len(t_list)
        if tlen < 1:
            continue
#
#--- header values are repeated for the given data period in the dictionary
#
        tlm_fmt  = head['TLM_FMT']
        datamode = head['DATAMODE']
        tlm_list = [tlm_fmt  for i in range(0, tlen)]
        dmd_list = [datamode for i in range(0, tlen)]
        nul_list = ['x' for i in range(0, tlen)]
        d_dict   = update_dict('tfmt',  d_dict, tlm_list)
        d_dict   = update_dict('dmode', d_dict, tlm_list)
#
#--- update none boolen hrc hk 0 data
#
        val_list = list(numpy.setdiff1d(hrchk_list, cast_list))
        d_dict   = update_d_dict(val_list, fdata, d_dict, tlen)
#
#--- update scidpren part
#
        d_dict   = update_scidpren(scidpren, fdata, d_dict)
#
#--- update cast list
#
        d_dict   = update_cast_list(cast_list, fdata, d_dict)
#
#--- updata mtrcmndr list
#
        d_dict   = update_mtrcmndr_list(mtrcmndr, fdata, d_dict)
#
#---- update mtrselect list
#
        d_dict   = update_stat_list(mtrselect, 'MTRSELCT', 3, fdata, d_dict)
#
#--- update mtritmp list
#
        d_dict   = update_stat_list(mtritmp, 'MTRITMP',  5, fdata, d_dict)
#
#--- update mlswenbl list
#
        d_dict   = update_stat_list(mlswenbl, 'MLSWENBL', 3, fdata, d_dict)
#
#--- update hvpsstat list
#
        d_dict   = update_stat_list(hvpsstat, 'HVPSSTAT', 0, fdata, d_dict)
#
#--- following two data are not checked in this script anymore
#
        try:
            out1 = d_dict['mdrvast']
            out2 = d_dict['fifoavr']
        except:
            out1 = []
            out2 = []
        d_dict['mdrvast'] = out1 + nul_list
        d_dict['fifoavr'] = out2 + nul_list
#
#--- clean up the dictonary entry; remove the entries which are the repeat of the previous
#
    c_dict = {}
    for ent in hkall_list:
        c_dict[ent] = [d_dict[ent][0]]

    t_list = d_dict['time']
    for k in range(1, len(t_list)):
#
#--- check wether any changes from the previous
#
        rchk = 0
        for m in range(1, len(hkall_list)):
            msid = hkall_list[m]
            if c_dict[msid][-1] != d_dict[msid][k]:
                rchk = 1
                break
#
#--- if changes are found, add to a clean dictionary
#
        if rchk == 1:
            for ent in hkall_list:
                try: 
                    out = d_dict[ent][k]
                    c_list = c_dict[ent]
                    c_list.append(out)

                    c_dict[ent] = c_list
                except:
                    pass

    return c_dict

#----------------------------------------------------------------------------------------
#-- update_d_dict: update data dictionary for given column names                       --
#----------------------------------------------------------------------------------------

def update_d_dict(alist, fdata, d_dict, tlen):
    """
    update data dictionary for given column names
    input:  d_dict  --- the current dictionary
            alist   --- a list of fits column names
            fdata   --- fits data array
            tlen    --- the length of list
    output: d_dict  --- updated dictionary
    """
    for ent in alist:
#
#--- check the data exists in the fits file
#
        try:
            dlist = list(fdata[ent])
        except:
            dlist = [-999 for i in range(0, tlen)]
#
#--- check whether the data already exists in the dictionary
#
        d_dict = update_dict(ent, d_dict, dlist)

    return d_dict

#----------------------------------------------------------------------------------------
#-- update_scidpren: update the data dictionary of scidpre part                        --
#----------------------------------------------------------------------------------------

def update_scidpren(col_list, fdata, d_dict):
    """
    update the data dictionary of scidpre part
    input:  col_list    --- a list of names of scidpren data column
            fdata       --- fits data array
            d_dict      --- data dictionary
    output: d_dict      --- updated data dictionary
    """
#
#--- there are two different part of the array that we want to extract data
#
    scol_list = col_list[0:4]
    d_dict    = update_stat_list(scol_list, 'SCIDPREN', 4, fdata, d_dict)

    scol_list = col_list[4:]
    d_dict    = update_stat_list(scol_list, 'SCIDPREN', 7, fdata, d_dict)

    return d_dict

#----------------------------------------------------------------------------------------
#-- update_cast_list: update "cast" entry data update in data dictionary               --
#----------------------------------------------------------------------------------------

def update_cast_list(col_list, fdata, d_dict):
    """
    update "cast" entry data update in data dictionary
    input:  col_list    --- a list of entry names
            fdata       --- fits data array
            d_dict      --- data dictionary
    output: d_dict      --- updated data dictionary
    """
    for col in col_list:
        data   = [row[7] for row in list(fdata[col])]
        d_dict = update_dict(col, d_dict, data)

    return d_dict

#----------------------------------------------------------------------------------------
#-- update_mtrcmndr_list: update mtrcmdr part of the data dictionary                   --
#----------------------------------------------------------------------------------------

def update_mtrcmndr_list(col_list, fdata, d_dict):
    """
    update mtrcmdr part of the data dictionary
    input:  col_list    --- a list of entry names
            fdata       --- fits data array
            d_dict      --- data dictionary
    output: d_dict      --- updated data dictionary
    """
#
#--- the first entry does not have the offset
#
    data   = [row[0] for row in list(fdata['MTRCMNDR'])]
    d_dict = update_dict(col_list[0], d_dict, data)
#
#--- all other has 1 offset
#
    scol_list = col_list[1:]
    d_dict    = update_stat_list(scol_list, 'MTRCMNDR', 1, fdata, d_dict)

    return d_dict

#----------------------------------------------------------------------------------------
#-- update_stat_list: update boolen entries of the dictonary                           --
#----------------------------------------------------------------------------------------

def update_stat_list(col_list, col_nam, offset, fdata, d_dict):
    """
    update boolen entries of the dictonary
    input:  col_list    --- a list of the entry names
            col_name    --- data column name
            offset      --- off set from the entry position in col_list to that of fdata
            fdata       --- fits data array
            d_dict      --- the data dictionary
    output: d_dict      --- updated data dictionary
    """
    for k in range(0, len(col_list)):
        col    = col_list[k]
        m      = k + offset
        data   = [row[m] for row in list(fdata[col_nam])]
        d_dict = update_dict(col, d_dict, data)

    return d_dict

#----------------------------------------------------------------------------------------
#-- update_dict: update dictonary for give col name and a new data list                --
#----------------------------------------------------------------------------------------

def update_dict(cname, d_dict, nlist):
    """
    update dictonary for give col name and a new data list
    input:  cname   --- column name
            d_dict  --- data dictionary
            nlist   --- a list of data
    output: d_idct  --- updated data dictionary
    """
#
#--- check whether the data already exists in the dictionary
#
    try:
        pdata = d_dict[cname]
    except:
        pdata = []
#
#--- update the dictionary
#
    ndata = pdata + nlist
    d_dict[cname] = ndata

    return d_dict

#----------------------------------------------------------------------------------------
#-- extract_data: using arc5gl, extract hk data                                        --
#----------------------------------------------------------------------------------------

def extract_data(start, stop):
    """
    using arc5gl, extract hk data
    input:  start   --- start time (any format accepted by arc5gl)
            stop    --- stop time  (any format accepted by arc5gl)
    output: extracted hrc hk0 fits data files
            flist   --- a list of hk0 fits files
    """
    line  = 'operation=retrieve\n'
    line  = line + 'dataset=flight\n'
    line  = line + 'detector=hrc\n'
    line  = line + 'level=0\n'
    line  = line + 'filetype=hrchk\n'
    line  = line + 'tstart=' + str(start) + '\n'
    line  = line + 'tstop='  + str(stop)  + '\n'
    line  = line + 'go\n'

    flist = mcf.run_arc5gl_process(line)

    return flist

#----------------------------------------------------------------------------------------
#-- read_bs_data: read backstop data                                                   --
#----------------------------------------------------------------------------------------

def read_bs_data(start, stop):
    """
    read backstop data
    input:  start   --- start time in seconds from 1998.1.1
            stop    --- stop time in seconds from 1998.1.1
    output: bs_dict --- dictionary of backstop command data (key: msid + time)
    """ 
#
#--- set values in the intialization value dict
#
    ifile  = house_keeping + 'bk_initialization'
    data   = mcf.read_data_file(ifile)
    i_dict = {}
    for ent in data:
        atemp = re.split('=', ent)
        cname = atemp[0].strip()
        val   = int(float(atemp[1].strip()))
        i_dict[cname] = val
    ptime = 48902399
    i_dict['time'] = 48902399
#
#--- to get a correct intialization values before start time is reached, set the 
#--- data collection period for 30 days before the starting time
#
    pstart = start - 86400 * 30 
#
#--- read backstop data
#
    ifile = data_dir + 'hrc_backstop_extracted2007'
    bdata = mcf.read_data_file(ifile)
#
#--- start working on the main data dict
#
    d_dict = {}
    for ent in hkall_list:
        d_dict[ent] = []

    for ent in bdata:
        [stime, msid, val] = get_bks_value(ent)

        if stime < pstart:
            continue
#
#--- while i_dict prep period, just update the i_dict
#
        if stime < start:
            if stime - ptime > 180:
                i_dict = update_i_dict(msid, val, i_dict)
                i_dict['time'] = stime
                ptime = stime
            continue
#
#--- period ended; stop
#
        elif stime > stop:
            break
#
#--- time is in the data collecting period
#
        for col in hkall_list:
            out = d_dict[col]
            out.append(i_dict[col])
            d_dict[col] = out
#
#--- keep updating i_dict
#
        if stime - ptime > 180:
            i_dict = update_i_dict(msid, val, i_dict)
            i_dict['time'] = stime
            ptime = stime
#
#--- clean up the dictonary entry; remove the entries which are the repeat of the previous
#
    c_dict = {}
    for ent in hkall_list:
        c_dict[ent] = [d_dict[ent][0]]

    t_list = d_dict['time']
    for k in range(1, len(t_list)):
        rchk = 0
        for m in range(1, len(hkall_list)):
            msid = hkall_list[m]
            if c_dict[msid][-1] != d_dict[msid][k]:
                rchk = 1
                break
        if rchk == 1:
            for ent in hkall_list:
                try: 
                    out = d_dict[ent][k]
                    c_list = c_dict[ent]
                    c_list.append(out)

                    c_dict[ent] = c_list
                except:
                    pass
    return c_dict

#----------------------------------------------------------------------------------------
#-- get_bks_value: extract time, msid, and value from the given data line              --
#----------------------------------------------------------------------------------------
        
def get_bks_value(iline):
    """
    extract time, msid, and value from the given data line
    input:  iline   --- data line
                example input line:
                    2019:227:04:07:00.409   682229290   hw  2S2HVOF
                    2019:227:04:07:01.409   682229291   hw  2S1STHV     2S1STHV2=0
    output: stime   --- time in seconds from 1998.1.1
            msid    --- msid
            val     --- value, if it is available
    """
    atemp = re.split('\s+', iline)
    stime = int(float(atemp[1]))

    msid  = atemp[3]
    if len(atemp) > 4:
        mc = re.search('=', atemp[4])
        if mc is not None:
            btemp = re.split('=', atemp[4].strip())
            val   = int(float(btemp[1]))
        else:
            val   = 'na'
    else:
        val = 'na'

    return [stime, msid, val]

#----------------------------------------------------------------------------------------
#-- update_i_dict: update i_dict                                                       --
#----------------------------------------------------------------------------------------

def update_i_dict(msid, val, i_dict):
    """
    update i_dict
    input:  msid    --- msid
            val     --- value for the msid (if exists)
            i_dict  --- data dictionary
    output: i_dict  --- udated data dictionary
    """
    if msid.upper()   == '215PCAOF':    i_dict['p15cast'] = 0     #--- +15V LVPS off 
    elif msid.upper() == '215PCAON':    i_dict['p15cast'] = 1     #--- +15V LVPS on 
    elif msid.upper() == '2SPTTHV':     i_dict['sptpast'] = val   #--- Spect Det Top Plate HV Step
    elif msid.upper() == '215NCAOF':    i_dict['n15cast'] = 0     #--- -15V LVPS off 
    elif msid.upper() == '215NCAON':    i_dict['n15cast'] = 1     #--- -15V LVPS on 
    elif msid.upper() == '2SPTBHV':     i_dict['spbpast'] = val   #--- Spect Det Bottom Plate HV Step
    elif msid.upper() == '224PCAOF':    i_dict['p24cast'] = 0     #--- +24V LVPS off                   
    elif msid.upper() == '224PCAON':    i_dict['p24cast'] = 1     #--- +24V LVPS on                   

    elif msid.upper() == '2IMTTHV':     i_dict['imtpast'] = val   #--- Imaging Det Top Plate HV Step
    elif msid.upper() == '2IMTBHV':     i_dict['imbpast'] = val   #--- Imaging Det Bottom Plate HV Step
#
#---- SCTHAST = 256* 2PSHBALD + 2PSLBALD, except if both values are 0, SCTHAST is unchanged.
#
    elif msid.upper() == '2PSHBALD':     
           i_dict['scthast_p1'] = 256 * val                       #--- Last Motor Step Load Upper Part
    elif msid.upper() == '2PSLBALD':     
           if (i_dict['scthast_p1'] != 0) or (val != 0): 
            i_dict['scthast'] = i_dict['scthast_p1'] +  val       #--- Last Motor Step Load
    elif msid.upper() == '2FCPUALV':    i_dict['fcpuast'] = val   #--- Forced Coarse Position - U
    elif msid.upper() == '2FCPVALV':    i_dict['fcpvast'] = val   #--- Forced Coarse Position - V
    elif msid.upper() == '2CBHUALV':    i_dict['cbhuast'] = val   #--- Blanking High Limit - U
    elif msid.upper() == '2CBLUALV':    i_dict['cbluast'] = val   #--- Blanking Low Limit - U
    elif msid.upper() == '2CBHVALV':    i_dict['cbhvast'] = val   #--- Blanking High Limit - V
    elif msid.upper() == '2CBLVALV':    i_dict['cblvast'] = val   #--- Blanking Low Limit - V
    elif msid.upper() == '2WDTHATH':    i_dict['wdthast'] = val   #--- Width Threshold Setting
    elif msid.upper() == '2S1STHV':     i_dict['s1hvst']  = val   #--- Shield PMT 1 HV Setting
    elif msid.upper() == '2S2STHV':     i_dict['s2hvst']  = val   #--- Shield PMT 2 HV Setting
    elif msid.upper() == '2ULDIATH':    i_dict['uldialv'] = val   #--- Upper Level Disc Setting
#
#--- backstop values and hrc hk values are different in the following sevral commands
#
    elif msid.upper() == '2LLDIATH':                              #--- Lower Level Disc Setting
        i_dict['lldialv'] = val
        if val == 8:    i_dict['lldialv'] = 131
        elif val == 21: i_dict['lldialv'] = 138
        elif val == 48: i_dict['lldialv'] = 152
        elif val == 64: i_dict['lldialv'] = 160
                    
    elif msid.upper() == '2GRDVAAM':                              #--- Grid Bias Setting
        i_dict['grdvalv'] = val
        if val == 28:   i_dict['grdvalv'] = 141
        elif val == 46: i_dict['grdvalv'] = 150
                  
    elif msid.upper() == '2RSRFAAM':                              #--- Range Switch Setting
        i_dict['rsrfalv'] = val
        if val == 115:   i_dict['rsrfalv'] = 185
        elif val == 125: i_dict['rsrfalv'] = 190
#
#--- matches of any values commented out below are not checked.
#
    elif msid.upper() == '2NYMTASL':                              #--- -Y SHUTTER MOTOR SELECT
        i_dict['nymtast'] = 1
        i_dict['almtast'] = 0
    elif msid.upper() == '2PYMTASL':                              #--- +Y SHUTTER MOTOR SELECT
        i_dict['pymtast'] = 1
        i_dict['almtast'] = 0
    elif msid.upper() == '2CLMTASL':                              #--- CALSRC MOTOR SELECT
        i_dict['clmtast'] = 1
        i_dict['almtast'] = 0
    elif msid.upper() == '2DRMTASL':                              #--- DOOR MOTOR SELECT
        i_dict['drmtast'] = 1
        i_dict['almtast'] = 0
    elif msid.upper() == '2ALMTADS':                              #--- ALL MOTORS DESELECT
        i_dict['almtast'] = 1
        i_dict['nymtast'] = 0
        i_dict['pymtast'] = 0
        i_dict['clmtast'] = 0
        i_dict['drmtast'] = 0
    elif msid.upper() == '2NSTBAEX':                              #--- MOVE N STEPS TWRD OPEN'MAX LS
        i_dict['mcnbamd'] = 1
        i_dict['mcnaamd'] = 0
        i_dict['mclbamd'] = 0
        i_dict['mclaamd'] = 0
        i_dict['mcpramd'] = 0    
    elif msid.upper() == '2NSTAAEX':                              #--- MOVE N STEPS TWRD CLOS'HOM LS
        i_dict['mcnaamd'] = 1
        i_dict['mcnbamd'] = 0
        i_dict['mclbamd'] = 0
        i_dict['mclaamd'] = 0
        i_dict['mcpramd'] = 0    
    elif msid.upper() == '2MVLBAEX':                              #--- MOVE TO OPEN'MAX LIMIT SWITCH
        i_dict['mclbamd'] = 1
        i_dict['mcnbamd'] = 0
        i_dict['mcnaamd'] = 0
        i_dict['mclaamd'] = 0
        i_dict['mcpramd'] = 0    
    elif msid.upper() == '2MVLAAEX':                              #--- MOVE TO CLOS'HOME LIM SWITCH
        i_dict['mclaamd'] = 1
        i_dict['mcnbamd'] = 0
        i_dict['mcnaamd'] = 0
        i_dict['mclbamd'] = 0
        i_dict['mcpramd'] = 0    
    elif msid.upper() == '2MVPSAEX':                              #--- STEP FM HOME TO POS CTR VALU  
        i_dict['mcpramd'] = 1
        i_dict['mcnbamd'] = 0
        i_dict['mcnaamd'] = 0
        i_dict['mclbamd'] = 0
        i_dict['mclaamd'] = 0    
    elif msid.upper() == '2MDRVAEN':     i_dict['mdrvast'] = 1    #--- MOTOR DRIVE ENABLE    
    elif msid.upper() == '2MDRVADI':     i_dict['mdrvast'] = 0    #--- MOTOR DRIVE DISABLE    
    elif msid.upper() == '2SMOTAEN':     i_dict['drotast'] = 1    #--- SELECTED MTR OVERTEM PROT ENAB 
    elif msid.upper() == '2SMOTADI':     i_dict['drotast'] = 0    #--- SELECTED MTR OVERTEM PROT DISA 
    elif msid.upper() == '2SMOIAEN':     i_dict['droiast'] = 1    #--- SELECTED MTR OVERCUR PROT ENAB 
    elif msid.upper() == '2SMOIADI':     i_dict['droiast'] = 0    #--- SELECTED MTR OVERCUR PROT DISA 
    elif msid.upper() == '2STFLAEN':     i_dict['sflgast'] = 1    #--- ENABLE STOP FLAGS'CLEAR STOP FLAGS 
    elif msid.upper() == '2STFLADI':     i_dict['sflgast'] = 0    #--- DISABLE STOP FLAGS'CLEAR STOP FLAGS 
    elif msid.upper() == '2OMSLAEN':     i_dict['oslsast'] = 1    #--- OPEN'MAX SECON LIM SW ENAB  
    elif msid.upper() == '2OMSLADI':     i_dict['oslsast'] = 0    #--- OPEN'MAX SECON LIM SW DISA  
    elif msid.upper() == '2OMPLAEN':     i_dict['oplsast'] = 1    #--- OPEN'MAX PRIMARY LIM SW ENAB  
    elif msid.upper() == '2OMPLADI':     i_dict['oplsast'] = 0    #--- OPEN'MAX PRIMARY LIM SW DISA  
    elif msid.upper() == '2CHSLAEN':     i_dict['cslsast'] = 1    #--- CLOS'HOME SECON LIM SW ENAB   
    elif msid.upper() == '2CHSLADI':     i_dict['cslsast'] = 0    #--- CLOS'HOME SECON LIM SW DISA   
    elif msid.upper() == '2CHPLAEN':     i_dict['cplsast'] = 1    #--- CLOS'HOME PRIMARY LIM SW ENAB  
    elif msid.upper() == '2CHPLADI':     i_dict['cplsast'] = 0    #--- CLOS'HOME PRIMARY LIM SW DISA
    elif msid.upper() == '2CLMDAOF':     i_dict['clmdast'] = 0    #--- CALIBRATION MODE DISABLE
    elif msid.upper() == '2CLMDAON':     i_dict['clmdast'] = 1    #--- CALIBRATION MODE ENABLE
    elif msid.upper() == '2FIFOAOF':     i_dict['fifoavr'] = 0    #--- DATA FIFO DISABLE
    elif msid.upper() == '2FIFOAON':     i_dict['fifoavr'] = 1    #--- DATA FIFO ENABLE 
    elif msid.upper() == '2OBSVASL':     i_dict['obnlasl'] = 0    #--- OBSERVING MODE SELECT    
    elif msid.upper() == '2NXILASL':     i_dict['obnlasl'] = 1    #--- NEXT-IN-LINE MODE SELECT    
    elif msid.upper() == '2SPNLASL':     i_dict['spmdasl'] = 1    #--- SPECT DETECTOR SPECT
    elif msid.upper() == '2SPIMASL':     i_dict['spmdasl'] = 0    #--- IMG MODE SELECT   
    elif msid.upper() == '2EBLKADI':     i_dict['eblkavr'] = 0    #--- EBLK VALIDITY DISABLE        
    elif msid.upper() == '2EBLKAEN':     i_dict['eblkavr'] = 1    #--- EBLK VALIDITY ENABLE         
    elif msid.upper() == '2CBLKADI':     i_dict['cblkavr'] = 0    #--- CBLK VALIDITY DISABLE         
    elif msid.upper() == '2CBLKAEN':     i_dict['cblkavr'] = 1    #--- CBLK VALIDITY ENABLE          
    elif msid.upper() == '2ULDIADI':     i_dict['uldiavr'] = 0    #--- ULD VALIDITY DISABLE           
    elif msid.upper() == '2ULDIAEN':     i_dict['uldiavr'] = 1    #--- ULD VALIDITY ENABLE            
    elif msid.upper() == '2WDTHADI':     i_dict['wdthavr'] = 0    #--- WIDTH VALIDITY DISABLE          
    elif msid.upper() == '2WDTHAEN':     i_dict['wdthavr'] = 1    #--- WIDTH VALIDITY ENABLE           
    elif msid.upper() == '2SHL1ADI':     i_dict['shldavr'] = 0    #--- SHIELD VALIDITY DISABLE         
    elif msid.upper() == '2SHL1AEN':     i_dict['shldavr'] = 1    #--- SHIELD VALIDITY ENABLE          
    elif msid.upper() == '2SPHVOF':                               #--- SPECT DET HV OFF
        i_dict['sponst']  = 0
        i_dict['sptpast'] = 0
        i_dict['spbpast'] = 0
    elif msid.upper() == '2SPHVON':     i_dict['sponst']  = 1     #--- SPECT DET HV ON               
    elif msid.upper() == '2SPCLDS':     i_dict['spclst']  = 0     #--- SPECT ILIM DISABLE
    elif msid.upper() == '2SPCLEN':     i_dict['spclst']  = 1     #--- SPECT ILIM ENABLE        
    elif msid.upper() == '2S1HVOF':     i_dict['s1onst']  = 0     #--- SHIELD A HV OFF
    elif msid.upper() == '2S1HVON':     i_dict['s1onst']  = 1     #--- SHIELD A HV ON                
    elif msid.upper() == '2IMHVOF':                               #--- IMAGING DET HV OFF
        i_dict['imonst']  = 0
        i_dict['imtpast'] = 0
        i_dict['imbpast'] = 0
    elif msid.upper() == '2IMHVON':     i_dict['imonst']  = 0     #--- IMAGING DET HV ON             
    elif msid.upper() == '2IMCLDS':     i_dict['imclst']  = 0     #--- IMAGING ILIM DISABLE
    elif msid.upper() == '2IMCLEN':     i_dict['imclst']  = 1     #--- IMAGING ILIM ENABLE       
    elif msid.upper() == '2S2HVOF':     i_dict['s2onst']  = 0     #--- SHIELD B HV OFF
    elif msid.upper() == '2S2HVON':     i_dict['s2onst']  = 1     #--- SHIELD B HV ON                
#
#---- software command related value assignments...
#
    if msid.upper() == 'AODSDITH':                                #---- dither disable
        i_dict['dither']  = 0    
        i_dict['enab91']  = 0
        if i_dict['enab87'] == 1: 
            i_dict['imtpast'] = 42
            i_dict['imbpast'] = 53
            i_dict['imclst']  = 1
        elif i_dict['enab88'] == 1: 
            i_dict['sptpast'] = 43
            i_dict['spbpast'] = 54
            i_dict['spclst']  = 1

    if (msid.upper() == 'AOENDITH') and (i_dict['dither'] == 0): #---- dither enable
        i_dict['dither']  = 1
        if i_dict['enab91'] == 1: 
            if i_dict['enab87'] == 1: 
                i_dict['imtpast'] = 42
                i_dict['imbpast'] = 53
                i_dict['imclst']  = 1
            elif i_dict['enab88'] == 1: 
                i_dict['sptpast'] = 43
                i_dict['spbpast'] = 54
                i_dict['spclst']  = 1
            elif i_dict['enab89'] == 1: 
                i_dict['imbpast'] = 89
                i_dict['imtpast'] = 77
                i_dict['imclst']  = 1
            elif i_dict['enab90'] == 1: 
                i_dict['sptpast'] = 93
                i_dict['spbpast'] = 105
                i_dict['spclst']  = 1
            
    if msid.upper() == 'CODISASX': 
        if val == 87:   i_dict['enab87'] = 0
        elif val == 88: i_dict['enab88'] = 0
        elif val == 89: i_dict['enab89'] = 0
        elif val == 90: i_dict['enab90'] = 0
        elif val == 91: i_dict['enab91'] = 0
        elif val == 92: i_dict['enab92'] = 0
        elif val == 93: i_dict['enab93'] = 0
    
    if msid.upper() == 'COENASX': 
        if val == 87: 
            i_dict['enab87']  = 1
            i_dict['enab87a'] = 1
            i_dict['enab88a'] = 0
        elif val == 88: 
            i_dict['enab88']  = 1
            i_dict['enab87a'] = 0
            i_dict['enab88a'] = 1
        elif val == 89: 
            i_dict['enab89']  = 1
            i_dict['enab89a'] = 1
            i_dict['enab90a'] = 0
        elif val == 90: 
            i_dict['enab90']  = 1
            i_dict['enab89a'] = 0
            i_dict['enab90a'] = 1
        elif val == 91: 
            i_dict['enab91'] = 1
        elif val == 92: 
            i_dict['enab92'] = 1
        elif val == 93: 
            i_dict['enab93'] = 1
        
    if msid.upper() == 'COACTSX': 
        if (i_dict['enab87'] == 1) and (val ==87): 
            i_dict['imtpast'] = 42
            i_dict['imbpast'] = 53
            i_dict['imclst']  = 1
        elif (i_dict['enab88'] == 1) and (val == 88): 
            i_dict['sptpast'] = 43
            i_dict['spbpast'] = 54
            i_dict['spclst']  = 1
        elif (i_dict['enab89'] == 1) and (val == 89): 
            i_dict['imbpast'] = 89
            i_dict['imtpast'] = 77
            i_dict['imclst']  = 1
        elif (i_dict['enab90'] == 1) and (val == 90): 
            i_dict['sptpast'] = 93
            i_dict['spbpast'] = 105
            i_dict['spclst']  = 1
        elif val == 91: 
            if i_dict['dither'] == 1: 
                if i_dict['enab87'] == 1: 
                    i_dict['imtpast'] = 42
                    i_dict['imbpast'] = 53
                    i_dict['imclst']  = 1
                
                if i_dict['enab89'] == 1: 
                    i_dict['imbpast'] = 89
                    i_dict['imtpast'] = 77
                    i_dict['imclst']  = 1
                
                if i_dict['enab88'] == 1: 
                    i_dict['sptpast'] = 43
                    i_dict['spbpast'] = 54
                    i_dict['spclst']  = 1
                
                if i_dict['enab90'] == 1: 
                    i_dict['sptpast'] = 93
                    i_dict['spbpast'] = 105
                    i_dict['spclst']  = 1
            
        elif (i_dict['enab92'] == 1) and (val == 92): 
            i_dict['imtpast'] = 0
            i_dict['imbpast'] = 0
            i_dict['imonst']  = 1
            i_dict['imclst']  = 1
        elif (i_dict['enab93'] == 1) and (val == 93): 
            i_dict['sptpast'] = 0
            i_dict['spbpast'] = 0
            i_dict['sponst']  = 1
            i_dict['spclst']  = 1
            
    return i_dict

#----------------------------------------------------------------------------------------
#-- compare_two_dicts: find a location of command in two dictionaries and compare them --
#----------------------------------------------------------------------------------------

def compare_two_dicts(c1_dict, c2_dict, offset1= 0, offset2=180):
    """
    find a location of command in two dictionaries and compare them
    input:  c1_dict --- dictionary 1
            c2_dict --- dictionary 2
            offset1 --- how match before the ceter time we should include in the test
            offset2 --- how match after the center time we should include in the test
    output: save    --- a list of command mismatch data
    """
    c1_time = c1_dict['time']
    hlen    = len(c1_time)
    c2_time = c2_dict['time']
    blen    = len(c2_time)

    save  = []
    psave = []
    for m in range(0, hlen):
        stime = c1_time[m]
        start = stime - offset1
        stop  = stime + offset2
        chk   = 0
        psave = []
        for k in range(0, blen):
#
#--- time period matches
#
            if c2_time[k] >= start and c2_time[k] < stop:
                psave.append(k)
                chk = 1
            if chk > 0 and stop < c2_time[k]:
                break
#
#--- the case the time of the mian set starts or finishes outside of the compared list
#
        if chk == 0:
            if stime < c2_time[0]:
                psave.append(0)
            elif stime > c2_time[blen-1]:
                psave.append(blen-1)
            else:
                continue

        if len(psave) > 0:
            if m > 0:
#
#--- check which msids are different from the previous round and use only changed msids
#
                chk_list = []
                for ent in hkall_list[1:]:
                    if c1_dict[ent][m-1] != c1_dict[ent][m]:
                        chk_list.append(ent)
            else:
                continue

            if len(chk_list) > 0:
                out = check_cmd_match(c1_dict, c2_dict, chk_list, m, psave[0])
                save = save + out

    save = clean_m_list(save)

    return save 

#----------------------------------------------------------------------------------------
#-- check_cmd_match: check commands mis-match                                          --
#----------------------------------------------------------------------------------------

def check_cmd_match(c1_dict, c2_dict, chk_list,  m, n):
    """
    check commands mis-match
    input:  c1_dict --- dictionary of command set 1
            c2_dict --- dictionary of command set 2
            chk_list    --- msid to be checked
            m       --- position of the command in c1_dict
            n       --- position of the command in c2_dict
    output: save    --- a list of mismatch data in a form of:
                            <time>:<msid>:<base command><compared command>
    """
    time1 = c1_dict['time'][m]
#
#--- check whether the next command change happens less than 180 seconds
#
    save = []
    for msid in chk_list:
        if c1_dict[msid][m] != c2_dict[msid][n]:
            line = str(time1) + ':' + msid + ':' 
            line = line + str(c1_dict[msid][m]) + ':' + str(c2_dict[msid][n])
            save.append(line)

    return save

#----------------------------------------------------------------------------------------
#-- clean_m_list: remove a similar mismatch from the list                              --
#----------------------------------------------------------------------------------------

def clean_m_list(alist):
    """
    remove a similar mismatch from the list
    input:  alist   --- a list of <time>:<msid>:<val1>:<val2>
    output: save    --- a cleaned up list
    """
    if len(alist) < 1:
        return alist
#
#--- save the first list
#
    save  = [alist[0]]
#
#--- start from the second entry
#
    for ent in alist[1:]:
        try:
            atemp = re.split(':', ent)
            stime = float(atemp[0])
            msid  = atemp[1]
            val1  = float(atemp[2])
            val2  = float(atemp[3])
        except:
            continue
#
#--- compare the current line with these already saved in the list
#
        for comp in save:
            atemp = re.split(':', comp)
            ctime = float(atemp[0])
            cid   = atemp[1]
            try:
                cval1 = float(atemp[2])
            except:
                cval1 = atemp[2]
            try:
                cval2 = float(atemp[3])
            except:
                cval2 = atemp[3]

            if msid != cid:
                chk = 1
                continue
            else:
                chk = 0
                diff = abs(ctime - stime)
#
#--- if the value is same, if there is more than 6 hr difference, record it
#
                if diff > 21600:
                    chk = 1
                    continue
                else:
                    if (val1 == cval1) and (val2 == cval2):
                        break
                    else:
                        chk = 1
        if chk == 1:
            save.append(ent)

    return save

#----------------------------------------------------------------------------------------
#-- create_notification: create command mismatch notification content                  --
#----------------------------------------------------------------------------------------

def create_notification(mlist, pchk):
    """
    create command mismatch notification content
    input:  mlist   --- a list of mismatch
            pchk    --- whether this is against backstop or hrc hk commands
    ouput:  line    --- a content of the notificaiton
    """
#
#--- read the description of the msids and create a dict
#
    dfile = house_keeping + 'description'
    out   = mcf.read_data_file(dfile)
    d_dict = {}
    for ent in out:
        atemp = re.split('#', ent)
        msid  = atemp[0].strip().lower()
        d_dict[msid] = ent.upper()

    if pchk == 1:
        line =        '-----------------------------------------------------\n'
        line = line + 'Backstop Command without Matched HRC HK State Change:\n'
        line = line + '-----------------------------------------------------\n'
    else:
        line =        '----------------------------------------------------------\n'
        line = line + 'HRC HK State Change without Corresponding Backstop Command\n'
        line = line + '----------------------------------------------------------\n'

    ptime = 0
    if len(mlist) < 1:
        return ''
#
#--- if the two consequtive "periods" are inside of 3 mins, drop the second set
#--- for a cleaner display of the results
#
    achk = 0
    for ent in mlist:
        atemp = re.split(':', ent)
        stime = float(atemp[0])
        if stime == ptime:
            if achk > 0:
                continue
        elif stime != ptime:
            diff = stime - ptime
            ptime = stime
            if diff < 180:
                achk = 1
                continue
            else:
                achk = 0

            ltime = Chandra.Time.DateTime(stime).date
            btemp = re.split('\.', ltime)
            ltime = btemp[0]
            line  = line + '\nTime: ' + ltime + ' (' + str(atemp[0]) + ')\n'

            if pchk == 1:
                line  = line + 'backstop      hrchk     description\n'
            else:
                line  = line + 'hrchk         backstop  description\n'

        line  = line + adjust_ent(atemp[2]) + '\t\t' + adjust_ent(atemp[3]) + '\t'
        line  = line + d_dict[atemp[1].lower()] + '\n'

    return line

#----------------------------------------------------------------------------------------
#-- send_out_notification: write out command mismatch eamil for email                  --
#----------------------------------------------------------------------------------------

def send_out_notification(line1, line2):
    """
    write out command mismatch eamil for email
    input:  line1   --- mismatch of given backstop commands
            line2   --- mismatch of given hrchk commands
    output: <exc_dir>/cmd_diff_cont
    """
    if line1 != '':
        if line2 != '':
            line = line1 + '\n\n' + line2
        else:
            line = line1
    else:
        line = line2
    
    ifile = exc_dir + 'cmd_diff_cont'
    with open(ifile, 'w') as fo:
        fo.write(line)

#----------------------------------------------------------------------------------------
#-- adjust_ent: adjust the length of entry to 4 letters                                --
#----------------------------------------------------------------------------------------

def adjust_ent(val):
    """
    adjust the length of entry to 4 letters
    input:  val --- the value to be adjucted
    output: val --- the adjusted value
    """
    if val == 'True':
        val = 1
    elif val == 'False':
        val = 0

    val  = str(val).strip()
    vlen = len(val)
    for k in range(vlen, 5):
        val = val + ' ' 

    return val

#----------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 3:
        start = sys.argv[1]
        stop  = sys.argv[2]
#
#--- make sure that in put is in Chandra Time format
#
        mc = re.search(':', start)
        if mc is not None:
            start = Chandra.Time.DateTime(start).secs
            stop  = Chandra.Time.DateTime(stop).secs
        else:
            start = float(start)
            stop  = float(stop)
    else:
        start = ''
        stop  = ''

    check_cmd_diff(start, stop)
