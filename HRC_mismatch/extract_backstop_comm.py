#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#########################################################################################
#                                                                                       #
#   extract_backstop_comm.py: extract expected HRC HW/SW commands from mplogs           #
#                                                                                       #
#       The relevant SW command MSIDs for SCSs are                                      #
#        COACTSX  - activate SCS, the SCS number is given by COACTS1                    #
#        COENASX  - enable SCS, the SCS number is given by  COENAS1                     #
#        CODISASX - disable SCS, the SCS number is given by CODISAS1                    #
#        COTERMSX - terminate SCS, the SCS number is given by COTERMS1                  #
#                                                                                       #
#       For the HRC SW the relevant SCS numbers are:                                    #
#        87 - HRC-I MCP HV ramp down                                                    #
#        88 - HRC-S MCP HV ramp down                                                    #
#        89 - HRC-I MCP HV ramp up                                                      #
#        90 - HRC-S MCP HV ramp up                                                      #
#        91 - HRC MCP HV Dither Control                                                 #
#        92 - HRC-I MCP HV On                                                           #
#        93 - HRC-S MCP HV On                                                           #
#                                                                                       #
#       The dither MSIDs are:                                                           #
#        AODSDITH - disable dither                                                      #
#        AOENDITH - enable dither                                                       #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Aug 14, 2019                                                   #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import time
import Chandra.Time
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project2/Script3.6/house_keeping/dir_list'
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

mon_list = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]

#----------------------------------------------------------------------------------------
#-- extract_backstop_comm: extract expected HRC HW/SW commands from mplogs             --
#----------------------------------------------------------------------------------------

def extract_backstop_comm():
    """
    extract expected HRC HW/SW commands from mplogs
    input:  none  but read from /data/mpcrit1/mplogs/
    output: updated <data_dir>/hrc_backstop_extracted2007
    """
#
#--- get a backstop data from the last two updated data files
#
    data = get_backstop_data()
#
#--- go through the backstop data file
#
    save  = {}
    tlist = []
    for ent in data:
#
#--- check hw commad or sw command
#
        mc1 = re.search('COMMAND_HW', ent)
        mc2 = re.search('COMMAND_SW', ent)
#
#--- extract HRC Hardware command (TLMSID= 2 indicates HRC)
#
        if mc1 is not None:
            nc1 = re.search('TLMSID= 2',  ent)
            if nc1 is not None:
                atemp = re.split('\s+', ent)
                msid  = atemp[12].replace(',', '')

                nc2 = re.search('SCS', atemp[13])
                if nc2 is None:
                    value = atemp[13]
                else:
                    value = ''
#
#--- convert time to Chandra time
#
                stime = int(Chandra.Time.DateTime(atemp[0]).secs + 0.5)
#
#--- prepend the time to the data line and save
#
                if len(msid) < 8:
                    line = str(atemp[0]) + '\t' + str(stime) + '\thw\t' + msid + '\t\t' + str(value)
                else:
                    line = str(atemp[0]) + '\t' + str(stime) + '\thw\t' + msid + '\t'   + str(value)

                save[stime] = line
                tlist.append(stime)
#
#--- software command cases
#
        elif mc2 is not None:
            nc1 = re.search('COACTSX',  ent)
            nc2 = re.search('COENASX',  ent)
            nc3 = re.search('CODISASX', ent)
            nc4 = re.search('COTERMSX', ent)
            nc5 = re.search('AODSDITH', ent)
            nc6 = re.search('AOENDITH', ent)
            if nc1 is not None:
                [line, stime] = check_scs_val('COACTSX', 'COACTS1', ent)
                if len(line) > 0:
                    save[stime] = line
                    tlist.append(stime)

            elif nc2 is not None:
                [line, stime] = check_scs_val('COENASX', 'COENAS1', ent)
                if len(line) > 0:
                    save[stime] = line
                    tlist.append(stime)

            elif nc3 is not None:
                [line, stime] = check_scs_val('CODISASX', 'CODISAS1', ent)
                if len(line) > 0:
                    save[stime] = line
                    tlist.append(stime)

            elif nc4 is not None:
                [line, stime] = check_scs_val('COTERMSX', 'COTERMS1', ent)
                if len(line) > 0:
                    save[stime] = line
                    tlist.append(stime)

            elif nc5 is  not None:
                atemp = re.split('\s+', ent)
                stime = Chandra.Time.DateTime(atemp[0]).secs
                stime = int(stime + 0.5)
                line  = atemp[0] + '\t' + str(stime) + '\tSW\tAODSDITH\t\t\t<---->'
                save[stime] = line
                tlist.append(stime)

            elif nc6 is  not None:
                atemp = re.split('\s+', ent)
                stime = Chandra.Time.DateTime(atemp[0]).secs
                stime = int(stime + 0.5)
                line  =  atemp[0] + '\t' + str(stime) + '\tSW\tAOENDITH\t\t\t<---->'
                save[stime] = line
                tlist.append(stime)
#
#--- sort the data with time and remove duplicates
#
    tlist = list(set(tlist))
    tlist = sorted(tlist)
#
#--- read the past backstop data and find the last entry time
#
    ifile    = data_dir + 'hrc_backstop_extracted2007'
    backstop = mcf.read_data_file(ifile)
    atemp    = re.split('\s+', backstop[-1])
    le_time  = int(atemp[1])
#
#--- append the result to the data file
#
    line  = ''
    for ent in tlist:

        if ent <= le_time:
            continue
        out  = save[ent]

        line = line + out + '\n'

    if len(line) > 0:
        ofile = data_dir + 'hrc_backstop_extracted2007'

        cmd   = 'cp -f ' + ofile + ' ' + ofile + '~'
        os.system(cmd)

        with open(ofile, 'a') as fo:
            fo.write(line)

#----------------------------------------------------------------------------------------
#-- check_scs_val: find a scs value and create an output data line                     --
#----------------------------------------------------------------------------------------

def check_scs_val(cnam1, cnam2, inline):
    """
    find a scs value and create an output data line
    input:  cnam1   --- command name
            cnam2   --- value name
            inline  --- input line from /data/mpcrit1/mplogs/* data
    output: line    --- an output data line
            stime   --- time stamp of the data in seconds from 1998.1.1
    """
    splt    = cnam2 + '='
    atemp   = re.split(splt, inline)
    btemp   = re.split('\s+', atemp[1])
    scs_val = float(btemp[0])

    if (scs_val >= 87) and (scs_val <= 93):
        ctemp = re.split('\s+',inline)
        stime = Chandra.Time.DateTime(ctemp[0]).secs
        stime = int(stime + 0.5)
        line  = ctemp[0] + '\t' + str(stime) + '\tsw\t'+ cnam1 
        if len(cnam1) < 8:
            line  = line     + '\t\t'+ cnam2 +' =' + str(scs_val)
        else:
            line  = line     + '\t'+ cnam2 +' =' + str(scs_val)

        if scs_val == 91:
            line = line + '\t<----'
    else:
        line  = ''
        stime = 0

    return [line, stime]

#----------------------------------------------------------------------------------------
#-- get_backstop_data: extract backstop data for an appropriate time periods           --
#----------------------------------------------------------------------------------------
    
def get_backstop_data():
    """
    extract backstop data for an appropriate time periods
    input:  none but read from /data/mpcrit1/mplogs/
    output: data    --- a list of data
    """
#
#--- find today's date
#
    out = time.strftime('%Y:%m', time.gmtime())
    atemp = re.split(':', out)
    this_year = int(atemp[0])
    this_mon  = int(atemp[1])
#
#--- check for the last three months of data file names
#
    ym_save   = [[this_year, this_mon]]
    s_mon     = this_mon
    s_year    = this_year
    for k in range(1, 3):
        s_mon -= 1
        if s_mon  < 1:
            s_mon  += 12
            s_year -= 1
        ym_save.append([s_year, s_mon])

    ym_save.reverse()
#
#--- collect data file names
#
    for ym in ym_save:
        [year, dmon] = ym
        smon = dmon - 1
        month = mon_list[smon]
        cmd = 'ls -tr /data/mpcrit1/mplogs/' + str(year) + '/' + month 
        cmd = cmd   + '*/ofls/*.backstop >>' + zspace
        os.system(cmd)

    ldata = mcf.read_data_file(zspace, remove=1)
#
#---- use the last two data files 
#
    data  = mcf.read_data_file(ldata[-2]) + mcf.read_data_file(ldata[-1])

    return data

#----------------------------------------------------------------------------------------

if __name__ == "__main__":

    extract_backstop_comm()
