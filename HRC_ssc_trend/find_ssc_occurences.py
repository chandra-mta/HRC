#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#############################################################################################
#                                                                                           #
#   find_ssc_occurences.py: find ssc bad data counts for each month and update data files   #
#                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                           #
#           last update: Aug 15, 2019                                                       #
#                                                                                           #
#############################################################################################

import sys
import os
import string
import re
import math
import time
import random

import Ska.engarchive.fetch as fetch
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project10/Scripts/house_keeping/dir_list'

with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var   = atemp[1].strip()
    line  = atemp[0].strip()
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
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

d_list = ['fifo_data', 'temp_data']
masks  = ['0x7f',      '0x0400']

#---------------------------------------------------------------------------------
#-- find_ssc_occurences: find ssc bad data counts for each month and update data file
#---------------------------------------------------------------------------------

def find_ssc_occurences():
    """
    find ssc bad data counts for each month and update data file
    input:  none but read from ska database
    output: <data_dir>/fifo_data <data_dir>/temp_data
    """
    for k in range(0, 2):
        lname = d_list[k]
        ifile = data_dir + lname
        mask  = masks[k]
#
#--- create list of data collection periods in a month long unit
#
        [y_list, m_list] = create_time_list(ifile)
#
#--- read each month data and count bad data
#
        for m in range(0, len(y_list)):
            [start, stop] = set_collection_period(y_list[m], m_list[m])

            print("Processing: " + str(y_list[m]) + ':' + str(m_list[m]))

            cnt  = collect_bad_data(start, stop, mask)
#
#--- print out the result in the data file
#
            update_data_file(y_list[m], m_list[m], cnt, ifile)

#---------------------------------------------------------------------------------
#-- create_time_list: create lists of data collection starting periods in month unit
#---------------------------------------------------------------------------------

def create_time_list(ifile):
    """
    create lists of data collection starting periods in month unit
    input:  ifile   --- data file name
    output: y_list  --- a list of year when each data collection starts
            m_list  --- a list of month when each data collection stats
    """
#
#--- find this month's date
#
    out   = time.strftime('%Y:%m', time.gmtime())
    atemp = re.split(':', out)
    lyear = int(float(atemp[0]))
    lmon  = int(float(atemp[1]))
#
#--- find the last entry date from the data file
#
    [syear, smon] = find_last_entry(ifile)
#
#--- create lists of year and month for each month
#
    y_list = []
    m_list = []
    for year in range(syear, lyear):
        for mon in range(1, 13):
            if year == syear:
                if mon <= smon:
                    continue
            elif year == lyear:
                if mon >= lmon:
                    break

            y_list.append(year)
            m_list.append(mon)

    return[y_list, m_list]

#---------------------------------------------------------------------------------
#-- set_collection_period: set a data collection period in the form of <yyyy>:<jjj>:00:00:00
#---------------------------------------------------------------------------------

def set_collection_period(year, mon):
    """
    set a data collection period in the form of <yyyy>:<jjj>:00:00:00
    input:  year    --- year of the data collection starts
            mon     --- month of the data colleciton starts
    output: start   --- data collection starting time in <yyyy>:<jjj>:00:00:00
            stop    --- data collection stopping time in <yyyy>:<jjj>:00:00:00
    """
#
#--- set finishing date
#
    syear = year
    smon  = mon + 1
    if smon > 12:
        smon = 1
        syear += 1
#
#--- convert time format
#
    start = str(year)  + ':' + mcf.add_leading_zero(mon)  + ':01'
    start = time.strftime('%Y:%j:00:00:00', time.strptime(start, '%Y:%m:%d'))

    stop  = str(syear) + ':' + mcf.add_leading_zero(smon) + ':01'
    stop  = time.strftime('%Y:%j:00:00:00', time.strptime(stop,  '%Y:%m:%d'))

    return [start, stop]

#---------------------------------------------------------------------------------
#-- find_last_entry: find the last entry date                                   --
#---------------------------------------------------------------------------------

def find_last_entry(ifile):
    """
    find the last entry date
    input:  ifile   --- data file name
    output: lyear   --- year of the last entry
            lmon    --- month of the last entry
    """
    out   = mcf.read_data_file(ifile)
    try:
        atemp = re.split('\s+', out[-1])
        lyear = int(float(atemp[0]))
        lmon  = int(float(atemp[1]))
    except:
        lyear = 1999
        lmon  = 8
        
    return [lyear, lmon]

#---------------------------------------------------------------------------------
#-- collect_bad_data: find bad data                                             --
#---------------------------------------------------------------------------------

def collect_bad_data(start, stop, mask):
    """
    find bad data
    input:  start   --- data collection starting time (usually <yyyy>:<jjj>:<hh>:<mm>:<ss>
            stop    --- data collection stopping time
            mask    --- mask to find the data (only two options: '0x7f', '0x0400'
    output: cnt     --- counts of bad data
    """
#
#--- get data from ska database
#
    dat = fetch.Msid('HRC_SS_HK_BAD', start, stop)
#
#--- select bad data
#
    if mask == '0x7f':
        bad  = (dat.vals & 0x7f)   > 0
    else:
        bad  = (dat.vals & 0x0400) > 0
#
#--- count bad data occarnaces
#
    cnt  = count_bad_cases(bad)

    return cnt

#---------------------------------------------------------------------------------
#-- count_bad_cases: count bad data                                             --
#---------------------------------------------------------------------------------

def count_bad_cases(data):
    """
    count bad data
    input:  data    --- a list of data
    output: cnt     --- count of bad data
    """
    cnt  = 0
    for ent in data:
        if ent:
            cnt += 1

    return cnt

#---------------------------------------------------------------------------------
#-- update_data_file: append the new data line to the data file                 --
#---------------------------------------------------------------------------------

def update_data_file(year, mon, cnt, ifile):
    """
    append the new data line to the data file
    input:  year    --- year
            mon     --- month
            cnt     --- count 
            ifile   --- data file name
    output: updated data file
                data format exampe: 2014    6   22 (<year>\t<month>\t<count>)
    """
    line = str(year) + '\t' + str(mon) + '\t' + str(cnt) + '\n'
    with open(ifile, 'a') as fo:
        fo.write(line)

#--------------------------------------------------------------------------------

if __name__ == '__main__':

    find_ssc_occurences()
