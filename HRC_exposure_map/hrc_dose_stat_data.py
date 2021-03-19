#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#       hrc_dose_stat_data.py: extract statistics from HRC S and I files                        #
#                                                                                               #
#       author: t. isobe (tisobe@cfa.harvard.edu)                                               #
#                                                                                               #
#       last update: Mar 19, 2021                                                               #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import time
import random
#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release', shell='tcsh')
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project11/Scripts/house_keeping/dir_list'
with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
#
#--- append path to a private folder
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)

import mta_common_functions as mcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())  
zspace = '/tmp/zspace' + str(rtail)

#------------------------------------------------------------------------------------------
#--- comp_stat: compute statistics for the hrc image and print out the result            --
#------------------------------------------------------------------------------------------

def comp_stat(ifile, year, month, ofile):
    """
    compute statistics for the hrc image and print out the result 
    input: ifile    --- hrc image file 
            year    --- year 
            month   --- month, 
            ofile   --- output file name.
    output: ofile   --- a file with stat results updated
    """
    if os.path.isfile(ifile):
#
#--- to avoid getting min value from the outside of the frame edge of a CCD, set threshold
#
        try:
            cmd1 = "/usr/bin/env PERL5LIB="
            cmd2 = ' /bin/nice -n15 dmimgthresh infile=' + ifile 
            cmd2 = cmd2 + ' outfile=zcut.fits  cut="0:1.e10" value=0 clobber=yes'

            cmd  = cmd1 + cmd2
            bash(cmd,  env=ascdsenv)
            cmd1 = "/usr/bin/env PERL5LIB="
            cmd2 = ' dmstat  infile=zcut.fits  centroid=no >' + zspace
            cmd  = cmd1 + cmd2
            bash(cmd,  env=ascdsenv)

            mcf.rm_files('./zcut.fits')
            data = mcf.read_data_file(zspace)
        except:
            data = []
        
        val = 'NA'
        for ent in data:
            ent.lstrip()
            m = re.search('mean', ent)
            if m is not None:
                atemp = re.split('\s+|\t', ent)
                val   = atemp[1]
                break

        if val != 'NA':
#
#--- output of readStat is:
#--- (mean,  dev,  dmin,  dmax , min_pos_x,  min_pos_y,  max_pos_x,  max_pos_y)
#
            out = readStat(zspace)
            mcf.rm_files(zspace)
            (sig1, sig2, sig3) = find_sigma_values(ifile)

        else:
            out = ('NA','NA','NA','NA','NA','NA','NA','NA')
            (sig1, sig2, sig3) = ('NA', 'NA', 'NA')

    else:
        out = ('NA','NA','NA','NA','NA','NA','NA','NA')
        (sig1, sig2, sig3) = ('NA', 'NA', 'NA')
#
#--- print out the results
#
    if out[0] == 'NA':
        line = '%d\t%d\t' % (year, month)
        line =  line + 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n'
    else:
        line = '%d\t%d\t' % (year, month)
        line = line +  '%5.6f\t'   % float(out[0])
        line = line +  '%5.6f\t'   % float(out[1])
        line = line +  '%5.1f\t'   % float(out[2])
        line = line +  '(%d,%d)\t' % (float(out[4]), float(out[5]))
        line = line +  '%5.1f\t'   % float(out[3])
        line = line +  '(%d,%d)\t' % (float(out[6]), float(out[7]))
        line = line +  '%5.1f\t'   % float(sig1)
        line = line +  '%5.1f\t'   % float(sig2)
        line = line +  '%5.1f\n'   % float(sig3)

    if os.path.isfile(ofile):
        with open(ofile, 'a') as fo:
            fo.write(line)
    else:
        with open(ofile, 'w') as fo:
            fo.write(line)

#------------------------------------------------------------------------------------------
#--- readStat:  dmstat output file and extract data values.                             ---
#------------------------------------------------------------------------------------------

def readStat(ifile):

    """
    read dmstat output file and extract data values. 
    input:  ifile    --- input file name
    output: (mean, dev, min, max, min_pos_x, min_pos_y, max_pos_x, max_pos_y)
    """
    mean      = 'NA'
    dev       = 'NA'
    dmin      = 'NA'
    dmax      = 'NA'
    min_pos_x = 'NA'
    min_pos_y = 'NA'
    max_pos_x = 'NA'
    max_pos_y = 'NA'

    data = mcf.read_data_file(ifile)

    for ent in data:
        ent.lstrip()
        atemp = re.split('\s+|\t+', ent)
        m1 = re.search('mean',  ent)
        m2 = re.search('sigma', ent)
        m3 = re.search('min',   ent)
        m4 = re.search('max',   ent)

        if m1 is not None:
            mean = float(atemp[1])

        elif m2 is not None:
            dev  = float(atemp[1])

        elif m3 is not None:
            dmin  = float(atemp[1])
            btemp = re.split('\(', ent)
            ctemp = re.split('\s+|\t+', btemp[1])
            min_pos_x = float(ctemp[1])
            min_pos_y = float(ctemp[2])

        elif m4 is not None:
            dmax  = float(atemp[1])
            btemp = re.split('\(', ent)
            ctemp = re.split('\s+|\t+', btemp[1])
            max_pos_x = float(ctemp[1])
            max_pos_y = float(ctemp[2])

    return (mean, dev, dmin, dmax, min_pos_x, min_pos_y, max_pos_x, max_pos_y)

#------------------------------------------------------------------------------------------
#-- find_sigma_values: find 2 sigma, 3sigma, and 4 sigma values of the given data        --
#------------------------------------------------------------------------------------------
        
def find_sigma_values(fits):
    """
    find 2 sigma, 3sigma, and 4 sigma values of the given data
    input:  fits    --- image fits file name
    output: (sigma1, sigma2, sigma3)
    """
#
#-- make histgram
#
    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmimghist infile=' + fits 
    cmd2 = cmd2 + '  outfile=outfile.fits hist=1::1 strict=yes clobber=yes'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmlist infile=outfile.fits outfile=' + zspace + ' opt=data'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)

    data = mcf.read_data_file(zspace, remove=1)
#
#--- read bin # and its count rate
#
    hbin = []
    hcnt = []
    vsum = 0

    for ent in data:
        atemp = re.split('\s+|\t+', ent)
        if mcf.is_neumeric(atemp[0]):
            hbin.append(float(atemp[1]))
            val = int(atemp[4])
            hcnt.append(val)
            vsum += val
#
#--- checking one sigma and two sigma counts
#
    if len(hbin) > 0:
        v68 = int(0.68  * vsum)
        v95 = int(0.95  * vsum)
        v99 = int(0.997 * vsum)
        sigma1 = -999
        sigma2 = -999
        sigma3 = -999
        acc= 0
        for i in range(0, len(hbin)):
            acc += hcnt[i]
            if acc > v68 and sigma1 < 0:
                sigma1 = hbin[i]

            elif acc > v95 and sigma2 < 0:
                sigma2 = hbin[i]

            elif acc > v99 and sigma3 < 0:
                sigma3 = hbin[i]
                break
    
        return (sigma1, sigma2, sigma3)
    
    else:
        return(0, 0, 0)

#------------------------------------------------------------------------------------------
#--- hrc_dose_extract_stat_data_month: compute HRC statistics                           ---
#------------------------------------------------------------------------------------------

def hrc_dose_extract_stat_data_month(year='NA', month='NA'):
    """
    compute HRC statistics
    input   year    --- year
            month   --- month
    output: <stat_dir>/<inst>_<sec>_acc_out
            <stat_dir>/<inst>_<sec>_dff_out
    """
    if year == 'NA' or month == 'NA':
        year  = raw_input('Year: ')
        month = raw_input('Month: ')
    
    year  = int(year)
    month = int(month)

    syear  = str(year)
    smonth = str(month)
    if month < 10:
        smonth = '0' + smonth

    for i in range(0,10):
        ifile = data_s_dir +  'Cumulative/HRCS_08_1999_' + smonth + '_' 
        ifile = ifile + syear + '_' + str(i) +  '.fits.gz'

        out   = stat_s_dir + '/hrcs_' + str(i) + '_acc_out'
        comp_stat(ifile, year, month, out)

        ifile = data_s_dir  +  'Month/HRCS_' + smonth + '_' 
        ifile = ifile + syear + '_' + str(i) +  '.fits.gz'

        out   = stat_s_dir + '/hrcs_' + str(i) + '_dff_out'
        comp_stat(ifile, year, month, out)

    for i in range(0,9):
        ifile = data_i_dir +  'Cumulative/HRCI_08_1999_' + smonth + '_' 
        ifile = ifile + syear + '_' + str(i) +  '.fits.gz'

        out   = stat_i_dir + '/hrci_' + str(i) + '_acc_out'
        comp_stat(ifile, year, month, out)

        ifile = data_i_dir +  'Month/HRCI_' + smonth + '_' 
        ifile = ifile + syear + '_' + str(i) +  '.fits.gz'

        out   = stat_i_dir + '/hrci_' + str(i) + '_dff_out'
        comp_stat(ifile, year, month, out)

#--------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    year  = int(sys.argv[1])
    month = int(sys.argv[2])

    hrc_dose_extract_stat_data_month(year, month)
