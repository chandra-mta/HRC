#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#       hrc_dose_get_data.py: obtain HRC Evt 1 data for a month and create              #
#                                cumulative data fits files in multiple image files     #
#                                                                                       #
#       author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                       #
#       last updated: Mar 21, 2021                                                      #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
import numpy
import astropy.io.fits  as pyfits
#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param ', shell='tcsh')
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project11/Scripts/house_keeping/dir_list'
f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

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
#
#--- converTimeFormat contains MTA time conversion routines
#
#
#---- mta common functions
#
import mta_common_functions as mcf
#
#--- setting sections for subdeviding image
#
#---for HRC S
#
xstart_s = (   0,    0,     0,     0,     0,     0,     0,     0,     0,     0)
xend_s   = (4095, 4095,  4095,  4095,  4095,  4095,  4095,  4095,  4095,  4095)
ystart_s = (   1, 4916,  9832, 14748, 19664, 24580, 29496, 34412, 39328, 44244)         
yend_s   = (4915, 9831, 14747, 19663, 24579, 29495, 34411, 39327, 44243, 49159)
#
#--- for HRC I
#
xstart_i  = (   1,    1,     1,  5462,  5462,  5462, 10924, 10924, 10924)
xend_i    = (5461, 5461 , 5461, 10923, 10923, 10923, 16385, 16385, 16385)
ystart_i  = (   1, 5462, 10924,     1,  5462, 10924,     1,  5562, 10942)
yend_i    = (5461,10923, 16385,  5461, 10923, 16385,  5461, 10923, 16385)

#---------------------------------------------------------------------------------------
#--- run_hrc_dose: run main hrc dose script by giving data extraction period    
#---------------------------------------------------------------------------------------

def run_hrc_dose(startYear, startMonth, stopYear, stopMonth):
    """
    run main hrc dose script by giving data extraction period
    input:  startYear   --- year of data collection start
            startMonth  --- month of data colleciton start
            stopYear    --- year of data collection stop
            stopMonth   --- month of data coloeciton stop
    output: <data_dir>/Month/Hrc<I/S>_<month>_<year>_<sec>.fits.gz        
            <data_dir>/Cumulative/Hrc<I/S>_08_1999_<month>_<year>_<sec>.fits.gz        
    """
#
#--- start extracting the data for the year/month period
#
    for year in range(startYear, stopYear+1):
#
#---create a list of month appropriate for the year
#
        month_list =  make_month_list(year, startYear, stopYear, startMonth, stopMonth)

        for month in month_list:
#
#--- extract the data for the month/year
#
            hrc_dose_get_data(year, month)

#---------------------------------------------------------------------------------------
#-- hrc_dose_get_data: extract HRC evt1 data from a month and create cumulative data fits file
#---------------------------------------------------------------------------------------

def hrc_dose_get_data(year, month):
    """
    extract HRC evt1 data from a month and create cumulative data fits file. 
    input:  year    --- year of data collection
            month   --- month of data colleciton
    output: <data_dir>/Month/Hrc<I/S>_<month>_<year>_<sec>.fits.gz        
            <data_dir>/Cumulative/Hrc<I/S>_08_1999_<month>_<year>_<sec>.fits.gz        
    """
    smonth = str(month)
    if month < 10:
        smonth = '0' + smonth
#
#--- output file name settings
#
    outfile_i = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    for i in range(0, 9):
        outfile_i[i] = 'HRCI_' + str(smonth) + '_' + str(year) + '_' + str(i) + '.fits'

    outfile_s = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    for i in range(0, 10):
        outfile_s[i] = 'HRCS_' + str(smonth) + '_' + str(year) + '_' + str(i) + '.fits'
#
#--- using ar5gl, get file names
#
    syear = str(year)
    start = syear + '-' + smonth + '-01T00:00:00'

    nextMonth = month + 1
    if nextMonth > 12:
        lyear     = year + 1
        nextMonth = 1
    else:
        lyear  = year                

    smonth = str(nextMonth)
    if nextMonth < 10:
        smonth = '0' + smonth

    syear = str(lyear)
    stop  = str(lyear) + '-' + smonth + '-01T00:00:00'

    line = 'operation=browse\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector=hrc\n'
    line = line + 'level=1\n'
    line = line + 'filetype=evt1\n'
    line = line + 'tstart=' + start + '\n'
    line = line + 'tstop=' +  stop  + '\n'
    line = line + 'go\n'
    f    = open('./zspace', 'w')
    f.write(line)
    f.close()
    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 =  ' /proj/sot/ska/bin/arc5gl -user isobe -script ./zspace > ./zout'
    cmd  = cmd1 + cmd2
    os.system(cmd2)
    mcf.rm_files('./zspace')

    fitsList = mcf.read_data_file('./zout', remove=1)
#
#--- extract each evt1 file, extract the central part, and combine them into a one file
#
#---counters for how many hrc-i and hrc-s are extracted
#
    hrciCnt   = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    hrcsCnt   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  

    chk = 0
    for line in fitsList:
        m = re.search('fits', line)
        if m is None:
            continue
        chk = 1

        try:
            atemp = re.split('\s+', line)
            fitsName  = atemp[0]
            line = 'operation=retrieve\n'
            line = line + 'dataset=flight\n'
            line = line + 'detector=hrc\n'
            line = line + 'level=1\n'
            line = line + 'filetype=evt1\n'
            line = line + 'filename=' + fitsName + '\n'
            line = line + 'go\n'
            f    = open('./zspace', 'w')
            f.write(line)
            f.close()
            cmd1 = "/usr/bin/env PERL5LIB="
            cmd2 =  ' /proj/sot/ska/bin/arc5gl -user isobe -script ./zspace > ./zout'
            cmd  = cmd1 + cmd2
            os.system(cmd2)
            mcf.rm_files('./zspace')

            cmd = 'gzip -d ' + fitsName + '.gz'
            os.system(cmd)
        except:
            continue
#
#---checking which HRC (S or I)
#
        detector = whichHRC(fitsName) 
        print("Processing: " + fitsName + ' : ' + detector)
#
#---- full image is devided into 10 sections
#
        for i in range(0, 10):
#
#--- HRC-I has only 9 sections HRC-S has 10 sections
#
            if detector == 'HRC-I' and i == 9:
                break
#
#--- set command line
#
            line = set_cmd_line(fitsName, detector,  i)
#
#--- create an image file
#
            ichk =  create_image(line, 'ztemp.fits')
            #print("Section: " + str(i) + ' : Status: ' + str(ichk))
#
#--- for HRC S
#
            if (detector == 'HRC-S') and (ichk > 0):
#
#--- add ztemp.fits to fits, if there if no fits, mv ztemp.fits to fits
#
                fits  = 'total_s' + str(i) + '.fits'
                combine_image('ztemp.fits', fits)
                hrcsCnt[i] += 1
#
#--- for HRC I
#
            elif (detector == 'HRC-I') and (ichk > 0):
                fits  = 'total_i' + str(i) + '.fits'
                combine_image('ztemp.fits', fits)
                hrciCnt[i] += 1

        mcf.rm_files('out.fits')
        mcf.rm_files(fitsName)

    for i in range(0, 10):    
        if hrcsCnt[i] > 0:
            cmd = 'mv total_s' + str(i) + '.fits ' + data_s_dir
            cmd = cmd + 'Month/' + outfile_s[i]
            os.system(cmd)
            cmd = 'gzip ' + data_s_dir + 'Month/*.fits'
            os.system(cmd)
        
        createCumulative(year, month, 'HRC-S', data_s_dir, i)

    for i in range(0, 9):    
        if hrciCnt[i] > 0:
            cmd = 'mv total_i' + str(i) + '.fits ' + data_i_dir
            cmd = cmd + 'Month/' + outfile_i[i]
            os.system(cmd)

            cmd = 'gzip ' + data_i_dir + 'Month/*.fits'
            os.system(cmd)
        
        createCumulative(year, month, 'HRC-I', data_i_dir, i)
#
#--- clean up 
#
    cmd = 'rm -rf ./*fits'
    os.system(cmd)

#-----------------------------------------------------------------------------------------
#-- set_cmd_line: generate image creating command line for dmcopy                      ---
#-----------------------------------------------------------------------------------------

def set_cmd_line(fitsName, detector, i):
    """
    generate image creating command line for dmcopy 
    input:  fitsName    --- fits file name
            detector    --- detector name either HRC-S or HRC-I
            i           --- section #
    output: line        --- command line to extract data using dmcopy 
    """
    if detector == 'HRC-S':
        xstart = str(xstart_s[i])
        xend   = str(xend_s[i])
        ystart = str(ystart_s[i])
        yend   = str(yend_s[i])
        mem    = str(200)

    elif detector == 'HRC-I':
        xstart = str(xstart_i[i])
        xend   = str(xend_i[i])
        ystart = str(ystart_i[i])
        yend   = str(yend_i[i])
        mem    = str(200)

    line = fitsName + '[EVENTS][bin rawx='+ xstart + ':' + xend 
    line = line + ':1, rawy=' + ystart + ':' + yend + ':1]'
    line = line + '[status=xxxxxx00xxxxxxxxx000x000xx00xxxx][option type=i4,mem='+ mem + ']'

    return line

#---------------------------------------------------------------------------------
#--- createCumulative: create cumulative hrc data                               --
#---------------------------------------------------------------------------------

def createCumulative(year, month, detector, arch_dir, i=0):
    """
    create cumulative hrc data for a given year and month
    input:  year        --- year
            month       --- month
            detector    --- detector: HRC_I/HRC_S
            arch_dir    --- the directory to save the data
            i           --- section of the image
    output: <data_dir>/Cumulative/Hrc<I/S>_08_1999_<month>_<year>_<sec>.fits.gz        
    """
#
#--- find the previous period
#
    pyear = year
    pmonth = month -1

    if pmonth < 1:
        pmonth = 12
        pyear -= 1

    syear  = str(year)
    smonth = str(month)

    if month < 10:
        smonth = '0' + smonth

    spyear  = str(pyear)
    spmonth = str(pmonth)

    if pmonth < 10:
        spmonth = '0' + spmonth

    if detector == 'HRC-I':
        inst = 'HRCI'
    else:
        inst = 'HRCS'
#
#--- set file names
#
    hrc   = inst + '_'         + smonth  + '_' + syear  + '_' + str(i) + '.fits.gz'
    chrc  = inst + '_08_1999_' + spmonth + '_' + spyear + '_' + str(i) + '.fits.gz'
    chrc2 = inst + '_08_1999_' + smonth  + '_' + syear  + '_' + str(i) + '.fits'
    dfile = arch_dir + 'Month/' + hrc
    pfile = arch_dir + 'Cumulative/' + chrc
    ofile = arch_dir + 'Cumulative/' + chrc2

#
#--- if the monthly file exists, reduce the size of the file before combine it 
#--- into a cumulative data
#
    if os.path.isfile(dfile):
        line = dfile + '[opt type=i2,null=-99]'
        cmd1 = "/usr/bin/env PERL5LIB="
        cmd2 = ' dmcopy infile="' + line + '"  outfile="./ztemp.fits"  clobber="yes"'
        cmd  = cmd1 + cmd2
        bash(cmd,  env=ascdsenv)

        if os.path.isfile(pfile):
            cmd1 = "/usr/bin/env PERL5LIB="
            cmd2 = ' dmimgcalc infile=' + pfile
            cmd2 = cmd2 + ' infile2=ztemp.fits outfile =' 
            cmd2 = cmd2 + ofile + ' operation=add clobber=yes'
            cmd  = cmd1 + cmd2
            bash(cmd,  env=ascdsenv)
            mcf.rm_files('./ztemp.fits')
        else:
            cmd  = 'mv -f ztemp.fits ' + ofile 
            os.system(cmd)


        cmd  = 'gzip ' + ofile 
        os.system(cmd)
#
#--- if the monthly fie does not exist, just copy the last month's cumulative data
#
    else:
        try:
            cmd = 'cp ' + arch_dir + 'Cumulative/' + chrc + ' '  + arch_dir 
            cmd = cmd   + 'Cumulative/'  + chrc2 + '.gz'
            os.system(cmd)
        except:
            pass

#---------------------------------------------------------------------------------
#--- whichHRC: determine HRC I or HRC S observation from HRC event file        ---
#---------------------------------------------------------------------------------

def whichHRC(ifile):
    """
    determine HRC I or HRC S observation from HRC event file; input: HRC event file
    input:  ifile       --- fits file name
    output: detector    --- detector name
    """
    flist     = pyfits.open(ifile)
    try:
        detector  = flist[1].header['DETNAM']
    except:
        detector  = flist[0].header['DETNAM']
    flist.close()

    return detector

#---------------------------------------------------------------------------------
#--- combine_image: combine two fits image files.
#---------------------------------------------------------------------------------

def combine_image(fits1, fits2):
    """
    combine two fits image files. input :fits1 fits2.  a combined fits file 
    will be moved to fits2.
    input:  fits1   --- fits file name1
            fits2   --- fits file name2; the combined fits file will 
                        be saved under this name
    output: fits2   --- updated fits file
    """
    if os.path.isfile(fits2):
        try:
            t1 = pyfits.open(fits1)
            t2 = pyfits.open(fits2)
     
            img1 = t1[0].data
            img2 = t2[0].data
            (x1, y1) = img1.shape
            (x2, y2) = img2.shape
    
            if (x1 != x2) or (y1 != y2):
                chk = 0
     
            else:
                new_img = img1 + img2
                header  = pyfits.getheader(fits1)
                pyfits.writeto('./mtemp.fits', new_img, header)
     
            t1.close()
            t2.close()
    
            mcf.rm_files(fits1)
#
#--- rename the combined fits image to "fits2"
#
            cmd = 'mv -f mtemp.fits ' + fits2
            os.system(cmd)
        except:
            mcf.rm_files(fits1)

    else:
        cmd =  'mv ' + fits1 + ' ' + fits2
        os.system(cmd)

#---------------------------------------------------------------------------------
#--- create_image: create image file according to instruction                  ---
#---------------------------------------------------------------------------------

def create_image(line, outfile):
    """
    create image file according to instruction "line".
    input line: instruction,, outfile: output file name
    """
    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmcopy "' + line + '" out.fits option=image clobber=yes'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)

    try:
        cmd1 = "/usr/bin/env PERL5LIB="
        cmd2 = ' dmstat out.fits centroid=no > stest'
        cmd  = cmd1 + cmd2
        bash(cmd,  env=ascdsenv)
    except:
        pass
#
#--- check actually the image is created
#
    sdata = mcf.read_data_file('./stest', remove=1)

    val = 'NA'
    for lent in sdata:
        m = re.search('mean', lent)
        if m is not None:
            atemp = re.split('\s+|\t+', lent)
            val = atemp[1]
            break
         
    if val != 'NA' and float(val) > 0:
        cmd = 'mv out.fits ' + outfile
        os.system(cmd)

        return 1                        #--- the image file was created
    else:
        return 0                        #--- the image file was not created


#---------------------------------------------------------------------------------
#-- make_month_list: create an appropriate month list for a given conditions  ----
#---------------------------------------------------------------------------------

def make_month_list(year, startYear, stopYear, startMonth, stopMonth):
    """
    create an appropriate month list for a given conditions
    input: year, startYear, stopYear, startMonth, stopMonth
    """
#
#--- fill up the month list
#
    month_list = []

    if startYear == stopYear:
#
#--- the period is in the same year
#
        month_list = range(startMonth, stopMonth+1)
    else:
#
#--- if the period is over two or more years, we need to set three sets of month list
#
        if year == startYear:
            month_list = range(startMonth, 13)
        elif year == stopYear:
            month_list = range(1,stopMonth+1)
        else:
            month_list = range(1,13)

    return month_list

#------------------------------------------------------------------------------------

if __name__ == '__main__':
    
    if len(sys.argv) > 4:
        start_year = int(float(sys.argv[1]))
        start_mon  = int(float(sys.argv[2]))
        stop_year  = int(float(sys.argv[3]))
        stop_mon   = int(float(sys.argv[4]))
        run_hrc_dose(start_year, start_mon, stop_year, stop_mon)

    elif len(sys.argv) > 2:
        year = int(float(sys.argv[1]))
        mon  = int(float(sys.argv[2]))
        hrc_dose_get_data(year, mon)

