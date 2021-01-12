#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#####################################################################################################
#                                                                                                   #
#       hrc_gain_fit_voigt.py: extract hrc evt2 files and fit a normal distribution on pha values   #
#                                                                                                   #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                                #
#                                                                                                   #
#           Last Update:    Jan 12, 2021                                                            #
#                                                                                                   #
#####################################################################################################

import os
import sys
import re
import string
import random
import operator
import math
import numpy
import time
import random
import astropy.io.fits  as pyfits
import unittest
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts2/house_keeping/dir_list_py'

with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
#
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import mta_common_functions  as mcf    #---- contains other functions commonly used in MTA scripts
import fit_voigt_profile     as voigt
import adjust_arlac_position as aap
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

working_dir = exc_dir + '/Working_dir/'

#
#--- voigt fitting occasionally issue warning; ignore it
import warnings
warnings.simplefilter("ignore")

#---------------------------------------------------------------------------------------------------
#--- hrc_gain_fit_voigt:  extract hrc evt2 file and create pha distribution                      ---
#---------------------------------------------------------------------------------------------------

def  hrc_gain_fit_voigt(candidate_list):
    """
    extract hrc evt2 file, find the brightest object and create pha distribution
    Input:  candidate_list      --- if it is obsid, use it as a input
                            otherwise, a list of new candidates will be create based database
    Output: <header>_pha.dat    --- pha distribution data
            <header>_gfit.png   --- a plot of pha data
            fitting_results     --- a table of fitted results
    """
#
#--- keep the results in save_line
#
    save_line = ''
#
#--- if an obsid is provided, analyize that, else get new obsids from databases
#
    if len(candidate_list) > 0:
        for obsid in candidate_list:

            print(str(obsid))

            hfile = extract_hrc_evt2(obsid)
            if hfile == 'na':
                continue
#
#--- get a file name header for the later use
#
            temp  = re.split('N', hfile)
            hname = temp[0]
#
#--- extract information from the fits file header
#
            hrdr = pyfits.getheader(hfile, 1)
#
#--- current AR Lac position adjucted with proper motion
#
            date     = hrdr['date-obs']
            [ra, dec] = aap.adjust_arlac_position(date)
#
#--- find the diffrence between real AR Lac position and nominal postion 
#--- so that we can determin how much area we should include 
#

            ra_diff  = abs(ra  - float(hrdr['ra_pnt']))  * 60.0
            dec_diff = abs(dec - float(hrdr['dec_pnt'])) * 60.0
            rad_diff = math.sqrt(ra_diff * ra_diff + dec_diff * dec_diff)

            if rad_diff < 10.0:
                fit_rad = 60.0
            else:
                fit_rad = 200.0
#
#--- find a location of the brightest object (assume it is AR Lac) in sky coordinates
#
            try:
                [x, y] = find_center(hfile)
            except:
                continue
#
#--- extract pha values in the given area
#
            pha = extract_pha(hfile, x, y, fit_rad)

#--- create pha count distribution
#
            pmax = max(pha) + 1
            pha_bin  = [x for x in range(0, pmax)]
            pha_hist = [0 for x in range(0, pmax)]

            for ent in pha:
                pha_hist[ent] += 1
#
#--- print out the distirbution results
#
            line = ''
            for i in range(0, pmax):
                line = line + str(pha_bin[i]) + '\t' + str(pha_hist[i]) + '\n'

            outfile = data_dir + hname + '_pha.dat'
            with  open(outfile, 'w') as fo:
                fo.write(line)
#
#--- gzip file 
#
            cmd = 'gzip -f ' + outfile
            os.system(cmd)
#
#--- find median point
#
            med  = find_med(pha_hist)
#
#--- fit a voigt distribution on the data
#
            try:
                [center, width, amp, alphaD, alphaL, I, a_back, b_back]  \
                    = voigt.fit_voigt_profile(pha_bin, pha_hist, type='voigt', plot_title=hfile)
            except:
                [center, width, amp, alphaD, alphaL, I, a_back, b_back]  \
                    = [0, 0, 0, 0, 0, 0, 0, 0]
#
#--- rename a plotting file
#            
            outfile = plot_dir  + 'Indivisual_Plots/' + hname + '_gfit.png'
            cmd     = 'mv out.png ' + outfile
            os.system(cmd)

            line = str(obsid) + '\t' + str(hrdr['date-obs']) + '\t' + str(hrdr['tstart']) + '\t' 
            line = line + hrdr['detnam']         + '\t'   + str(hrdr['ra_pnt'])   + '\t' 
            line = line + str(hrdr['dec_pnt'])   + '\t\t' + str(round(ra_diff,3)) + '\t' 
            line = line + str(round(dec_diff, 3))+ '\t'   + str(round(rad_diff,3))+ '\t' 
            line = line + str(med)               + '\t\t' + str(round(center, 3)) + '\t' 
            line = line + str(round(amp, 3))     + '\t'   + str(round(width, 3))  + '\t'
            line = line + str(hrdr['roll_pnt'])  + '\t'   + str(hrdr['foc_len'])  + '\t' 
            line = line + str(hrdr['defocus'])   + '\t'   + str(hrdr['sim_x'])    + '\t' 
            line = line + str(hrdr['sim_y'])     + '\t'   + str(hrdr['sim_z'])    + '\t'
            line = line + str(round(alphaD,4))   + '\t'   + str(round(alphaL,4))  + '\t' 
            line = line + str(round(center,3))   + '\t'   + str(round(I,2))       + '\t' 
            line = line + str(round(a_back,2))   + '\t'   + str(round(b_back,2))  + '\n'

            save_line = save_line + line
#
#--- remove the evt2 file
#
            mcf.rm_files(hfile)
#
#--- if there is any new data, print it out
#
    if save_line != '':
#
#--- print out the fitting result
#
        outfile = data_dir + 'fitting_results'

        copied_file = outfile + '~'
        cmd = 'cp ' + outfile + ' ' + copied_file
        os.system(cmd)
    
        with open(outfile, 'a') as fo:
            fo.write(save_line)

#---------------------------------------------------------------------------------------------------
#-- find_med: find median point of pha postion                                                   ---
#---------------------------------------------------------------------------------------------------

def find_med(x):
    """
    find median point of pha postion
    Input:   x ---  a list of pha counts
    OUtput: position of the estimated median
    """
    total = 0
    for ent in x:
        total += ent
    
    half = int(0.5 * total)
    sum =  0
    for i in range(0, len(x)):
        sum += x[i]
        if sum > half:
            pos = i
            break
    
    return i - 1

#---------------------------------------------------------------------------------------------------
#--- extract_hrc_evt2: extract hrc evt2 file                                                     ---
#---------------------------------------------------------------------------------------------------

def extract_hrc_evt2(obsid):
    """
    extract hrc evt2 file 
    Input: obsid    --- obsid of the data
    Output: hrcf<obsid>*evt2.fits.gz
            file name if the data is extracted. if not 'na'
    """
#
#--- write  required arc5gl command
#
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector=hrc\n'
    line = line + 'level=2\n'
    line = line + 'filetype=evt2\n'
    line = line + 'obsid=' + str(obsid) + '\n'
    line = line + 'go\n'

    with  open(zspace, 'w') as fo:
        fo.write(line)

    cmd =  ' /proj/sot/ska/bin/arc5gl  -user isobe -script ' + zspace + ' > fits_list'
#
#--- run arc5gl
#
    os.system(cmd)
    mcf.rm_files(zspace)
#
#--- check the data is actually extracted
#
    data = mcf.read_data_file('fits_list', remove=1)

    for ent in data:
        mc = re.search('fits.gz', ent)
        if mc is not None:
            cmd = 'gzip -d ' + ent 
            os.system(cmd)
            fits = ent.replace('.gz', '')
            return fits
    else:
        return 'na'

#---------------------------------------------------------------------------------------------------
#--- find_center: find the brightest object position from the given event file                   ---
#---------------------------------------------------------------------------------------------------

def find_center(ifile):
    """
    find the brightest object position from the given event file
    Input:      ifile   ---- evnt2 fits file
    Ouput:      xv, yv  ---- sky coordinates of the brightest object
    """
    data = pyfits.getdata(ifile)
    chipx = data.field('X')
    chipy = data.field('Y')
#
#--- because the array is too large to handle in one swipe, divide it into 8x8 segments
#
    xmin   = min(chipx)
    ymin   = min(chipy)
    xmax   = max(chipx)
    ymax   = max(chipy)
    xstep  = int((xmax-xmin) / 8 )
    ystep  = int((ymax-ymin) / 8 )
#
#--- find  the interval which contains largest samples 
#
    cposx = 0
    cposy = 0
    cmax  = 0
    for i in range (0, 8):
        xstart = xstep * i + xmin
        xstop  = xstart + xstep
        for j in range (0, 8):
            ystart = ystep * j + ymin
            ystop  = ystart + ystep

            mask = (data.field('X') >= xstart) & (data.field('X') < xstop) & (data.field('Y') \
                            >= ystart) & (data.field('Y') < ystop)
            temp = data[mask]
            chipx_p = temp.field('X')
            chipy_p = temp.field('Y')
            if len(chipx_p) > cmax:
                cmax = len(chipx_p)
                cposx = i
                cposy = j
#
#--- extract the area of the highest count
#
    xpos_list = []
    ypos_list = []
    maxv_list = []
    xstart = xstep * cposx + xmin
    xstop  = xstart + xstep

    ystart = ystep * cposy + ymin
    ystop  = ystart + ystep

    mask = (data.field('X') >= xstart) & (data.field('X') < xstop) & (data.field('Y') \
                            >= ystart) & (data.field('Y') < ystop)
    temp = data[mask]
    chipx_p = temp.field('X')
    chipy_p = temp.field('Y')
#
#--- count up the events. bin to 2x2 so that we get enough count in each bin
#
    xmin = min(chipx_p)
    xmax = max(chipx_p)
    xdim = int(0.5 * (xmax - xmin)) + 1
    ymin = min(chipy_p)
    ymax = max(chipy_p)
    ydim = int(0.5 * (ymax - ymin)) + 1

    cbin = [[0 for y in range(0, ydim)] for x in range(0, xdim)]
    for j in range(0, len(chipy_p)):
        xpos = int(0.5 * (chipx_p[j]-xmin))
        ypos = int(0.5 * (chipy_p[j]-ymin))
        cbin[xpos][ypos] += 1
#
#--- now find max position
#
    vmax = 0
    xx   = 0
    yy   = 0
    for m in range(0, xdim):
        for n in range(0, ydim):
            if cbin[m][n] > vmax:
                vmax = cbin[m][n]
                xx = m
                yy = n
#
#--- take the mddle of the bin as the brightest spot
#
    xv = int(xx * 2.0 + 1.0  + xmin)
    yv = int(yy * 2.0 + 1.0  + ymin)

    return [xv, yv]

#---------------------------------------------------------------------------------------------------
#--- extract_pha: extract pha data for a given area                                             ----
#---------------------------------------------------------------------------------------------------

def extract_pha(ifile, x, y, drange):
    """
    extract pha data for a given area
    Input:  ifile   ---- input fits file name
            x, y    ---- the center sky coordinates of the area to extract
            range   ---- size of the area +/- x or y
    Outpu:  pha     ---- a list of pha
    """
    xmin = x - drange
    xmax = x + drange
    ymin = y - drange
    ymax = y + drange

    data = pyfits.getdata(ifile)
    mask = (data.field('X') >= xmin) & (data.field('X') < xmax) \
                    & (data.field('Y') >= ymin) & (data.field('Y') < ymax)
    area = data[mask]

    pha  = area.field('PHA')
    pha  = map(int, pha)
    pha  = list(pha)

    return pha

#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):

    """
    testing functions
    """

#------------------------------------------------------------

    def test_find_center(self):

        ifile    = 'hrcf14313N001_evt2.fits'

        position = find_center(ifile)

        test_position = [14459, 14428]
        self.assertEquals(position, test_position)

#------------------------------------------------------------

    def test_extract_pha(self):

        file    = 'hrcf14313N001_evt2.fits'
        [x, y]  = [14459, 14428]
        fit_rad = 60

        pha = extract_pha(file, x, y, fit_rad)
        self.assertEquals(max(pha), 229)

#------------------------------------------------------------

    def test_extract_hrc_evt2(self):

        obsid = '14313'

        ifile  = extract_hrc_evt2(obsid)
        self.assertEquals(ifile, 'hrcf14313N001_evt2.fits')

#--------------------------------------------------------------------

if __name__ == '__main__':

    unittest.main()


