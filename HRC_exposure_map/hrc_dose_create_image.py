#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#       hrc_dose_create_image.py: convert hrc fits files to png image files                     #
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
import fnmatch 
#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; \
                   source /home/mta/bin/reset_param', shell='tcsh')
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
#--- append path to a privte folder
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- this convert fits files to image files
#
import mta_convert_fits_to_image as mtaimg

#--------------------------------------------------------------------------------------------
#---  create_hrc_maps: create HRC image maps for given year and month                    ----
#--------------------------------------------------------------------------------------------

def create_hrc_maps(year, month, manual=0, chk=0):
    """
    create HRC image maps for given year and month 
    input:  year    --- year
            month   --- month
            manual  --- if > 0, ds9 is used to create images. 
            chk     --- if it is >0, 99.5% cut will be applied for the data in manual mode case
    output: <img_dir>/<Inst>/<Month>/Hrc<inst>_<month>_<year>.png
            <img_dir>/<Inst>/<Cumulative>/Hrc<inst>_08_1999_<month>_<year>.png
    """
#
#--- image for sections for full images
#
    for inst in ['Hrc_I', 'Hrc_S']:
        if inst == 'Hrc_I':
            bdir  = data_i_dir 
            idir  = img_i_dir 
            idir2 = web_img_i_dir
        else:
            bdir  = data_s_dir 
            idir  = img_s_dir 
            idir2 = web_img_s_dir
    
        mdir  = bdir  + 'Month/'
        odir  = idir  + 'Month/'
        odir2 = idir2 + 'Month/'
        if manual == 0:
            hrc_dose_conv_to_png(mdir, odir, year, month)
        else:
            hrc_dose_conv_to_png_manual(mdir, odir, odir2, year, month, chk=0)

        cdir  = bdir  + 'Cumulative/'
        odir  = idir  + 'Cumulative/'
        odir2 = idir2 + 'Cumulative/'
        if manual == 0:
            hrc_dose_conv_to_png(cdir, odir, year, month)
        else:
            hrc_dose_conv_to_png_manual(cdir, odir, odir2,  year, month, chk=0)
    
    
#--------------------------------------------------------------------------------------------
#--- hrc_dose_conv_to_png: convet fits files into png images                              ---
#--------------------------------------------------------------------------------------------

def hrc_dose_conv_to_png(indir, outdir, year, month):
    """
    convet fits files into png images
    input:  indir   --- a directory where to find the data
            outdir  --- image output directory
            yeear   --- year
            month   --- month
    output: <img_dir>/<Inst>/<Month>/Hrc<inst>_<month>_<year>.png
            <img_dir>/<Inst>/<Month>/Hrc<inst>_08_1999_<month>_<year>.png
    """
    syear = str(year)
    smon  = str(month)
    if month < 10:
        smon = '0' + smon

    hname =  'HRC*' + smon + '_' + syear + '*.fits*'

    for ifile in os.listdir(indir):

        if fnmatch.fnmatch(ifile, hname):

            btemp   = re.split('\.fits', ifile)
            out     = btemp[0]
            outfile = outdir + out

            file_p  = indir + ifile

            mtaimg.mta_convert_fits_to_image(file_p, outfile, 'linear', '125x125', 'sls', 'png')
            cmd = 'convert -trim ' + outfile + '.png  ztemp.png'
            os.system(cmd)
            cmd = 'mv ztemp.png ' + outfile + '.png'
            os.system(cmd)
        else:
            pass

#--------------------------------------------------------------------------------------------
#--- hrc_dose_conv_to_png_manual: convet fits files into png images using ds9             ---
#--------------------------------------------------------------------------------------------

def hrc_dose_conv_to_png_manual(indir, outdir, outdir2, year, month, scale='sqrt', color='sls', chk=0):
    """
    convet fits files into png images using ds9. this must be used manually on screen
    input:  indir   --- a directory where to find the data
            outdir  --- image output directory
            outdir2 --- image output directory for html page
            yeear   --- year
            month   --- month
            sclae   --- scale, such sqrt, log, or linear
            color   --- color map name
            chk     --- if it is >0, 99.5% cut will be applied for the data
    output: <img_dir>/<Inst>/<Month>/Hrc<inst>_<month>_<year>.png
            <img_dir>/<Inst>/<Month>/Hrc<inst>_08_1999_<month>_<year>.png
    """
    syear = str(year)
    smon  = str(month)
    if month < 10:
        smon = '0' + smon

    hname =  'HRC*' + smon + '_' + syear + '*.fits*'

    for ifile in os.listdir(indir):

        if fnmatch.fnmatch(ifile, hname):

            btemp    = re.split('\.fits', ifile)
            out      = btemp[0]
            outfile  = outdir  + out + '.png'
            outfile2 = outdir2 + out + '.png'

            ifits   = indir + ifile

            cmd = "/usr/bin/env PERL5LIB= "
            cmd = cmd + ' ds9 ' + ifits + ' -geometry 760x1024 -zoom to fit '
            if chk > 0:
                cmd = cmd + '-scale mode 99.5  -scale ' + scale  +' -cmap ' + color
            else:
                cmd = cmd + '-scale ' + scale  +' -cmap ' + color

            cmd = cmd + ' -colorbar yes -colorbar vertical -colorbar numerics yes -colorbar space value '
            cmd = cmd + ' -colorbar fontsize 12  -saveimage png ' + outfile + ' -exit'

            bash(cmd,  env=ascdsenv)
            cmd = 'cp -f ' + outfile + ' ' + outfile2
            os.system(cmd)
        else:
            pass

#--------------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv)  == 3:
        year  = int(float(sys.argv[1]))
        month = int(float(sys.argv[2]))
        create_hrc_maps(year, month)

    elif len(sys.argv) == 4:
        year   = int(float(sys.argv[1]))
        month  = int(float(sys.argv[2]))
        manual = int(float(sys.argv[3]))
        create_hrc_maps(year, month, manual)

    elif len(sys.argv) == 5:
        year   = int(float(sys.argv[1]))
        month  = int(float(sys.argv[2]))
        manual = int(float(sys.argv[3]))
        chk    = int(float(sys.argv[4]))
        create_hrc_maps(year, month, manual, chk)
    else:
        print("Usage: hrc_dose_create_image.py <year> <month> <manual:optional> <chk:optional>")


