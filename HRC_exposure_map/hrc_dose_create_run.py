#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################################
#                                                                                               #
#       hrc_dose_create_run.py: a master scripts to run HRC scripts                             #
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
#--- append path to a privte folder
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- this convert fits files to image files
#
import hrc_dose_get_data            as hdgd         #--- extract data
import hrc_dose_create_image        as hdci         #--- create png files
import hrc_dose_stat_data           as hdsd         #--- compute statistics
import hrc_dose_html_updates        as hdhu         #--- update html pages
import hrc_dose_plot_exposure_stat  as hdpes        #--- create stat plot

admin = 'tisobe@cfa.harvard.edu'

#--------------------------------------------------------------------------------------------
#-- hrc_dose_create_run: a master scripts to run HRC scripts                               --
#--------------------------------------------------------------------------------------------

def hrc_dose_create_run(year='', month=''):
    """
    a master scripts to run HRC scripts
    input:  year --- year
            month   --- month
            if they are not given, the script will use the previous month's year/month
    output: HRC data in <data_dir> and png image files in <img_dir>
    """
#
#-- find today's date
#
    if year == '':
        today = time.strftime("%Y:%m:%d", time.gmtime())
        atemp = re.split(':', today)
        year  = int(float(atemp[0]))
        month = int(float(atemp[1]))
        mday  = int(float(atemp[2]))
    else:
        year  = int(float(year))
        month = int(float(month))
        mday  = 30
#
#--- if it is the first 10 days of the month, update for the previous month
#
    if mday >1 and mday < 10:
        month -= 1
        if month < 1:
            month = 12
            year -= 1
#
#--- run scripts
#
    hdgd.hrc_dose_get_data(year, month)
    hdci.hrc_dose_create_image(year, month)
    hdsd.hrc_dose_extract_stat_data_month(year, month)
    hdpes.hrc_dose_plot_exposure_stat()
    hdhu.hrc_dose_make_data_html()
    hdhu.update_main_html()
    hdhu.create_img_html('','')
#
#--- update the date links of the image html page of the one month before the current ones
#
    lyear   = year
    lmonth  = month -1
    if lmonth < 1:
        lmonth = 12
        lyear -= 1
    hdhu.create_img_html(lyear, lmnth)

    cmd = 'rm -rf *fits zout'
    os.system(cmd)
#
#--- send out emial to admin
#
    text = 'HRC Exposure maps for ' + str(year) + '/' + str(month) + ' was prcoessed.\n'
    text = text + "You still need to run:\n"
    text = text + "/data/aschrc6/wilton/isobe/Project11/Scripts/hrc_dose_create_image.py "
    text = text + str(year) + ' ' + str(month) + '1\n'

    with open('./ztemp', 'w') as fo:
        fo.write(text)

    cmd  = 'cat ./ztemp | mailx -s"Subject: HRC Exposure Map Processed" ' + admin
    os.system(cmd)

    os.system('rm -rf .ztemp')

#--------------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv)  == 3:
        year  = int(float(sys.argv[1]))
        month = int(float(sys.argv[2]))
    else:
        year  = ''
        month = ''

    hrc_dose_create_run(year, month)









