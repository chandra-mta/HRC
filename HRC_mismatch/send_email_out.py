#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#           send_email_out.py: check any mismatches today and, if so, send out email    #
#                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                       #
#           Last Update: Mar 19, 2021                                                   #
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

email  = 'swolk@cfa.harvard.edu,msobolewska@cfa.harvard.edu,dpatnaude@cfa.harvard.edu,tisobe@cfa.harvard.edu'
###email = 'tisobe@cfa.harvard.edu'

#-----------------------------------------------------------------------------------
#-- sent_email_out: check any mismatches today and, if so, send out email         --
#-----------------------------------------------------------------------------------

def sent_email_out():
    """
    check any mismatches today and, if so, send out email
    input: none but read from <exc_dir>/cmd_diff_cont|voltage_cont|dither_cont
    output: email sent out
    """
#
#--- command mis match 
#
    ifile = exc_dir + 'cmd_diff_cont'
    if os.path.isfile(ifile):
        with open(ifile, 'r') as f:
            out = f.read()
        mcf.rm_files(ifile)
    else:
        out = ''
#
#--- voltage mis match
#
    ifile = exc_dir + 'voltage_cont'
    if os.path.isfile(ifile):
        with open(ifile, 'r') as f:
            sout = f.read()
    
        if sout != '':
            if out == '':
                out = sout
            else:
                out = out + '\n\n' + sout
        mcf.rm_files(ifile)
#
#--- dither mis match
#
    ifile = exc_dir + 'dither_cont'
    if os.path.isfile(ifile):
        with open(ifile, 'r') as f:
            sout = f.read()
    
        if sout != '':
            if out == '':
                out = sout
            else:
                out = out + '\n\n' + sout
        mcf.rm_files(ifile)

    if out != '':
        with open(zspace, 'w') as fo:
            fo.write(out)
#
#--- copy the email to the reocrd dir
#
        today = time.strftime('%m_%d_%Y', time.gmtime())
        cmd = 'cp -f ' + zspace + ' ' + record_dir + 'email_' + today
        os.system(cmd)
#
#--- seond out email
#
        cmd = 'cat ' + zspace + '| mailx -s "Subject: HRC Command Mismatch" ' + email
        os.system(cmd)

        mcf.rm_files(zspace)

#-----------------------------------------------------------------------------------

if __name__ == "__main__":

    sent_email_out()
