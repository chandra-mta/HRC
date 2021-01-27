#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################################
#                                                                                               #
#   sendout_email.py: send out email to users when new HZ43 observations are added to the data  #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Jan 21, 2021                                                           #
#                                                                                               #
#################################################################################################

import sys
import os
import string
import re
import time
import random
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

with open(path, 'r') as f:
    data = [line.strip() for line in f.readlines()]

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec("%s = %s" %(var, line))
html_top = html_top.replace('#', ':')
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

admin  = 'tisobe@cfa.harvard.edu'

#-----------------------------------------------------------------------------------------
#-- sendout_email: send out email to users when new HZ43 observations are added to the data 
#-----------------------------------------------------------------------------------------

def sendout_email():
    """
    send out email to users when new HZ43 observations are added to the data
    input:  none, but read from <house_keeping>/chk_save
    output: email send out
    """

    cfile = house_keeping + 'chk_save'
    with open(cfile, 'r') as f:
        out = int(f.read())

    if out > 0:
        text = 'New HZ43 observations are added to:\n\n'
        text = text + 'https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/HZ43/hz43.html'
        text = text + '\n'

        with open(zspace, 'w') as fo:
            fo.write(text)

        cmd = 'cat ' + zspace + ' |mailx -s "Subject: New HZ43 Observation" vkashyap@cfa.harvard.edu'
        #####os.system(cmd)

        cmd = 'cat ' + zspace + ' |mailx -s "Subject: New HZ43 Observation" ' + admin
        os.system(cmd)

        mcf.rm_files(zspace)
#
#--- a couple of other things when new observations are in
#
        cmd = 'chgrp -R mtagroup /proj/web-cxc-dmz/htdocs/contrib/cxchrc/HRC_trendings/HZ43/*'
        os.system(cmd)
        cmd = 'chgrp -R mtagroup /data/aschrc6/wilton/isobe/Project8/HZ43/Data/*'
        os.system(cmd)
#
#--- change chk_save to 0
#
        cfile = house_keeping + 'chk_save'
        with open(cfile, 'w') as fo:
            fo.write('0\n')

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    sendout_email()
