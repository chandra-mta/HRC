#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

import os
import sys
import re
import string


for year in range(1999,2021):
    print("YEAR: " + str(year))
    cmd = './hrc_plotting_maps.py ' + str(year)
    os.system(cmd)
