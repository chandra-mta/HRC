#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#################################################################################
#                                                                               #
#           plot_ssc_data.py: plot ssc data                                     #
#                                                                               #
#               author: t. isobe (tisobe@cfa.harvard.edu)                       #
#                                                                               #
#               last update: Mar 19, 2021                                       #
#                                                                               #
#################################################################################
import sys
import os
import string
import re
import time
import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines
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

#---------------------------------------------------------------------------------
#-- plot_ssc_data: plot ssc data                                                --
#---------------------------------------------------------------------------------

def plot_ssc_data():
    """
    plot ssc data
    input: none but read from <data_dir>/fifo_data and <data_dir>/temp_data
    output: <html_dir>/fifo_data_plot.png
            <html_dir>/temp_data_plot.png
    """
#
#--- find this year and set the ploting max to the next year
#
    yend = int(time.strftime('%Y', time.gmtime())) + 1

    for fname in ['fifo_data', 'temp_data']:
#
#--- read data
#
        ifile = data_dir + fname
        data  = mcf.read_data_file(ifile)
#
#--- separate data into 3 column data (year, mon, and data)
#
        out   = mcf.separate_data_to_arrays(data)
#
#--- convert the year and month to fractional year
#
        fyear = convert_to_fyear(out[0], out[1])

        ydata = out[2]
        oname = html_dir + fname + '_plot.png'
#
#--- plot data
#
        plot_data(fyear, ydata, 1999, yend, oname)

#---------------------------------------------------------------------------------
#-- convert_to_fyear: convert lists of year and month into a list of fractional year
#---------------------------------------------------------------------------------

def convert_to_fyear(year, mon):
    """
    convert lists of year and month into a list of fractional year
    input:  year    --- a list of year
            mon     --- a list of month
    output: fyear   --- a list of fractional year
    """
    fyear = []
    for k in range(0, len(year)):
        out = year[k] + mon[k] / 12.0
        fyear.append(out)

    return fyear

#---------------------------------------------------------------------------------
#-- plot_data: plot given data set                                              --
#---------------------------------------------------------------------------------

def plot_data(t_list, data, xmin, xmax, outname):
    """
    plot given data set
    input:  t_list  --- x axis data
            data    --- y axis data
            xmin    --- min of x axis
            xmax    --- max of x axis
            outname --- the name of output plot file name
    output: outname
    """
    plt.close('all')

    ymin = 0
    ymax = int(1.2 * max(data))

    ax = plt.subplot(111)
#
#--- make sure that the plotting ranges are not set automatically
#
    ax.set_autoscale_on(False)
    ax.set_xbound(xmin, xmax)
    ax.set_xlim(left=xmin,   right=xmax, auto=False)
    ax.set_ylim(bottom=ymin, top=ymax,   auto=False)

    plt.plot(t_list, data, color='blue', marker='o', markersize='2', lw='1')
    
    plt.ylabel('Occarances per Month')
    plt.xlabel('Time (Year)')

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 6.0)
    plt.savefig(outname, format='png', dpi=200)

    plt.close('all')

#---------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_ssc_data()
