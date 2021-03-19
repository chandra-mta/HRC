#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#       hrc_dose_plot_exposure_stat.py:  plotting trendings of avg, min, max, 1 sigma   #
#                                       2 sigma, and 3 sigma trends                     #
#                                                                                       #
#       author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                       #
#       last update: Mar 19, 2021                                                       #
#                                                                                       #
#########################################################################################

import sys
import os
import string
import re
#
#--- pylab plotting routine related modules
#
import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as lines
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
sys.path.append(mta_dir)
sys.path.append(bin_dir)

import mta_common_functions     as mcf

#------------------------------------------------------------------------------------------------
#--- hrc_dose_plot_exposure_stat: read hrc database, and plot history of exposure             ---
#------------------------------------------------------------------------------------------------

def hrc_dose_plot_exposure_stat():
    """
    read hrc database, and plot history of exposure. 
    input:  none, but read from <stat_i_dir> and <stat_s_dir>
    output: <plt_i_dir>/xxxx.png   <plt_s_dir>/xxxx.png
    """
    for detector in ('hrci', 'hrcs'):
        if detector == 'hrcs':
            indir  = stat_s_dir
            outdir = plt_s_dir
            cstop  = 10
        else:
            indir  = stat_i_dir
            outdir = plt_i_dir
            cstop  = 9

        for section in range(0, cstop):
            inst  = detector + '_' + str(section)
            ofile = outdir + inst + '.png'
#
#--- data list: [date, mean, min, max, std1, std2, std3]
#
            acc_data =  read_hrc_stat_data(indir, inst, 'acc')
            dff_data =  read_hrc_stat_data(indir, inst, 'dff')
#
#--- plot data
#
            try:
                plot_hrc_dose(acc_data, dff_data, ofile)
            except:
                pass

#------------------------------------------------------------------------------------------------
#-- read_hrc_stat_data: read hrc stat data                                                     --
#------------------------------------------------------------------------------------------------

def read_hrc_stat_data(indir, inst, part):
    """
    read hrc stat data
    input:  indir   --- data directory path
            inst    --- hrcs or hrci with section #; e.g.: hrci_0
            part    --- monthly (dff) or cumulative(acc)
    output: date    --- a list of date in fractional year
            avg     --- a list of mean values
            smin    --- a list of min values
            smax    --- a list of max values
            s1      --- a list of 1 sigma values
            s2      --- a list of 2 sigma values
            s3      --- a list of 3 sigma values
    """
    date = []
    avg  = []
    smin = []
    smax = []
    s1   = []
    s2   = []
    s3   = []

    ifile = indir + inst + '_' + part + '_out'
    data  = mcf.read_data_file(ifile)

    for ent in data:
#
#--- if there is no data (monthly data often do not have data), it is listed as NA.
#--- save it 0 instead
#
        mc  = re.search('na', str(ent).lower())
        if mc is not None:
            atemp = re.split('\s+', ent)
            time  = float(atemp[0]) + float(atemp[1])/12.0 + 0.5
            date.append(time)
            avg.append(0)
            smin.append(0)
            smax.append(0)
            s1.append(0)
            s2.append(0)
            s3.append(0)
#
#--- normal case: date, mean, min, max, std1, std2, std3
#
        else:
            atemp = re.split('\s+', ent)
            time  = float(atemp[0]) + float(atemp[1])/12.0 + 0.5
            date.append(time)
            avg.append(float(atemp[2]))
            smin.append(float(atemp[4]))
            smax.append(float(atemp[6]))
            v = float(atemp[8])
            if v < 0.0:
                s1.append(smax[-1])
            else:
                s1.append(v)
            v = float(atemp[9])
            if v < 0.0:
                s2.append(smax[-1])
            else:
                s2.append(v)
            v = float(atemp[10])
            if v < 0.0:
                s3.append(smax[-1])
            else:
                s3.append(v)

    return [date,  avg, smin, smax, s1, s2, s3]

#------------------------------------------------------------------------------------------------
#--- plot_hrc_dose: plot 6 panels of hrc quantities.                                           --
#------------------------------------------------------------------------------------------------

def plot_hrc_dose(acc_data, dff_data, ofile):
    """
    plot 6 trend plots of monthly mean, min, max, and cumulative mean, min, and max
    input : acc_data    --- a list of lists of data of cumulative stats
            dff_data    --- a list of lists of data of monthly stats
            ofile       --- output png file name
    output: ofile       --- a png file
    """
#
#--- open the data lists
#
    [date, amean, amin, amax, accs1, accs2, accs3] = acc_data
    [date, dmean, dmin, dmax, dffs1, dffs2, dffs3] = dff_data
#
#---- set a few parameters
#
    plt.close('all')
    mpl.rcParams['font.size'] = 9
    props = font_manager.FontProperties(size=6)
    plt.subplots_adjust(hspace=0.05)
    plt.subplots_adjust(wspace=0.12)
#
#--- mean
#
    ax1 = plt.subplot(3,2,1)
    plot_panel(date, dmean, 'Average', ax1)
#
#--- mean cumulative
#
    ax2 = plt.subplot(3,2,2)
    plot_panel(date, amean, 'Average Cumulative', ax2)
#
#--- max
#
    ax3 = plt.subplot(3,2,3)
    plot_panel(date, dmax, 'Maximum', ax3)
#
#--- max cumulative
#
    ax4 = plt.subplot(3,2,4)
    plot_panel(date, amax, 'Maximum Cumulative', ax4)
#
#--- 68, 95, and 99.6% levels
#
    labels = ["68% Value ", "95% Value", "99.7% Value"]
    ax5 = plt.subplot(3,2,5)
    plot_panel2(date, dffs1, dffs2, dffs3,  labels, ax5)
#
#--- 68, 95, and 99.6% cumulative
#
    ax6 = plt.subplot(3,2,6)
    plot_panel2(date, accs1, accs2, accs3, labels, ax6)
#
#--- plot x axis tick label only at the bottom ones
#
    for ax in ax1, ax2, ax3, ax4, ax5, ax6:
        if ax != ax5 and ax != ax6:
            for label in ax.get_xticklabels():
                label.set_visible(False)
        else:
            pass
#
#--- putting axis names
#
        ax3.set_ylabel('Counts per Pixel')
        ax5.set_xlabel('Year')
        ax6.set_xlabel('Year')
#
#--- set the size of the plotting area in inch (width: 10.0in, height 10.0in)
#   
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 10.0)
#
#--- save the plot in png format
#   
    plt.savefig(ofile, format='png', dpi=200)

    plt.close('all')

#------------------------------------------------------------------------------------------------
#---   plot_panel: plotting each panel for a given "ax"                                       ---
#------------------------------------------------------------------------------------------------

def plot_panel(x, y, label, ax):

    """
    plotting each panel for a given "ax". 
    input:  x       --- a list of x values
            y       --- a list of y values
            label   --- title label
            ax      ---designation of the plot
    output: a trend plot on pnael <ax>
    """
#
#--- x axis setting: here we assume that x is already sorted
#
    xmin = x[0]
    xmax = x[len(x) -1]
    diff = xmax - xmin
    xmin = int(xmin - 0.05 * diff)
    xmax = int(xmax + 0.05 * diff) + 1
    xbot = int(xmin + 0.05 * diff)
#
#--- y axis setting
#
    ymin = 0
    ymax = max(y)
    try:
        mag  = find_magnitude(ymax)
        mag1 = mag -1
    except:
        mag1 = 0

    ymax = int(ymax /10**mag1 + 2.5) * 10**mag1
    ytop = 0.88 * ymax
#
#--- setting panel 
#

    ax.set_autoscale_on(False)         #---- these three may not be needed for the new pylab, but 
    ax.set_xbound(xmin,xmax)           #---- they are necessary for the older version to set

    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)
#
#--- plot line
#
    plt.plot(x, y, color='blue', lw=1, marker='+', markersize=1.5)

    plt.text(xbot, ytop, label)

#------------------------------------------------------------------------------------------------
#-- plot_panel2: plotting each panel for a given "ax" with three data sets                     --
#------------------------------------------------------------------------------------------------

def plot_panel2(x, s1, s2, s3, labels, ax):

    """
    plotting each panel for a given "ax" with three data sets
    input:  x       --- a list of x values
            s1      --- a list of y values for the first group
            s2      --- a list of y values for the second group
            s3      --- a list of y values for the third group
            label   --- title label
            ax      ---designation of the plot
    output: a trend plot on pnael <ax>
    """
#
#--- x axis setting: here we assume that x is already sorted
#
    xmin = x[0]
    xmax = x[len(x) -1]
    diff = xmax - xmin
    xmin = int(xmin - 0.05 * diff)
    xmax = int(xmax + 0.05 * diff) + 1
    xbot = int(xmin + 0.05 * diff)
#
#--- y axis setting
#
    ymin = 0
    ymax = max(s3)
    try:
        mag  = find_magnitude(ymax)
        mag1 = mag -1
    except:
        mag1 = 0

    ymax = int(ymax /10**mag1 + 2.5) * 10**mag1
    ytop = 0.88 * ymax
#
#--- setting panel 
#
    ax.set_autoscale_on(False)         #---- these three may not be needed for the new pylab, but 
    ax.set_xbound(xmin,xmax)           #---- they are necessary for the older version to set

    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)
#
#--- plot line
#
    p1, = plt.plot(x, s1, color='blue',  lw=1, marker='', markersize=0.0)
    p2, = plt.plot(x, s2, color='green', lw=1, marker='', markersize=0.0)
    p3, = plt.plot(x, s3, color='orange',lw=1, marker='', markersize=0.0)

    legend([p1, p2, p3], [labels[0], labels[1], labels[2]], loc=2, fontsize=9)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

def find_magnitude(val):
    """
    find a magnitude of the value
    input:  val --- numeric value
    output: magnitude
    """

    return math.floor(math.log10(val))

#------------------------------------------------------------------------------------

if __name__ == '__main__':

    hrc_dose_plot_exposure_stat()
