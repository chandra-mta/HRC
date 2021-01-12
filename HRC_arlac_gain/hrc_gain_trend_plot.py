#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#########################################################################################
#                                                                                       #
#   hrc_gain_trend_plot.py: create time trend of Voigt profiles fit on HRC PHA data     #
#                                                                                       #
#               author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                       #
#               last update: Nov 05, 2019                                               #
#                                                                                       #
#########################################################################################

import os
import sys
import re
import string
import random
import operator
import time
import numpy
import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as lines

import mpld3
from mpld3 import plugins, utils
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
#--- append a path to a private folder to python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- converTimeFormat contains MTA time conversion routines
#
import mta_common_functions as mcf
#
#--- least sq fitting routine (see https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html)
#
from kapteyn import kmpfit
#
#--- temp writing file name
#
import random
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)
#
#--- some information used in several different locations
#
label_list  = ['Dist < 5',  '5 < Dist < 10', '10< Dist < 15', '15< Dist']
yMinSets    = [0, 0, 0]
yMaxSets    = [250, 250, 100]
yname       = 'PHA'
entLabels   = ["PHA Median", "PHA Voigt Peak Position", "PHA FWHM"]
colorList   = ('blue', 'green', 'red', 'aqua', 'lime', 'fuchsia', 'maroon', 'black', 'yellow', 'olive')
#
#--- html addresses
thtml      = "https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/"
html_trend = 'Trend_plots/'

#---------------------------------------------------------------------------------------------------
#-- hrc_gain_trend_plot: create time trend of Gaussian profiles fit on HRC PHA data              ---
#---------------------------------------------------------------------------------------------------

def hrc_gain_trend_plot():
    """
    create time trend of Gaussian profiles fit on HRC PHA data. It also create trend along 
    radial distnace for each year
    input:  none but the data is read from <house_keeping>/fitting_results
    outut:  time trend plots in <plot_dir> / hrc_i_time_trend.png hrc_s_time_trend.png
            radial_distance plots            hrc_i_radial_dist_year<year>.png
    """
#
#--- time trend plots
#
    xmin  = 1999
    ctemp = int(time.strftime('%Y', time.gmtime()))
    end_year =  int(ctemp) + 1       #--- set end of this year
    xmax  = end_year
    xname = 'Time (Year)'
#
#--- read voigt profile fitting to the pha distribution of each data from the table.
#
    ifile = data_dir + 'fitting_results'
    data  = mcf.read_data_file(ifile)

    obsid_i  = []
    time_i   = []
    dist_i   = []
    med_i    = []
    center_i = []
    width_i  = []
    
    obsid_s  = []
    time_s   = []
    dist_s   = []
    med_s    = []
    center_s = []
    width_s  = []
    
    for ent in data:
        if ent[0] == '#':           #--- avoid comments
            continue

        atemp = re.split('\s+', ent)
#
#--- if the width is < 1 or amp is smaller than 5, probably it did not find a correct source position
#
        if float(atemp[12]) < 1:
            continue

        if float(atemp[11]) < 5:
            continue
#
#--- separate the data to hrc i and hrc s
#
        fyear = mcf.chandratime_to_fraq_year(float(atemp[2]))

        if atemp[3] == 'HRC-I':
            obsid_i.append(atemp[0])
            time_i.append(fyear)
            dist_i.append(float(atemp[8]))
            med_i.append(float(atemp[9]))
            center_i.append(float(atemp[10]))
            width_i.append(float(atemp[12]))
        else:
            obsid_s.append(atemp[0])
            time_s.append(fyear)
            dist_s.append(float(atemp[8]))
            med_s.append(float(atemp[9]))
            center_s.append(float(atemp[10]))
            width_s.append(float(atemp[12]))
#
#--- create an interactive plot for HRC I
#
    xSets    = [time_i, time_i, time_i]
    ySets    = [med_i, center_i, width_i]
    dist     = dist_indexing(dist_i)

    intlabel = create_label_html(obsid_i, dist_i)
    [fig_hi, fit1_hi, fit2_hi] = manage_plots(xmin, xmax, yMinSets, yMaxSets, \
                xSets, ySets, dist,  xname, yname,  entLabels, intlabel=intlabel, inst = 'hrc_i')
#
#--- create an interactive plot for HRC S
#
    xSets    = [time_s, time_s, time_s]
    ySets    = [med_s, center_s, width_s]
    dist     = dist_indexing(dist_s)

    intlabel = create_label_html(obsid_s,  dist_s)
    [fig_hs, fit1_hs, fit2_hs] =  manage_plots(xmin, xmax, yMinSets, yMaxSets, \
                xSets, ySets, dist,  xname, yname,  entLabels, intlabel=intlabel, inst='hrc_s')
#
#---- update the main html page
#
    update_html_page(fig_hi, fig_hs, fit1_hi, fit2_hi, fit1_hs, fit2_hs, time_i)
#
#---- create radio html pages
#
    create_radial_html_page(time_i, obsid_i, dist_i, med_i, center_i, width_i, end_year, inst='i')
    create_radial_html_page(time_s, obsid_s, dist_s, med_s, center_s, width_s, end_year, inst='s')
#
#--- create energy gain fitting result table page
#
    gain_fitting_table()

#---------------------------------------------------------------------------------------------------
#-- dist_indexing: convert radial distance interval into index                                   ---
#---------------------------------------------------------------------------------------------------

def dist_indexing(dist):
    """
    convert radial distance interval into index
            dist < 5         ---> 0
            5< dist < 10    ----> 1
            10 < dist < 15  ----> 2
            15 < dist       ----> 3
    """
    dist_index = []
    for ent in dist:
        if ent < 5.0:
            dist_index.append(0)
        elif ent >=5 and ent < 10:
            dist_index.append(1)
        elif ent >= 10 and ent < 15:
            dist_index.append(2)
        else:
            dist_index.append(3)

    return dist_index

#---------------------------------------------------------------------------------------------------
#-- manage_plots: create three panel interactive plots                                           ---
#---------------------------------------------------------------------------------------------------

def manage_plots(xmin, xmax, yMinSets, yMaxSets, xSets, ySets, dist,  xname,\
                 yname, entLabels, intlabel='', extraNote='', inst='hrc_i'):
    """
    create three panel interactive plots
    Input:  xmin        --- xmin
            xmax        --- xmax
            yMinSets    --- a list of ymin 
            yMaxSets    --- a list of ymax
            xSets       --- a list of lists containing x-axis data
            ySets       --- a list of lists containing y-axis data
            dist:       --- category of the data (see: dist_indexing)
            xname       --- x axis name
            yname       --- y axis name
            entLabels   --- a list of the names of each data
            intlabel    --- the list of the informaiton displayed on the interactive plot
            extraNote   --- an extra Note you want to add

    Output: fout        --- interactive plot
            flist1      --- the line fitted result before 2012
            flist2      --- the line fitted results after 2012
    """
    fout = ''
    flist1 = []
    flist2 = []
    for k in range(0, len(xSets)):
        ymin = yMinSets[k]
        ymax = yMaxSets[k]
        xdata = xSets[k]
        ydata = ySets[k]
        elabel = entLabels[k]

        [fig, fit1, fit2] = plotPanel(xmin, xmax, ymin, ymax, xdata, ydata,\
                                dist,  xname, yname, elabel, intlabel, extraNote, inst)

#
#--- converting mp plot data into html format
#
        fout   = fout + mpld3.fig_to_html(fig)

        flist1 = flist1 + fit1
        flist2 = flist2 + fit2

    return [fout, flist1, flist2]

#---------------------------------------------------------------------------------------------------
#--- plotPanel: create a single pnael plot                                                       ---
#---------------------------------------------------------------------------------------------------

def plotPanel(xmin, xmax, ymin, ymax, xdata, ydata,  dist,  xname, yname,\
              elabel, intlabel='', extraNote='', inst='hrc_i'):
    """
    create a single pnael plot
    Input:  xmin        --- xmin
            xmax        --- xmax
            ymin        --- ymin 
            ymax        --- ymax
            xdata       --- a list of x-axis data
            ydata       --- a list of y-axis data
            dist        --- category of the data (see: dist_indexing)
            xname       --- x axis name
            yname       --- y axis name
            elabel      --- a list of the names of each data
            intlabel    --- the list of the informaiton displayed on the interactive plot
            extraNote   --- an extra Note you want to add

    Output: fig         --- interactive plot
            fit1        --- the line fitted result before 2012
            fit2        --- the line fitted results after 2012
    """
#
#--- set line propetry etc 
#
    marker_list = ['s', 'o', 'D', 'v']
    marker_size = [60, 70, 60, 70]
#
#--- this css is used for the pop up page
#
    css = """
        body{
            width:600px;
            height:300px;
        }
        p{
            text-align:center;
        }
        """
#
#--- clean up the plotting device
#
    plt.close('all')
#
#--- set plotting frames

    fig, ax = plt.subplots()
#
#---- set a few parameters
#
    mpl.rcParams['font.size'] =12 
    props = font_manager.FontProperties(size=12)
    plt.subplots_adjust(hspace=0.08)
#
#--- separate them according to the distances
#
    xv   = [[],[],[],[]]
    yv   = [[],[],[],[]]
    labv = [[],[],[],[]]
#
#--- shift x position slightly so that the data points plotted won't overlap too closely
#
    for k in range(0, len(dist)):
        dm = dist[k]
        if dm == 0:
            xval = xdata[k] - 0.10
        elif dm == 1:
            xval = xdata[k]
        elif dm == 2:
            xval = xdata[k] + 0.10
        else:
            xval = xdata[k] + 0.20

        xv[dm].append(xval)
        yv[dm].append(ydata[k])
        labv[dm].append(str(intlabel[k]))
#
#---- data plotting 
#
    fit1 = []
    fit2 = []
    for m in range(0, 4):

        pv = ax.scatter(xv[m], yv[m], color=colorList[m], marker=marker_list[m],\
                           s=marker_size[m], lw =0)
#
#--- connect the data points to the descriptions
#
        plugins.connect(fig, mpld3.plugins.PointHTMLTooltip(pv, labv[m],  css=css, hoffset=-140))

        shifted1 = []
        yv1      = []
        shifted2 = []
        yv2      = []
        b_point  = 2012                     #--- separate the data before and after 2012.1.1
        for n in range(0, len(xv[m])):
            if xv[m][n] < b_point and xv[m][n] >= 2000:
#
#--- to get the better fit, shift the data to start from 1999 rather than year 0
#
                shifted1.append(xv[m][n] - 1999)
                yv1.append(yv[m][n])

            elif xv[m][n] >= b_point:
                shifted2.append(xv[m][n] - 1999)
                yv2.append(yv[m][n])
#
#--- fitting the first half
#
        paramsinitial = [80, 0.3]
        ferr   = [0] * len(shifted1)
        fitobj = fit_line(paramsinitial, shifted1, yv1, ferr, 'linear')
        (a, b) = fitobj.params
        (ae,be)= fitobj.xerror
        fit1.append([a, b, be])
        start  = a + b * ((xmin-1999) - 0.5)
        stop   = a + b * 14
        ax.plot([xmin, b_point], [start, stop],  color=colorList[m-1],\
                    marker=marker_list[m-1], lw=2, label=label_list[m-1])
#
#--- fitting the latter half
#
        ferr   = [0] * len(shifted2)
        fitobj = fit_line(paramsinitial, shifted2, yv2, ferr, 'linear')
        (a, b) = fitobj.params
        (ae,be)= fitobj.xerror
        fit2.append([a, b, be])
        start  = a + b * 14
        stop   = a + b * ((xmax-1999) - 0.5)
        ax.plot([b_point, xmax], [start, stop],  color=colorList[m-1],\
                    marker=marker_list[m-1], lw=2)
#
#---- set plotting frames
#
    ax.legend(loc='upper left',  ncol=2, bbox_to_anchor=(0.4, 0.95), fancybox=True, shadow=True)
    xdiff = xmax - xmin
    xtext = xmin + 0.05 * xdiff

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.text(xtext, y_text(ymin, ymax), elabel, fontsize=12)
#
#--- add x label on the bottom panel
#
    ax.set_ylabel(yname)
    ax.yaxis.set_label_coords(-0.10, 0.5)
    ax.set_xlabel('Time (year)')
#
#--- set the size of plot
#
    fig.set_size_inches(10.0, 5.0)
    fig.tight_layout()

    plt.close('all')

    return [fig, fit1, fit2]

#---------------------------------------------------------------------------------------------------
#-- manage_plots2: plots multiple data in separate panels: this is for the radial distribution    --
#---------------------------------------------------------------------------------------------------

def manage_plots2(xmin, xmax, yMinSets, yMaxSets, xSets, ySets, xname, yname, entLabels, plabels, year):
    """
    This function plots multiple data in separate panels: this is for the radial distribution
    Input:  xmin        --- xmin
            xmax        --- xmax
            yMinSets    --- a list of ymin 
            yMaxSets    --- a list of ymax
            xSets       --- a list of lists containing x-axis data
            ySets       --- a list of lists containing y-axis data
            xname       --- x axis name
            yname       --- y axis name
            entLabels   --- a list of the names of each data
            plabel      --- the list of the informaiton displayed on the interactive plot
            year        --- the year of the dat

    Output: fout        --- interactive plot
            flist       --- the line fitted result 
    """
    fout = ''
    flist = []
    for k in range(0, len(xSets)):
        ymin = yMinSets[k]
        ymax = yMaxSets[k]
        xdata = xSets[k]
        ydata = ySets[k]
        elabel = entLabels[k]

        [fig, fit] = plotPanel2(xmin, xmax, ymin, ymax, xdata, ydata, xname,\
                                yname, elabel,plabels, year, k)
        fout = fout + mpld3.fig_to_html(fig)
        flist = flist + fit

    return [fout, flist]

#---------------------------------------------------------------------------------------------------
#--- plotPanel2: crate a single panel plot: radial distribution version                          ---
#---------------------------------------------------------------------------------------------------

def plotPanel2(xmin, xmax, ymin, ymax, xdata, ydata, xname, yname, elabel, plabels, year, ik):

    """
    create a single panel plot: this is for the radial distribution
    Input:  xmin        --- xmin
            xmax        --- xmax
            ymin        --- ymin 
            ymax        --- ymax
            xdata       --- a list of x-axis data
            ydata       --- a list of y-axis data
            xname       --- x axis name
            yname       --- y axis name
            elabel      --- a list of the names of each data
            plabel      --- the list of the informaiton displayed on the interactive plot
            year        --- the year of the dat
    Output: fig         --- interactive plot
            fit:        --- the line fitted result 
    """
#
#--- set line  properties
#
    marker_list = ['s', 'o', 'D', 'v']
    marker_size = [60, 70, 60, 70]
#
#--- this css is used for the pop up page
#
    css = """
        body{
            width:600px;
            height:300px;
     x   }
        p{
            text-align:center;
        }
        """
    xdiff  = xmax - xmin
    xtext  = xmin + 0.05 * xdiff
    xtext2 = xmin + 0.50 * xdiff
#
#--- clean up the plotting device
#
    plt.close('all')
#
#--- set plotting frames

    fig, ax = plt.subplots()
#
#---- set a few parameters
#
    mpl.rcParams['font.size'] =12 
    props = font_manager.FontProperties(size=12)
    plt.subplots_adjust(hspace=0.08)

    ytext  = ymax -0.1 * (ymax -ymin)
#
#---- data plotting 
#
    fit1 = []
    pv = ax.scatter(xdata, ydata, color=colorList[ik], marker=marker_list[ik], 
                       s=marker_size[ik], lw =0)
#
#---- connect the data points to the descriptions
#
    plugins.connect(fig, mpld3.plugins.PointHTMLTooltip(pv, plabels,  css=css, hoffset=-140))

    paramsinitial = [100, 1.0]
    if len(xdata) > 2:
        ferr   = [0] * len(xdata)
        try:
            fitobj = fit_line(paramsinitial, xdata, ydata, ferr, 'linear')
            (a, b) = fitobj.params
            (ae,be)= fitobj.xerror
            fit1.append([a, b, be])
            start  = a + b * xmin
            stop   = a + b * xmax
            ax.plot([xmin, xmax], [start, stop],  color=colorList[ik], lw=2)
#
#--- adding the linear fit result on the plotting 
#
            ca  = "%2.2f" % (round(a,  2))
            cb  = "%2.2f" % (round(abs(b),  2))
            cbe = "%2.2f" % (round(be, 2))
            if b < 0:
                line = 'Y = ' + ca + ' - (' + cb + '+/-' +  cbe + ') * X'
            else:
                line = 'Y = ' + ca + ' + (' + cb + '+/-' +  cbe + ') * X'
            ax.text(xtext2, ytext, line, fontsize=12)
        except:
            pass
#
#---- set plotting frames
#
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.text(xtext, y_text(ymin, ymax), elabel, fontsize=12)
#
#--- add x label on the bottom panel
#
    ax.set_ylabel(yname)
    ax.yaxis.set_label_coords(-0.10, 0.5)
    ax.set_xlabel('Radial Distance (Arcmins)')
#
#--- set the size of plot
#
    fig.set_size_inches(10.0, 5.0)
    fig.tight_layout()

    plt.close('all')

    return [fig, fit1]

#---------------------------------------------------------------------
#-- y_text: set y position of text                                 ---
#---------------------------------------------------------------------

def y_text(ymin, ymax):
    """
    set y position of text
    input:  ymin    --- min of y
            ymax    --- max of y
    output  ytext   --- y position of the text
    """
    ydiff = ymax - ymin
    ytext = ymin + 0.10 * ydiff

    return ytext

#---------------------------------------------------------------------------------------------------
#-- create_label_html: create a html page to be displayed for the data point                       --
#---------------------------------------------------------------------------------------------------

def create_label_html(obsid, dist, pos=1):
    """
    create a html page to be displayed for the data point
    input:  obsid   --- a list of obsids
            dist    --- a list of the distances from the center
            pos     --- two possible realtive paths to the indivisual plots
    output: hlist   --- a list of the html page content
    """
#
#--- set link information
#
    fdir      = plot_dir + 'Indivisual_Plots/'
    if pos == 1:
        html_plot = './' + html_trend + 'Indivisual_Plots/'
    else:
        html_plot = '../Indivisual_Plots/'

    hlink     = '<p> <img src="' + html_plot
#
#--- create a page for each obsid
#
    hlist = []
    for k in range(0, len(obsid)):
        title = '<h3 style="background-color:yellow">ObsId: ' + str(obsid[k])
        title = title + ' Distance: ' + str(dist[k]) + '</h3>'

        cobsid = adjust_obsid_digit(obsid[k])
        hfile  = 'hrcf' + cobsid + '_gfit.png'
#
#--- check whether the plot exists
#
        chk = fdir +  hfile
        if os.path.isfile(chk):
            alink = title + hlink  +  hfile + '" width=600px> </p>'
        else:
            alink = title + '<p>No Plot </p>'

        hlist.append(alink)

    return hlist
            
#---------------------------------------------------------------------------------------------------
#-- adjust_obsid_digit: adjust obsid digit for display                                            --
#---------------------------------------------------------------------------------------------------

def adjust_obsid_digit(obsid):
    """
    adjust obsid digit for display
    obsid   --- obsid
    cobsid  --- digit adjusted obsid
    """
    cval  = float(obsid)

    if cval < 10:
        cobsid = '0000' + str(obsid)
    elif cval < 100:
        cobsid = '000' + str(obsid)
    elif cval < 1000:
        cobsid = '00' + str(obsid)
    elif cval < 10000:
        cobsid = '0' + str(obsid)
    else:
        cobsid = str(obsid)

    return cobsid

#---------------------------------------------------------------------------------------------------
#-- update_html_page: update the main html page                                                    -
#---------------------------------------------------------------------------------------------------

def update_html_page(fig_hi, fig_hs,fit1_hi, fit2_hi, fit1_hs, fit2_hs, stime):
    """
    update the main html page
    input;  fig_hi      --- interactive plot for hrc i 
            fig_hs      --- interactive plot for hrc s
            fits1_hi    ---- line fitting result for hrc i before 2012
            fits2_hi    ---- line fitting result for hrc i after 2012
            fits1_hs    ---- line fitting result for hrc s before 2012
            fits2_hs    ---- line fitting result for hrc s after 2012
            stime       --- a list of data time in year

    output: <web_dir>/arlac_energy_trend.html
    """
#
#--- read the template for the top part of the page
#
    syear = sorted(stime)
    lyear = int(float(syear[-1])) + 1
    infile = house_keeping + 'hrc_trend_top_part'

    with open(infile, 'r') as f:
        out    = f.read()
#
#--- hrc i
#
    out    = out + '<h3 style="padding-top:10px;">HRC I</h3>\n'
    out    = out + '<p style="padding-bottom:40px">\n'
#
#--- interactive plot of hrc i
#
    out    = out + '<div style="padding-left:40px;">\n'

    #out    = out + mpld3.fig_to_html(fig_hi)
    out    = out + fig_hi

    out    = out + '</div>\n'

    out    = out + '</p>\n'

    out    = out + '<h3>Fitted Lines (Intercept is at Year 1999)</h3>'
#
#--- add line fitting results
#
    out    = out + '<table border=1 cellpadding=4 cellspacing=4 '
    out    = out + 'style="padding-top:20px;padding-bottom:20px;">\n'
    out    = out + '<tr>\n'
    out    = out + '<th>&#160;</th><th colspan=3>2000-2012</th><th colspan=3>2012-Current</th>\n'
    out    = out + '</tr>\n'
    out    = out + '<th>&#160;</th><th>Intercept</th><th>Slope</th><th>Slope Error</th>'
    out    = out + '<th>Intercept</th><th>Slope</th><th>Slope Error</th>\n'
    out    = out + '</tr>\n'


    for m in range(0, 4):
        out    = out + '<tr>\n'
        out    = out + '<th>' + label_list[m] + '</th>'
        for l in range(0, 3):
            val = '%2.2f' % round(fit1_hi[m][l], 2)
            out    = out +  '<td style="text-align:center">' + val + '</td>'
        for l in range(0, 3):
            val = '%2.2f' % round(fit2_hi[m][l], 2)
            out    = out +  '<td style="text-align:center">' + val + '</td>'
        out    = out + '</tr>\n'
    out    = out + '</table>\n'
#
#--- add link to the radial distribution plots
#
    out    = out + '<h3>HRC I: Radial Distribution of Each Year (Select year to open the plot)</h3>\n'
    out    = out + '<table border=1 cellpadding=4 cellspace=2>\n'
    out    = out + '<tr>\n'
    diff   = lyear - 1999
    for  k in range(0, diff):
        year = 1999 + k
        #out = out + '<td><a href="https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/'
        out = out + './'
        out = out + html_trend + 'Dist_Html/hrc_i_radial_dist_year'
        out = out + str(year) + '.html">' + str(year) + '</a></td>\n'
        if (k+1) % 14 == 0:
            out = out + '</tr><tr>\n'
    k = diff-1
    while((k+1) % 14 != 0):
        out = out + "<td>&#160;</td>\n"
        k += 1
    out    = out + '</tr>\n</table>\n'
#
#--- hrc s
#
    out    = out + '<h3 style="padding-top:10px;">HRC S</h3>\n'
    out    = out + '<p style="padding-bottom:40px">\n'

    out    = out + '<div style="padding-left:40px;">\n'

    #out    = out + mpld3.fig_to_html(fig_hs)
    out    = out +fig_hs

    out    = out + '</div>\n'
    out    = out + '</p>\n'

    out    = out + '<h3>Fitted Lines (Intercept is at Year 1999)</h3>'
#
#--- line fitted results table
#
    out    = out + '<table border=1 cellpadding=4 cellspacing=2 '
    out    = out + 'style="padding-top:20px;padding-bottom:20px;">\n'
    out    = out + '<tr>\n'
    out    = out + '<th>&#160;</th><th colspan=3>2000-2012</th><th colspan=3>2012-Current</th>\n'
    out    = out + '</tr>\n'
    out    = out + '<th>&#160;</th><th>Intercept</th><th>Slope</th><th>Slope Error</th>'
    out    = out + '<th>Intercept</th><th>Slope</th><th>Slope Error</th>\n'
    out    = out + '</tr>\n'

    for m in range(0, 4):
        out    = out + '<tr>\n'
        out    = out + '<th>' + label_list[m] + '</th>'
        for l in range(0, 3):
            val = '%2.2f' % round(fit1_hs[m][l], 2)
            out    = out +  '<td style="text-align:center">' + val + '</td>'
        for l in range(0, 3):
            val = '%2.2f' % round(fit2_hs[m][l], 2)
            out    = out +  '<td style="text-align:center">' + val + '</td>'
        out    = out + '</tr>\n'
    out    = out + '</table>\n'
#
#--- link to the radial distribution plots
#
    out    = out + '<h3>HRC S: Radial Distribution of Each Year (Select year to open the plot)</h3>\n'
    out    = out + '<table border=1 cellpadding=4 cellspace=2>\n'
    out    = out + '<tr>\n'
    diff   = lyear - 1999
    for  k in range(0, diff):
        year = 1999 + k
        #out = out + '<td><a href="https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/'
        out = out + './'
        out = out + html_trend + 'Dist_Html/hrc_i_radial_dist_year'
        out = out + str(year) + '.html">' + str(year) + '</a></td>\n'
        if (k+1) % 14 == 0:
            out = out + '</tr><tr>\n'
    k = diff-1
    while((k+1) % 14 != 0):
        out = out + "<td>&#160;</td>\n"
        k += 1
    out    = out + '</tr>\n</table>\n'
#
#--- rest...
#
    rfile  = house_keeping + 'hrc_trend_bottom'
    with  open(rfile, 'r') as fi:
        rest   = fi.read()

    out    = out + rest

    out    = out.replace('None', '')        #--- fixing a bug to remove un-wanted "None" appears on the web page
    ctime  = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    out    = out.replace('#CTIME#', ctime)  #--- update last modified date
#
#--- print out the html page
#
    outdir = web_dir + 'arlac_energy_trend.html'

    fo     = open(outdir, 'w')
    fo.write(out)
    fo.close()

#---------------------------------------------------------------------------------------------------
#-- create_radial_html_page: create html pages for radial distribution plots                     ---
#---------------------------------------------------------------------------------------------------

def create_radial_html_page(time_e, obsid_e, dist_e, med_e, center_e, width_e, end_year, inst='i'):
    """
    create html pages for radial distribution plots
    input:  time_e  --- a list of time (in frational year)
            obsid_e --- a list of obsids
            dist_e  --- a list of distances
            med_e   --- a list of median pha
            center_e    --- a list of voigt estimated center
            width_e     --- a list of fwhm
            end_year    --- the last year in the data
            inst        --- instrument either "i" or "s"
    output: <web_dir>/Trend_plots/Dist_Html/hrc_<inst>_radial_dist_year<year>.html'
    """
#
#--- trend along radial distance
#
    xmin  = 0
    xmax  = 25
    xname = 'Radial Distance(Arcmins)'
    stop  = end_year -1
#
#--- separate data by year
#
    for pyear in range(1999, end_year):
        dist   = []
        med    = []
        center = []
        width  = []
        obsid  = [] 
        pend = pyear + 1
        for j in range(0, len(time_e)):
            if  (time_e[j] >= pyear) and (time_e[j] < pend):
                obsid.append(obsid_e[j])
                dist.append(dist_e[j])
                med.append(med_e[j])
                center.append(center_e[j])
                width.append(width_e[j])

        tlen   = len(obsid)

        plabels  = create_label_html(obsid, dist, pos=2)
        xSets    = [dist, dist, dist]
        ySets    = [med, center, width]
#
#--- create the interactive plot
#
        if tlen > 0:
            [fig, fit] = manage_plots2(xmin, xmax, yMinSets, yMaxSets, xSets, ySets,\
                                    xname, yname, entLabels, plabels, pyear)
#
#--- create the html page
#
        out = '<!DOCTYPE html>\n'
        out = out + '<html>\n'
        out = out + '<head>\n'
        out = out + '<title> HRC Energy Trending: Radial Distribution</title>\n'
        out = out + '</head>\n'
        out = out + '<body style="background-color:#F5F5DC;width:95%;margin-left:10px; '
        out = out + 'margin-right;10px">\n'
        out = out + '<h3> Radial Distribution. Year: ' + str(pyear) + '</h3>\n'
        out = out + '<div style="text-align: right;">\n'
#
#--- add links to move forward, backward, and back to the top page
#
        if pyear > 1999 and pyear < stop:
            out = out + '<a href="' + html_trend + 'Dist_Html/hrc_' + inst 
            out = out + '_radial_dist_year' + str(pyear-1) + '.html">'
            out = out + '<b>&lt;&lt;Previous Year</b></a>'
            #out = out + '<span style="color:#F5F5DC">span</span>'
            out = out + ' &diams; '
            out = out + '<a href="' + html_trend + 'Dist_Html/hrc_' + inst 
            out = out + '_radial_dist_year' + str(pyear+1) + '.html">'
            out = out + '<b>Next Year&gt;&gt;</b></a>'
        elif pyear == 1999:
            out = out + '<a href="' + html_trend + 'Dist_Html/hrc_' + inst 
            out = out + '_radial_dist_year' + str(pyear+1) + '.html">'
            out = out + '<b>Next Year&gt;&gt;</b></a>'
        elif pyear == stop:
            out = out + '<a href="' + html_trend + 'Dist_Html/hrc_' + inst 
            out = out + '_radial_dist_year' + str(pyear-1) + '.html">'
            out = out + '<b>&lt;&lt;Previous Year</b></a>'

        out = out + '<span style="color:#F5F5DC">span</span>'
        out = out + '<a href="https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/'
        out = out + 'ArLac/arlac_energy_trend.html">'
        out = out + '<b>Back to Main Page</b></a>\n</div>\n'
#
#--- if there are data, add the interactive plot; else say "No Data"
#
        if tlen > 0:
            #out = out + mpld3.fig_to_html(fig)
            out = out + fig

            out = out + '<div style="text-align: right;padding-bottom:40px">\n'
            out = out + '<a href="https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/'
            out = out + 'ArLac/arlac_energy_trend.html">'
            out = out + '<b>Back to Main Page</b></a>\n</div>\n'
        else:
            out = out + '<h3 style="padding-top:50px; padding_bottom:50px;">No Data</h3>\n'

        out = out + '</body>\n'
        out = out + '</html>\n'

        ofile = plot_dir + '/Dist_Html/hrc_' + inst + '_radial_dist_year' 
        ofile = ofile   + str(pyear) + '.html'

        with  open(ofile, 'w') as fo:
            fo.write(out)

#-----------------------------------------------------------------------------------------
#-- gain_fitting_table: update voigt proflie fitting result page                       ---
#-----------------------------------------------------------------------------------------

def gain_fitting_table():
    """
    update voigt proflie fitting result page
    input: none, but read from <data_dir>/fitting_results
    output: <web_dir>/enegy_fitting_result.html
    """
    gfile = data_dir + 'fitting_results'
    data  = mcf.read_data_file(gfile)
    data  = sort_by_col(data, 2)
#
#--- create the html page
#
    efile = house_keeping + 'enegy_fitting_result_top'
    with open(efile, 'r') as f:
        out   = f.read()

    for ent in data:
        atemp = re.split('\s+', ent)
        out = out + '<tr>\n'
        for k in range(0, len(atemp)):
            if k == 0:
                out = out + '<td><a href="javascript:WindowOpener(\''  
                out = out + 'hrcf' + adjust_obsid_digit(atemp[0]) + '_gfit.png'
                out = out + '\')" style="text-align:right">' + atemp[0] + '</a></td>\n'
            elif k == 2:
                continue
            elif k >= 13:
                break
            elif k == 8:
                val = float(atemp[8])
                if val < 5:
                    color = 'aqua'
                elif  val >= 5 and val< 10:
                    color = 'lime'
                elif  val >= 10 and val < 15:
                    color = 'fuchsia'
                elif val > 15:
                    color = 'red'
                out = out + '<td style="background-color:' + color + ';">' + atemp[8] + '</td>\n'
            else:
                out = out + '<td>' + atemp[k] + '</td>\n'

        out = out + '</tr>\n'

    out = out + '</table>\n'

    outfile = plot_dir + '/enegy_fitting_result.html'
    with open(outfile, 'w') as fo:
        fo.write(out)

#-----------------------------------------------------------------------------------------
#-- sort_by_col: sort the list of list by a column                                      --
#-----------------------------------------------------------------------------------------

def sort_by_col(data,colnum):
    """
    sort the list of list by a column
    input:  data    --- a list of lists
            colnum  --- the position of the entry which you want to use for sorting
    output: out     --- the sorte list of list
    """

    sdict = {}
    clist = []
    for ent in data:
        atemp = re.split('\s+', ent)
        val   = float(atemp[colnum])
        clist.append(val)
        sdict[val] = ent

    clist.sort()

    out   = []
    for ent in clist:
        out.append(sdict[ent])

    return out

#-----------------------------------------------------------------------------------------
#-- fit_line: kmpfit calling function to fit the lines on data                          --
#-----------------------------------------------------------------------------------------

def fit_line(paramsinit, x, y, err, ptype):
    """
    kmpfit calling function to fit the lines on data
    input: paramsinit: initial guess for the parameters
    x, y: data
    ptype: linear or exp:
    output: two parameters (a, b) are returned
    """
    sx = []
    sy = []
    se = []
    avg  = mean(y)
    stdp = std(y)
    bot  = avg - 3.0 * stdp
    top  = avg + 3.0 * stdp
    i = 0
    for val in y:
        if (val >= bot) and (val <= top):
            sx.append(x[i])
            sy.append(y[i])
            se.append(err[i])
        i += 1
#
#--- make sure that the arrays are numpyed
#
    d = numpy.array(sx)
    v = numpy.array(sy)
    e = numpy.array(se)

    if ptype == 'linear':
#
#--- linear fit
#
        fitobj = kmpfit.Fitter(residuals=linear_res, data=(d, v, e))
        fitobj.fit(params0 = paramsinit)
    else:
#
#--- exp fit
#
        fitobj = kmpfit.Fitter(residuals=exp_res, data=(d, v, e))
        fitobj.fit(params0 = paramsinit)
    
    return fitobj

#-----------------------------------------------------------------------------------------
#-- linear_fit: linear model fit                                                       ---
#-----------------------------------------------------------------------------------------

def linear_fit(param, x):
    """
    linear model fit
    input: param: (a,b)
    x independent val
    ouptput: estimate
    """
    a, b = param

    return (a + b * x)

#-----------------------------------------------------------------------------------------
#-- linear_res: linear model resitual                                                   --
#-----------------------------------------------------------------------------------------

def linear_res(param, data):
    """
    linear model resitual
    input: param (a, b)
    data  (x, y)
    output: residual
    """
    a, b= param
    x, y, e = data
    
    res =  y - (a + b * x)
    return res

#-----------------------------------------------------------------------------------------
#-- exp_fit: exponential model                                                          --
#-----------------------------------------------------------------------------------------

def exp_fit(param, x):
    """
    exponential model
    input: param (a, b)
    x  independent variable
    output: estimate
    """
    a, b = param

    #return  (a * exp(-1.0 * b * x))
    return  (a * (1.0 - exp(-1.0 * b * x)))

#-----------------------------------------------------------------------------------------
#-- exp_res: exponential model residual                                                ---
#-----------------------------------------------------------------------------------------

def exp_res(param, data):
    """
    exponential model residual
    input param(a, b)
    data (x, y)
    output: residual
    """
    a, b= param
    x, y, e = data
    
    res = y - (a * exp(-1.0 * b * x))
    res = y - exp_fit(param, x)
    return res

#--------------------------------------------------------------------

if __name__ == '__main__':
    hrc_gain_trend_plot()


