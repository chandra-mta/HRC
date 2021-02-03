#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#################################################################################
#                                                                               #
#       create_html_page.py: update arlac monitoring html page                  #
#                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                           #
#                                                                               #
#           Last Update: Feb 03, 2021                                           #
#                                                                               #
#################################################################################

import sys
import os
import string
import re
import math
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
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts/house_keeping/dir_list'

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
sys.path.append(hrc_common)

import mta_common_functions     as mcf
import hrc_common_functions     as hcf
#
#--- temp writing file name
#
import random
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)
t_list = ['On Axis', 'Y-offset = +5 <-> 15 arcmin', 'Y-offset = -5 <-> -15  arcmin',\
          'Y-offset &ge; +15 arcmin or farther', 'Y-offset &le; -15 arcmin or farther']
p_list = ['i_0', 'i_10', 'i_m10', 'i_25', 'i_m25', 's_0', 's_10', 's_m10', 's_25', 's_m25']

admin  = 'tisobe@cfa.harvsard.edu'

#-----------------------------------------------------------------------------------------
#-- create_html_and_plot: update arlac monitoring html page                              --
#-----------------------------------------------------------------------------------------

def create_html_and_plot():
    """
    update arlac monitoring html page
    input:  none:
    output: <html_page>/arlac_vis_montior.html
            <html_dir>/Plots/*.png
    """
#
#--- read html page template
#
    template = house_keeping + '/arlac_template'
    with open(template, 'r') as f:
        text     = f.read()
#
#--- go though different settings; set input data, plot file name
#
    for k in range(0, len(p_list)):
        pos = p_list[k]
        dfile = data_dir + 'hrc_' + pos + '_results'
        out   = html_dir + 'Plots/hrc_' + pos + '_rate.png'
        try:
            data  = mcf.read_data_file(dfile)
        except:
            continue
#
#--- plot data
#
        if len(data) >= 3:
            plot_data(data, out, pos)
#
#--- create a table data then create a html page for the table
#
            if k < 5:
                inst = 'HRC-I: '
                m = 0
            else:
                inst = 'HRC-S: '
                m = 5
            title  = inst + t_list[k - m]
            line   = create_html_table(data, title)
            ofile  = html_dir + 'Sub_htmls/hrc_' + p_list[k] + '.html'
            with open(ofile, 'w') as fo:
                fo.write(line) 
#
#--- change the updated date
#
    ctime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    text = text.replace('#UPDATE#', ctime)
#
#--- print out the html page
#
    hfile = html_dir + 'arlac_vis_montior.html'
    with open(hfile, 'w') as fo:
        fo.write(text)
#
#--- send notification
#
    send_notification()

#-----------------------------------------------------------------------------------------
#-- plot_data: plotting data                                                            --
#-----------------------------------------------------------------------------------------

def plot_data(data, outname, pos):
    """
    plotting data
    input:  data    --- data
            outname --- output file name
            pos     --- indicator of position of aiming
    output  outname --- png plot of the data
    """
    date = []
    cnt  = []
    err  = []
    for ent in data:
        atemp = re.split('\s+', ent)
        if not mcf.is_neumeric(atemp[-1]):
            continue
        try:
            ytime = hcf.convert_time_to_fyear(atemp[2], tformat='%Y-%m-%dT%H:%M:%S')
    
            exp   = float(atemp[3])
            val   = float(atemp[4])
            if val < 10:
                continue
            val  /= exp
            sig   = float(atemp[5])/exp
    
            date.append(ytime)
            cnt.append(val)
            err.append(sig)
        except:
            continue
#
#--- fit a linear line
#
    [a, b, sa, sb] =  line_fit(date, cnt, err)
#
#--- set min max
#
    xmin = 1999.0
    ta   = time.localtime()
    xmax = ta.tm_year + 1
    ymin = 0.0
    ymax = 6
#
#--- start plotting
#
    plt.close('all')
    mpl.rcParams['font.size'] = 9
    props = font_manager.FontProperties(size=9)

    ax  = plt.subplot(111)
    ax.set_autoscale_on(False)
    ax.set_xbound(xmin,xmax)
    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)

    plt.errorbar(date, cnt, yerr=err, fmt='o', lw=1)
    plt.xlabel('Time (year)')
    plt.ylabel('Source Rate (cnt/s)')
#
#--- plot fitting line
#
    start = a + b * xmin
    stop  = a + b * xmax
    plt.plot([xmin, xmax], [start, stop], lw =1, color='blue')
#
#--- write the fitting equation
#
    xdiff = xmax - xmin
    ydiff = ymax - ymin
    xpos  = xmin + 0.1  * xdiff
    ypos  = ymax - 0.1  * ydiff
    ac    = '%1.2f' % (round(a,  2))
    if abs(b) < 0.01:
        bc    = '%1.4f' % (round(abs(b),  4))
        ec    = '%1.4f' % (round(sb, 4))
    else:
        bc    = '%1.2f' % (round(abs(b),  2))
        ec    = '%1.2f' % (round(sb, 2))
    if b > 0 :
        text  = '(Source Rate) = ' + ac  + ' + (' + bc + '+/-' + ec + ') * Time'
    else:
        text  = '(Source Rate) = ' + ac  + ' - (' + bc + '+/-' + ec + ') * Time'
    plt.text(xpos, ypos, text)
#
#--- set the size of the plotting area in inch (width: 10.0in, height 5 in)
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 5.0)

#
#--- save the plot in png format
#
    plt.savefig(outname, format='png', dpi=100)

#-----------------------------------------------------------------------------------------
#-- line_fit: fit a weighted linear line fit                                            --
#-----------------------------------------------------------------------------------------

def line_fit(x, y, e):
    """
    fit a weighted linear line fit
    input:  x       --- independent data
            y       --- dependent data
            e       --- y error
    output: a       --- intercept
            b       --- slope
            siga    --- error on the intercept
            sigb    --- error on the slope
    """
    suma  = 0
    sumx  = 0
    sumy  = 0
    sumx2 = 0
    sumy2 = 0
    sumxy = 0

    dlen = len(x)
    if dlen < 3:
        return [0, 0, 0, 0]

    for k in range(0, dlen):
        try:
            weight = 1.0 / e[k]**2
        except:
            weight = 1.0
        suma  += weight
        sumx  += weight * x[k]
        sumy  += weight * y[k]
        sumx2 += weight * x[k] * x[k]
        sumy2 += weight * y[k] * y[k]
        sumxy += weight * x[k] * y[k]

    delta = suma * sumx2 - sumx* sumx
    a     = (sumx2 * sumy - sumx * sumxy) / delta
    b     = (sumxy * suma - sumx * sumy ) / delta
    if dlen <= 2:
        siga = 0
        sigb = 0
    else:    
        var   = (sumy2 + a * a * suma + b * b * sumx2 \
                 - 2.0 *(a * sumy + b * sumxy - a * b * sumx)) / (len(x) -2)

        siga  = math.sqrt(var * sumx2 / delta)
        sigb  = math.sqrt(var * suma  / delta)

    return [a, b, siga, sigb]

#-----------------------------------------------------------------------------------------
#-- create_html_table: create a html table based on the given data                      --
#-----------------------------------------------------------------------------------------

def create_html_table(data, title):
    """
    create a html table page based on the given data
    input:  data    --- a table data
            title   --- a title of the table
    output: line    --- a html table  page
    """
    line = '<!DOCTYPE html>\n'
    line = line + '<html>\n<head>\n'
    line = line + '<title>Monitoring the UV/Ion Shield Health: Ar Lac Data Table</title>\n'
    line = line + '</head>\n'
    line = line + "<body style='background-color:#F5F5DC;width:95%;margin-left:10px;"
    line = line + "margin-right;10px'>"
    line = line + '<h2>' + title + '</h2>\n'
    line = line + '<div style="padding-bottom 10px;"></div>'

    line = line + '<table border=1 cellpadding=2 cellspacing=2>\n'
    line = line + '<tr>\n'
    line = line + '<th>ObsID</th>\n'
    line = line + '<th>Filename</th>\n'
    line = line + '<th>Date</th>\n'
    line = line + '<th>Exposure</th>\n'
    line = line + '<th>Net Counts</th>\n'
    line = line + '<th>Error</th>\n'
    line = line + '<th>DeadTime Correction</th>\n'
    line = line + '</tr>\n'
     
    out  = sorted_by_time(data)

    for ent in out:
        atemp = re.split('\s+', ent)
        line = line + '<tr>\n'
        line = line + '<td style="text-align:right">'  + str(atemp[0]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[1]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[2]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[3]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[4]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[5]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[6]) + '</td>\n'
        line = line + '</tr>\n'

    line = line + '</table>\n\n'

    lien = line + '</body>\n</html>\n'

    return line

#-----------------------------------------------------------------------------------------
#-- sorted_by_time: rearrange the data in time order                                    --
#-----------------------------------------------------------------------------------------

def sorted_by_time(data):
    """
    rearrange the data in time order
    input:  data    --- assume that the third entry is time in the form 
                        of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
    ouput:  out     --- time sorted data
    """
    tlist = []
    tdict = {}
    for ent in data:
        atemp = re.split('\s+', ent)
        tm = atemp[2]
        tm = tm.replace('-', '')
        tm = tm.replace('T', '')
        tm = tm.replace(':', '')
        tm = int(float(tm))
        tlist.append(tm)
        tdict[tm] = ent

    tlist.sort()

    out = []
    for ent in tlist:
        out.append(tdict[ent])

    return out

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def send_notification():

    text = '\nNew AR Lac Observations are added to:\n\n'
    text = text + 'https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_trendings/ArLac/arlac_vis_montior.html'
    text = text + '\n'

    with open(zspace, 'w') as fo:
        fo.write(text)

    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" vkashyap@cfa.harvard.edu'
    ####os.system(cmd)
    cmd = 'cat ' + zspace + ' |mailx -s "Subject: New AR Lac Observation" ' + admin
    os.system(cmd)

    mcf.rm_files(zspace)

#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    create_html_and_plot()
