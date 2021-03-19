#!/usr/bin/env /data/mta/Script/Python3.8/envs/ska3-shiny/bin/python

#########################################################################################
#                                                                                       #
#   hrc_dose_html_updates.py:   create  html data pages for a report                    #
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
import time
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
import mta_common_functions     as mcf

web_address = 'https://cxc.cfa.harvard.edu/contrib/cxchrc/HRC_Exposure/'

#------------------------------------------------------------------------------------------------
#--- hrc_dose_plot_exposure_stat: read hrc database, and plot history of exposure             ---
#------------------------------------------------------------------------------------------------

def hrc_dose_make_data_html():
    """
    read hrc database, and create html page: 
    input:  indir   --- data directory
            output  --- output
    output: <web_dir>/Sub_html/<inst>_<sec>.html; a table of stat results
    """
    for hrc in ('hrci', 'hrcs'):
#
#--- hrc i has 9 sections and hrc s has 10 sections
#
        if hrc == 'hrci':
            cstop = 9
            indir = stat_i_dir
        else:
            cstop = 10
            indir = stat_s_dir

        for sec in range(0, cstop):
            hrc_sec = hrc + '_' + str(sec)
#
#--- a trend data file contains following columns:
#---    year,month,mean_acc,std_acc,min_acc,min_apos, 
#---    max_acc,max_apos,asig1, asig2, asig3, mean_dff,std_dff,
#---    min_dff, min_dpos,max_dff,max_dpos, dsig1, dsig2, dsig3
#
            data = read_stat_data(indir, hrc_sec)
#
#--- create a HTML page to display histrical data
#
            create_html_page(data, hrc, sec, cstop)

#--------------------------------------------------------------------------------
#--  create_html_page: create HTML page to display HRC historical data       ----
#--------------------------------------------------------------------------------

def create_html_page(data, hrc, sec, ctop):
    """
    create HTML page to display HRC historical data.
    input:  data --- a list of lists of data:
                date,year,month,mean_acc,std_acc,min_acc,min_apos, 
                max_acc,max_apos,asig1, asig2, asig3, mean_dff,std_dff,
                min_dff, min_dpos,max_dff,max_dpos, dsig1, dsig2, dsig3
            hrc --- hrci or hrcs
            sec --- section
            ctop    --- max number of section
    output: <web_dir>/Sub_html/<inst>_<sec>.html
    """
#
#--- today's date
#
    today = time.strftime('%Y:%m:%d', time.gmtime())
    atemp = re.split(':', today)
    year  = int(float(atemp[0]))
    mon   = int(float(atemp[1]))
    day   = int(float(atemp[2]))

    smon = str(mon)
    if mon < 10:
        smon = '0' + smon

    sday = str(day)
    if day < 10:
        sday = '0' + sday
#
#--- start composing a html page
#
    aline  = read_template('sub_top1')    

    if hrc == 'hrci':
        hname = 'HRC I'
        wname = 'HRCI'
        inst  = 'i'
    else:
        hname = 'HRC S'
        wname = 'HRCS'
        inst  = 's'

    aline = aline + '<title>' + hname + ' Section: ' + str(sec) + ' History Data</title>\n'
    aline = aline + "</head>\n"

    aline = aline + '<body style="color:white;background-color:black">\n'
    aline = aline + '<div><a href="../hrc_exposure_map.html">Back to Top</a>\n'
    aline = aline + '<h2 style="text-align:center">Data: ' + hname 
    aline = aline + ' Section: ' + str(sec) +  '\n'

    if sec == 0:
        aline = aline + ' (<a href="./hrc' + inst + '_' + str(sec+1) + '.html">Next</a>)'
    elif sec == ctop-1:
        aline = aline + ' (<a href="./hrc' + inst + '_' + str(sec-1) + '.html">Prev</a>)'
    else:
        aline = aline + ' (<a href="./hrc' + inst + '_' + str(sec-1) + '.html">Prev</a> : '
        aline = aline + '<a href="./hrc' + inst + '_' + str(sec+1) + '.html">Next</a>)'

    aline = aline + '</h2>\n'

    aline = aline + "<div style='padding-bottom:30px'>\n"
    aline = aline + '<table border=1>\n'
    aline = aline + read_template('sub_col_header')
#
#--- open the data to indivisual lists
#
    [year,month,mean_acc,std_acc,min_acc,min_apos, \
     max_acc,max_apos,asig1, asig2, asig3, mean_dff,std_dff,\
     min_dff, min_dpos,max_dff,max_dpos,dsig1, dsig2, dsig3] = data

    dlen = len(year)

    for i in range(0, dlen):

        smonth = str(int(month[i]))
        if month[i] < 10:
            smonth = '0' + smonth

        syear = int(year[i])
#
#--- converting digit to letters, i.e. 1 to Jan
#
        cmonth = mcf.change_month_format(month[i])
#
#--- monthly HRC dose data
#
        if mean_dff[i] == 0 and std_dff[i] == 0:
            aline = aline + '<tr><td>%d</td><td>%d</td>' % (year[i], month[i])
            aline = aline + read_template('sub_nodata')
        else:
            aline = aline + '<tr><td>%d</td>' % year[i]
            aline = aline + '<td>%d</td>'     % month[i]

            try:
                aline = aline + '<td>%4.4f</td>'  % mean_dff[i]
                aline = aline + '<td>%4.4f</td>'  % std_dff[i]
                #aline = aline + '<td>%4.1f</td>'  % min_dff[i]
                #aline = aline + '<td>%s</td>'     % min_dpos[i]
                aline = aline + '<td>%4.1f</td>'  % max_dff[i]
                aline = aline + '<td>%s</td>'     % max_dpos[i]
                aline = aline + '<td>%4.1f</td>'  % dsig1[i]
                aline = aline + '<td>%s</td>'     % dsig2[i]
                aline = aline + '<td>%4.1f</td>\n'% dsig3[i]

#               if hrc == 'hrci':
#                   aline = aline + '<td><a href="' + data_i_dir + '/Month/'
#               else:
#                   aline = aline + '<td><a href="' + data_s_dir + '/Month/'
#               aline = aline +  wname + '_' + smonth + '_' + str(syear) + '.fits.gz">fits</a></td>\n'

                aline = aline + '<td><a href="' + web_address + 'Image/' + wname + '/Month/'
                aline = aline +  wname + '_' + smonth + '_' + str(syear) 
                aline = aline +  '_' + str(sec) + '.html">map</a></td>\n'
            except:
#
#--- for the case there is no data, print 'na'
#
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>'
                aline = aline + '<td>na</td>\n'
                aline = aline + '<td>na</td>\n'
#
#---- cumulative HRC dose data
#
        aline = aline + '<td>%4.4f</td>'    % mean_acc[i]
        aline = aline + '<td>%4.4f</td>'    % std_acc[i]
        #aline = aline + '<td>%4.1f</td>'    % min_acc[i]
        #aline = aline + '<td>%s</td>'       % min_apos[i]
        aline = aline + '<td>%4.1f</td>'    % max_acc[i]
        aline = aline + '<td>%s</td>'       % max_apos[i]
        aline = aline + '<td>%4.1f</td>'    % asig1[i]
        aline = aline + '<td>%s</td>'       % asig2[i]
        aline = aline + '<td>%4.1f</td>\n'  % asig3[i]

#        if hrc == 'hrci':
#            aline = aline + '<td><a href="' + data_i_dir + '/Cumulative/'
#        else:
#            aline = aline + '<td><a href="' + data_s_dir + '/Cumulative/'
#        aline = aline +  wname + '_' + smonth + '_' + str(syear) + '.fits.gz">fits</a></td>\n'

        aline = aline + '<td><a href="' + web_address + 'Image/' + wname + '/Cumulative/'
        aline = aline +  wname + '_08_1999_' + smonth + '_' + str(syear) + '_' + str(sec) +  '.html">map</a></td>\n'

#
#--- put header every new year so that we can read data easier
#
        if month[i] % 12 == 0 and i != (dlen-1):
            aline = aline + read_template('sub_col_header')

    aline = aline + '</table>\n\n'
    aline = aline + '<div style="padding_top:10px; padding_bottom:10px;">'
    aline = aline + '<a href="../hrc_exposure_map.html">Back to Top</a>\n'
#
#--- add today's date as update date
#
    today = time.strftime('%Y:%m:%d', time.gmtime())
    atemp = re.split(':', today)
    cmon  = mcf.change_month_format(int(float(atemp[1])))
    today = cmon + ' ' + atemp[2] + ', ' + atemp[0]

    aline = aline + read_template('sub_footer').replace('#UPDATE#', today)

    outdir = web_dir + 'Sub_html/' + hrc + '_' + str(sec) +  '.html'
    with open(outdir, 'w') as fo:
        fo.write(aline)

#------------------------------------------------------------------------------------
#-- update_main_html: update the top html page                                     --
#------------------------------------------------------------------------------------

def update_main_html():
    """
    update the top html page --- changes are just updated date
    input:  none, but read from <house_keeping>/exp_template
    output: <web_dir>/exposure.html
    """
    today = time.strftime('%m:%d:%Y', time.gmtime())
    atemp = re.split(':', today)
    tyear = int(float(atemp[2]))
    mon   = int(float(atemp[0]))
    day   = int(float(atemp[1]))
#
#--- display date
#
    cmon  = mcf.change_month_format(mon)
    today = cmon + ' ' + mcf.add_leading_zero(day) + ', ' + str(tyear)
#
#--- link date for the most recent plots
#
    [lyear, lmon] = find_last_entry_data()

    smon  = float(lmon)
    if smon < 10:
        lmon = '0' + lmon
    line = lmon + '_' + lyear

    data  = read_template('main_page')
    data  = data.replace('#LATEST#', line)
    data  = data.replace('#UPDATE#', today)

    ofile = web_dir + 'hrc_exposure_map.html'
    with open(ofile, 'w') as fo:
        fo.write(data)

#------------------------------------------------------------------------------------
#-- create_img_html: create htmls to display exposure map                          --
#------------------------------------------------------------------------------------

def create_img_html(year='', month=''):
    """
    create htmls to display exposure map.
    input:  year    --- year, if it is not given, the last month is used
            month   --- month, if it is not given, the last month is used
    output: <web_dir>/Image/HRC<inst>/Month/HRC<inst>_<mm>_<yyyy>_<sec>.html
            <web_dir>/Image/HRC<inst>/Cumulative/HRC<inst>_08_1999_<mm>_<yyyy>_<sec>.html
    """
#
#--- find the current year/month
#
    atemp  = time.strftime('%Y:%m', time.gmtime())
    [cyear, cmonth] = re.split(':', atemp)
    cyear  = int(float(cyear))
    cmonth = int(float(cmonth))
    chk = 0
#
#--- if the year and month were not passed, set them to those of the last month
#
    if mcf.is_neumeric(year):
        year  = int(float(year))
        month = int(float(month))
        if year < cyear:
            chk   = 1
    else:
        year  = cyear
        month = cmonth -1
        if month < 1:
            month = 12
            year -= 1
#
#--- set one month before and one month after
#
    pyear  = year
    pmonth = month - 1
    if pmonth < 1:
        pmonth =12 
        pyear -= 1

    nyear  = year
    nmonth = month + 1
    if nmonth > 12:
        nmonth = 1
        nyear += 1
#
#--- set link dates
#
    ldate  = mcf.add_leading_zero(month)  + '_' + str(year)
    pdate  = mcf.add_leading_zero(pmonth) + '_' + str(pyear)
    ndate  = mcf.add_leading_zero(nmonth) + '_' + str(nyear)

    monthly = read_template('mon_img_page')
    for  inst in ('S', 'I'):
        if inst == 'S':
            cstop = 10
            odir  = web_img_s_dir + 'Month/'
        else:
            cstop = 9
            odir  = web_img_i_dir + 'Month/'

        for sec in range(0, cstop):
#
#--- create image file link
#
            png = 'HRC' + inst + '_' + ldate + '_' + str(sec) + '.png'
            cpath = web_dir + 'Image/HRC' + inst + '/Month/' + png
            if os.path.isfile(cpath):
                pnglink = './' + png
            else:
                pnglink = './no_data.png'
#
#--- create link paths
#
            pfile  = './HRC' + inst + '_' + pdate + '_' + str(sec)   + '.html'
            nfile  = './HRC' + inst + '_' + ndate + '_' + str(sec)   + '.html'
            psfile = './HRC' + inst + '_' + ldate + '_' + str(sec-1) + '.html'
            nsfile = './HRC' + inst + '_' + ldate + '_' + str(sec+1) + '.html'
#
#--- section link
#
            if sec == 0:
                seclink = '<a href="' + nsfile + '">Next Section</a><br />'
            elif sec == cstop-1:
                seclink = '<a href="' + psfile + '">Prev Section</a><br />'
            else:
                seclink = '<a href="' + psfile + '">Prev Section</a>  '
                seclink = seclink + '<a href="' + nsfile + '">next Section</a><br />'
#
#--- time order link
#
            if year == 1999 and month == 8:
                tolink = '<a href="' + nfile + '">Next Month</a><br />'
            elif chk == 0:
                tolink = '<a href="' + pfile + '">Prev Month</a><br /> '
            else:
                tolink = '<a href="' + pfile + '">Prev Month</a>  '
                tolink = tolink + '<a href="' + nfile + '">Next Month</a><br />'
#
#--- section main link
#
            sublink = '../../../Sub_html/hrc' + inst.lower() + '_' + str(sec) + '.html'
#
#--- cumulative link
#
            cumlink = '<a href="../Cumulative/HRC' + inst + '_08_1999_' + ldate 
            cumlink = cumlink + '_' + str(sec) + '.html">Cumulative Plot</a>' 
#
#--- replace texts in the template
#
            otemp = monthly
            otemp = otemp.replace("#YEAR#",    str(year))
            otemp = otemp.replace("#MONTH#",   mcf.add_leading_zero(month))
            otemp = otemp.replace("#INST#",    inst)
            otemp = otemp.replace("#SEC#",     str(sec))
            otemp = otemp.replace("#PNGLINK#", pnglink)
            otemp = otemp.replace("#LATEST#",  ldate)
            otemp = otemp.replace("#CUMLINK#", cumlink)
            otemp = otemp.replace("#SECLINK#", seclink)
            otemp = otemp.replace("#TOLINK#",  tolink)
            otemp = otemp.replace("#SUBLINK#", sublink)
#
#--- set output fine name
#
            ofile = odir + 'HRC' + inst + '_' + ldate + '_' + str(sec) + '.html'
            with open(ofile, 'w') as fo:
                fo.write(otemp)
#
#--- cumulatvie page
#
    cumulative  = read_template('cum_img_page')
    for  inst in ('S', 'I'):
        if inst == 'S':
            cstop = 10
            odir  = web_img_s_dir + 'Cumulative/'
        else:
            cstop = 9
            odir  = web_img_i_dir + 'Cumulative/'

        for sec in range(0, cstop):
#
#--- create link paths
#
            pfile  = './HRC' + inst + '_08_1999_' + pdate + '_' + str(sec)   + '.html'
            nfile  = './HRC' + inst + '_08_1999_' + ndate + '_' + str(sec)   + '.html'
            psfile = './HRC' + inst + '_08_1999_' + ldate + '_' + str(sec-1) + '.html'
            nsfile = './HRC' + inst + '_08_1999_' + ldate + '_' + str(sec+1) + '.html'
#
#--- section link
#
            if sec == 0:
                seclink = '<a href="' + nsfile + '">next Section</a><br />'
            elif sec == cstop-1:
                seclink = '<a href="' + psfile + '">Prev Section</a><br />'
            else:
                seclink = '<a href="' + psfile + '">Prev Section</a>  '
                seclink = seclink + '<a href="' + nsfile + '">next Section</a><br />'
#
#--- time order link
#
            if year == 1999 and month == 8:
                tolink = '<a href="' + nfile + '">Next Month</a><br />'
            elif chk == 0:
                tolink = '<a href="' + pfile + '">Prev Month</a><br /> '
            else:
                tolink = '<a href="' + pfile + '">Prev Month</a>  '
                tolink = tolink + '<a href="' + nfile + '">Next Month</a><br />'
#
#--- section main link
#
            sublink = '../../../Sub_html/hrc' + inst.lower() + '_' + str(sec) + '.html'
#
#--- month link
#
            monlink = '<a href="../Month/HRC' + inst + '_' + ldate + '_' + str(sec) + '.html">'
            monlink = monlink + 'Month Plot</a>' 
#
#--- replace texts in the template
#
            otemp = cumulative
            otemp = otemp.replace("#YEAR#",    str(year))
            otemp = otemp.replace("#MONTH#",   mcf.add_leading_zero(month))
            otemp = otemp.replace("#INST#",    inst)
            otemp = otemp.replace("#SEC#",     str(sec))
            otemp = otemp.replace("#LATEST#",  ldate)
            otemp = otemp.replace("#MONLINK#", monlink)
            otemp = otemp.replace("#SECLINK#", seclink)
            otemp = otemp.replace("#TOLINK#",  tolink)
            otemp = otemp.replace("#SUBLINK#", sublink)
#
#--- set output fine name
#
            ofile = odir + 'HRC' + inst + '_08_1999_' + ldate + '_' + str(sec) + '.html'
            with open(ofile, 'w') as fo:
                fo.write(otemp)


#------------------------------------------------------------------------------------
#-- read_stat_data: read data from acis/hrc history data files                    ---
#------------------------------------------------------------------------------------

def read_stat_data(indir, inst):
    """
    read data from acis/hrc history data files
    input:  indir--- directory where the data locate
            inst --- instruments hrci_<section>, hrcs_<section>
    output: a list of lists
            [year,month,mean_acc,std_acc,min_acc,min_apos, max_acc,\
            max_apos,asig1, asig2, asig3, mean_dff,std_dff,min_dff, \
            min_dpos,max_dff,max_dpos,dsig1, dsig2, dsig3]
    """
    ifile = indir +  inst + '_' + 'acc_out'
    data  = mcf.read_data_file(ifile)
    save  = convert_to_columndata(data)

    ifile = indir +  inst + '_' + 'dff_out'
    data  = mcf.read_data_file(ifile)
    save2 = convert_to_columndata(data)
# 
#--- odata contains:
#---  [year,month,mean_acc,std_acc,min_acc,min_apos, max_acc,\
#---   max_apos,asig1, asig2, asig3, mean_dff,std_dff,min_dff, \
#---   min_dpos,max_dff,max_dpos,dsig1, dsig2, dsig3]
#
    odata = save + save2[2:]   #--- skipping date part from the second list
    return odata

#------------------------------------------------------------------------------------
#-- convert_to_columndata: convert a list of data into a list of lists            ---
#------------------------------------------------------------------------------------

def convert_to_columndata(data):
    """
    convert a list of data into a list of lists
    input:  data    --- a list of data
    output: save    --- a list of lists of data
    """

    slen = len(re.split('\s+', data[0]))
    save = []
    for k in range(0, slen):
        save.append([])

    for ent in data:
        out = re.split('\s+', ent)
        for k in range(0, slen):
            try:
                save[k].append(float(out[k]))
            except:
                save[k].append(out[k])


    return save

#------------------------------------------------------------------------------------
#-- read_template: read a template                                                 --
#------------------------------------------------------------------------------------

def read_template(part):
    """
    read a template
    input:  part    --- a file name which contain the template
    output: out     --- a text
    """
    ifile = house_keeping + 'Templates/' + part
    with open(ifile, 'r') as f:
        out = f.read()

    return out

#------------------------------------------------------------------------------------
#-- find_last_entry_data: find the last entry date from a data file                --
#------------------------------------------------------------------------------------

def find_last_entry_data():
    """
    find the last entry date from a data file
    input: none, but read from <stat_i_dir>/hrci_4_acc_out
    output: [<year>, <month>]
    """
    ifile = stat_i_dir + 'hrci_4_acc_out'
    data  = mcf.read_data_file(ifile)

    atemp = re.split('\s+', data[-1])

    return [atemp[0], atemp[1]]


#------------------------------------------------------------------------------------

if __name__ == '__main__':

    update_main_html()

    if len(sys.argv) > 2:
        year = int(float(sys.argv[1]))
        mon  = int(float(sys.argv[2]))
        create_img_html(year, mon)
    else:
        hrc_dose_make_data_html()
        update_main_html()
        create_img_html('','')

#    for year in range(1999, 2020):
#        for mon in range(1, 13):
#            if year == 1999 and mon < 8:
#                continue
#            print("TIME: " + str(year) + ' : ' + str(mon))
#            create_img_html(year, mon)



