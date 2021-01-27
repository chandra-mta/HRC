#!/usr/bin/env /data/mta/Script/Python3.9/bin/python3

#############################################################################################
#                                                                                           #
#   extract_hz43_stat.py: find out hz43 observations and update the data tables             #
#                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                           #
#           Last Update: Jan 21, 2021                                                       #
#                                                                                           #
#############################################################################################

import sys
import os
import string
import re
import math
import unittest
import time
import random
import numpy
import astropy.io.fits  as pyfits
from datetime import datetime
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
sys.path.append(hrc_common)

import mta_common_functions as mcf
import filter_evt_file      as fef
import hrc_common_functions as hcf
import adjust_hz43_position as ahp
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

NULL   = 'NULL'

outlist   = ['samp_p_list', 'samp_n_list', 'pi_p_list',  'pi_n_list']

#-----------------------------------------------------------------------------------------
#-- extract_hz43_stat: extract hz43 data, analyze and append the results to data table  --
#-----------------------------------------------------------------------------------------

def extract_hz43_stat(obsid, evt1):
    """
    extract hz43 data, analyze and append the results to data table
    input:  obsid   --- obsid
            evt1    --- evnt1 file name with full path
    output: see pinrt_out_<1/2>
    """
    parea = 125663.7         #---- 200 px radius area
#
#--- create the saving directory if it does not exist
#
    sdir = data_dir + 'Fittings/' + str(obsid)
    if not os.path.isdir(sdir):
        cmd = 'mkdir ' + sdir
        os.system(cmd)

    out = extract_data(obsid, evt1)
    if out == False:
        print(str(obsid) + " main data analysis stopped")

        fout = house_keeping + 'skip_obsids'
        with open(fout, 'a') as fo:
            fo.write(str(obsid) + '\n')

        cmd = 'rm *.fits*'
        os.system(cmd)

        return False
#
#---    out list below contains:
#---    hdata: [inst, date, tstart, tstop, exposure]
#---    others: [tg_m, tg_r, tg_d, pi, amp_sf, sumamps, samp] all of them are lists
#
    [hdata, cent_data, cent_bkg, arms_data, arms_bkg] = out

    inst  = hdata[0].lower()                        #--- hrc-i or hrc-s
    tspan = hdata[-1]                               #--- exposure time
#
#--- write down the center part results: [avg, sig, med, crate, brate, tspan]
#
    if (len(cent_data) > 0) and (len(cent_data[3]) >=50):

        pi_stats = get_bkg_subtracted_stats(cent_data[3], cent_bkg[3], tspan, 1, 1, cent=1)
        sp_stats = get_bkg_subtracted_stats(cent_data[6], cent_bkg[6], tspan, 1, 1, cent=1)

        save_ind_dist_center(sdir, cent_data[3], cent_data[6])

        oname = 'pi_center_list'
        print_out(obsid, hdata, pi_stats, oname)
    
        oname = 'samp_center_list'
        print_out(obsid, hdata, sp_stats, oname)
    else:
        fout = house_keeping + 'skip_obsids'
        with open(fout, 'a') as fo:
            fo.write(str(obsid) + '\n')

        cmd = 'rm *.fits*'
        os.system(cmd)

        return False
#
#---    bin_data below contains:
#---    [samp_p_list, samp_n_list, pi_p_list, pi_n_list, 
#---     bsamp_p_list, bsamp_n_list, bpi_p_list, bpi_n_list,
#---     area_p_list, area_n_list, barea_p_list, barea_n_list]
#---    all of them are lists of lits (about 20 lists in each list)
#
    bin_data    = separate_into_bin(arms_data, arms_bkg, inst)
#
#--- save the distributions
#
    save_ind_dist(sdir, bin_data[0], bin_data[1], bin_data[2], bin_data[3])
#
#--- write down the binned results to files
#--- outlist   = ['samp_p_list', 'samp_n_list', 'pi_p_list',  'pi_n_list']
#
    for k in range(0, 4):
        odata = bin_data[k]     #---- main data
        bdata = bin_data[k+4]   #---- bakcground data
        if k == 0 or k == 2:
            area  = bin_data[-4]
            barea = bin_data[-2]
        else:
            area  = bin_data[-3]
            barea = bin_data[-1]

        for m in range(0, len(odata)):
            if len(odata[m]) < 50:
                continue 

            stats =  get_bkg_subtracted_stats(odata[m], bdata[m], tspan, area[m], barea[m])
            oname = outlist[k] + '_' + str(m) 
            print_out(obsid, hdata, stats, oname)
#
#--- clean up 
#
    cmd = 'rm *.fits*'
    os.system(cmd)

    return True

#-----------------------------------------------------------------------------------------
#-- print_out: print out the center data                                              --
#-----------------------------------------------------------------------------------------

def print_out(obsid, head, stats, oname):
    """
    print out the center data
    input:  obsid   --- obsid
            heads   --- a list of the data; the first three are used;
                        [inst, date, tstart, tstop] 
            stats   --- stat resuts: [avg, sig, med, crate, brate, tspan]
            oname   --- head part of the output file name
    output: <data_dir>/samp_center_list_<i/s> pi_center_list_<i/s>
    """
    [inst, date, tstart, tstop, exposure]            = head
    [avg, sig, med, crate, brate, tspan, cerr, berr] = stats

    line = str(tstart) + '\t' + str(obsid)  + '\t' + str(date) + '\t' 
    line = line + str(int(exposure))        + '\t'
    line = line + "%3.3f" % round(avg, 3)   + '\t' 
    line = line + "%3.3f" % round(sig, 3)   + '\t' 
    line = line + str(int(med))             + '\t' 
    line = line + "%3.8f" % round(crate, 8) + '\t' 
    line = line + "%3.8f" % round(brate, 8) + '\t'
    line = line + "%3.8f" % round(cerr, 8)  + '\t' 
    line = line + "%3.8f" % round(berr, 8)  + '\n'

    if inst.lower() == 'hrc-i':
        outfile = data_dir + oname + '_i'
    else:
        outfile = data_dir + oname + '_s'

    with open(outfile, 'a') as fo:
        fo.write(line)

#-----------------------------------------------------------------------------------------
#-- extract_data: extract needed information for a given observation                    --
#-----------------------------------------------------------------------------------------

def extract_data(obsid, evt1):
    """
    extract needed information for a given observation 
    input:  obsid   --- obsid
            evt1    --- hrc evt1 fits file name
    output: a list of lists:
                    [inst, date, start,  crsv, crsu, amp_sf, pha, pi, sumamps, samp, chip_id, 
                     tg_m, tg_mlam, tg_srcid, tg_part, tg_smap]
            note: inst,date, start are strings, but all others are lists.
    """
#
#--- evt1 data --- already extracted
#
    if evt1 == "":
        print("no EVT1 fits file")
        return False
#
#--- extract evt1.5 data
#
    evta = hcf.run_arc5gl(0, 0,  obsid = obsid, operation='retrieve', level ='1.5', filetype='tgevt1')
    if evta == "":
        print("no EVT1.5 fits file")
        sfile = house_keeping + 'skip_obsids'
        with open(sfile, 'a') as fo:
            fo.write(str(obsid) + '\n')

        return False
#
#--- find which instrument
#
    hdr   = pyfits.getheader(evt1, 1)
    inst  = hdr['detnam']
    date  = hdr['date-obs']
    start = hdr['tstart']
    stop  = hdr['tstop']
    expo  = hdr['exposure']

    hdata = [inst, date, start, stop, expo]
#
#--- we need to drop a large off axis data. set 9 arcmin
#
    ra_pnt   = float(hdr['ra_pnt'])
    dec_pnt  = float(hdr['dec_pnt'])
    ra_targ  = float(hdr['ra_targ'])
    dec_targ = float(hdr['dec_targ'])

    adiff    = ra_targ - ra_pnt
    bdiff    = dec_targ - dec_pnt
    gout     = adiff*adiff + bdiff*bdiff
    adiff    = math.sqrt(gout)
    if adiff > 0.15:
        print("Offset too large (>9 amin); dropped")
        return False
#
#--- combine the data
#
    comb_evt1 = 'comb_evt.fits'
    col_list  = ['amp_sf', 'sumamps']
    add_coldata_to_fits(evta, evt1, col_list, comb_evt1)
#
#--- clean the data
#
    clean1 = 'evt1_clean.fits'
    carea  = 'clean_cent.fits'
    barea  = 'clean_cent_bkg.fits'
    fef.filter_evt_file(comb_evt1, clean1, obsid)
#
#--- get the center part
#
    chk = get_center_data(obsid, evt1, clean1, carea, barea)

    if chk:
        cent_data = get_data(carea, ['tg_m', 'tg_r', 'tg_d', 'pi','amp_sf', 'sumamps'])
        samp      = compute_samp(cent_data[5], cent_data[4], inst)
        cent_data.append(samp)

        cent_bkg  = get_data(barea, ['tg_m', 'tg_r', 'tg_d', 'pi','amp_sf', 'sumamps'])
        samp      = compute_samp(cent_bkg[5], cent_bkg[4], inst)
        cent_bkg.append(samp)
    else:
        cent_data = []
        cent_bkg  = []
#
#--- get the arm parts
#
    carea = 'clean_arms.fits'
    barea = 'clean_arms_bkg.fits'
    get_arm_data(clean1, carea, barea)

    arms_data = get_data(carea, ['tg_m', 'tg_r', 'tg_d', 'pi','amp_sf', 'sumamps'])
    samp      = compute_samp(arms_data[5], arms_data[4], inst)
    arms_data.append(samp)

    arms_bkg  = get_data(barea, ['tg_m', 'tg_r', 'tg_d', 'pi','amp_sf', 'sumamps'])
    samp      = compute_samp(arms_bkg[5], arms_bkg[4], inst)
    arms_bkg.append(samp)

    return [hdata, cent_data, cent_bkg, arms_data, arms_bkg]

#-----------------------------------------------------------------------------------------
#-- get_center_data: extract center part of the data and the background                 --
#-----------------------------------------------------------------------------------------

def get_center_data(obsid, evt1, clean, carea, barea):
    """
    extract center part of the data and the background
    input:  obsid   --- obsid
            evt1    --- event 1 data fits file name
            clean   --- cleaned fits file name
            carea   --- extracted main area output fits file
            barea   --- txtracted background output  fits file
    output: carea /barea
    """
#
#--- get det coordinates computed previously
#
    [detx, dety] = find_det_coordinates(evt1)
#
#--- try tgdetect to find the center part first
#
    dchk = 0
    try:
        cmd = ' tgdetect ' + evt1 + ' none outfile=zinfo.fits clobber=yes'
        os.system(cmd)

        if not os.path.isfile('./zinfo.fits'):
            return False

    
        t     = pyfits.open('zinfo.fits')
        tdata = t[1].data
        t.close()
        pos = 0
        test = tdata.field('x')[0]
        dchk = 1
#
    except:
        sfile = house_keeping + 'skip_obsids'
        with open(sfile, 'a') as fo:
            fo.write(str(obsid) + '\n')
        return False
#
#--- get the information about the center source area
#
    if dchk == 1:
        psave = []
        for ent in ('x', 'y', 'net_counts', 'net_counts_err', 'bkg_counts', \
                    'bkg_counts_err', 'net_rate', 'net_rate_err', 'bkg_rate',\
                    'bkg_rate_err'):
            try:
                data  = tdata.field(ent)
                out   = data[pos]
                psave.append(out)
            except:
                if detx == '':
                    return False

    mcf.rm_files('zinfo.fits')
#
#--- center part
#
    if detx != '':
        cmd   = 'dmcopy "' + clean + '[(detx,dety)=circle(' + str(detx) + ',' + str(dety) 
    else:
        cmd   = 'dmcopy "' + clean + '[(detx,dety)=circle(' + str(psave[0]) + ',' + str(psave[1]) 
    cmd   = cmd + ',200)]" outfile=' + carea + ' clobber=yes'
    os.system(cmd)
#
#--- background area; 400 pix above the center area
#
    ypos  = float(psave[1]) + 400
    cmd   = 'dmcopy "' + clean + '[(detx,dety)=circle(' + str(psave[0]) + ',' + str(ypos) 
    cmd   = cmd + ', 200)]" outfile=' + barea + ' clobber=yes'
    os.system(cmd)
    
    line  = str(psave[0])
    for  k in range(1, len(psave)):
        line  = line + '\t'  + str(psave[k])

    line = line + '\n'

    ofile = data_dir + 'Fittings/' + str(obsid) + '/center_info'
    with open(ofile, 'w') as fo:
        fo.write(line)

    return True

#-----------------------------------------------------------------------------------------
#-- find_det_coordinates: check whether manually determined det cooridates are avaliable -
#-----------------------------------------------------------------------------------------

def find_det_coordinates(evt1):
    """
    check whether manually determined det cooridates are avaliable
    input:  evt1            --- evt1 fits file name
    output: [detx, dety]    --- det coordinates
    """
    fout      = pyfits.open(evt1)
    fhead     = fout[1].header
    otime     = fhead['DATE-OBS']
    [ra, dec] = ahp.adjust_hz43_position(otime)
#
#--- convert coordinates from cel to det
#
    cmd = 'dmcoords ' + evt1 + ' opt=cel ra=' + str(ra) + ' dec=' + str(dec)
    cmd = cmd + ' verbose=1 > ' + zspace
    os.system(cmd)
#
#--- extract det coordindats
#
    detx = 0.0
    dety = 0.0
    info = mcf.read_data_file(zspace, remove=1)
    for ent in info:
        mc = re.search('DETX,DETY', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            detx = float(atemp[1])
            dety = float(atemp[2])
            break

    return [detx, dety]

#-----------------------------------------------------------------------------------------
#-- get_arm_data: run select_letg_arm  to get the letg covered area and background area  -
#-----------------------------------------------------------------------------------------

def get_arm_data(ofits, out1, out2):
    """
    run select_letg_arm  to get the letg covered area and background area
    input:  ofits   --- input fits file name
            out1    --- output file name of letg covered area
            out2    --- outpuf file name of background area
    output: out1 and out2
    """
#
#--- the data extraction of letg covered area
#
    select_letg_arm(ofits, out1)
#
#--- the data extraction of the background area
#
    select_letg_arm(ofits, out2, include=0)

#-----------------------------------------------------------------------------------------
#-- separate_into_bin: separate data into specified bin interval                       ---
#-----------------------------------------------------------------------------------------

def separate_into_bin(data, bkg, inst=''):
    """
    separate data into specified bin interval
    input:  data    --- a list of lists of:
                        [tg_m, tg_r, tg_d, pi, amp_sf, sumamps, samp]
            bkg     --- a list of background of above
            inst    --- hrc-i or hrc-s
    output: samp_p_list --- binned data of samp in positive side
            samp_n_list --- binned data of samp in negative side
            pi_p_list   --- binned data of pi in positive side
            pi_n_list   --- binned data of pi in negative side
    Note:
        tg_d in degree
        1 pixel = 4.265e-5 degrees
        1 tap   = 256 pixels = 1.88893A
        this create the interval to be: 0.058271186 in tg_d
        1.592e-04 deg/pix for hrc letg
        0.5731 arcsec/pix

    """
    step = 0.057802036      #--- 10A step in tg_d
    factor = 18667.0        #--- 1 tg_d x 1tg_d to 1 pix x 1pix area conversion
    d_p_p  = 1.592e-04
    factor = 1./d_p_p**2    #--- 1 tg_d x 1tg_d to 1 pix x 1pix area conversion

    cut  = 19
    if inst == 'hrc-i':         #--- if the instrument is HRC I, 20A step
            cut = 9
#
#--- separate arm data
#
    [tg_m, tg_r, tg_d, pi, amp_sf, sumamps, samp] = data

    samp_p_list  = []           #--- samp positive side of data
    samp_n_list  = []           #--- samp negative side of data
    pi_p_list    = []           #--- pi positive side of data
    pi_n_list    = []           #--- pi netative side of data
    tg_p_list    = []
    tg_n_list    = []
#
#--- all of them have 20 sub lists
#
    for k in range(0, 20):
        samp_p_list.append([])
        samp_n_list.append([])
        pi_p_list.append([])
        pi_n_list.append([])
        tg_p_list.append([])
        tg_n_list.append([])
#
#--- tg_r gives the value along the dispersion axis direction
#
    dcnt = 0
    for k in range(0, len(tg_m)):
            try:
                val = float(tg_r[k])
            except:
                continue
            if str(val) == 'nan':
                continue

            aval = abs(val)
            pos = int(aval / step)
            if pos < 5:
                continue
            if pos > cut:
                continue

            if inst == 'hrc-i':
                pos = int(0.5 * pos)

            if val>= 0:
                try:
                    samp_p_list[pos].append(samp[k])
                    pi_p_list[pos].append(pi[k])
                    tg_p_list[pos].append(tg_d[k])
                except:
                    pass

            else:
                try:
                    samp_n_list[pos].append(samp[k])
                    pi_n_list[pos].append(pi[k])
                    tg_n_list[pos].append(tg_d[k])
                except:
                    pass
#
#--- get the heights of each bin so that we can use it to get the
#--- same area size for the backgruond area
#
    height_p_list = []
    area_p_list   = []
    for k in range(0, len(tg_p_list)):
        try:
            dmin = min(tg_p_list[k])
            dmax = max(tg_p_list[k])
            h    = dmax - dmin
        except:
            h    = 0.0
        height_p_list.append(h)
        area_p_list.append(factor*h*step)
    height_n_list = []
    area_n_list   = []
    for k in range(0, len(tg_n_list)):
        try:
            dmin = min(tg_n_list[k])
            dmax = max(tg_n_list[k])
            h    = dmax - dmin
        except:
            h    = 0.0
        height_n_list.append(h)
        area_n_list.append(factor*h*step)
#
#--- separate background data
#
    [tg_m, tg_r, tg_d, pi, amp_sf, sumamps, samp] = bkg

    bsamp_p_list  = []           #--- samp positive side of data
    bsamp_n_list  = []           #--- samp negative side of data
    bpi_p_list    = []           #--- pi positive side of data
    bpi_n_list    = []           #--- pi netative side of data
    btg_p_list    = []
    btg_n_list    = []
#
#--- all of them have 20 sub lists
#
    for k in range(0, 20):
        bsamp_p_list.append([])
        bsamp_n_list.append([])
        bpi_p_list.append([])
        bpi_n_list.append([])
        btg_p_list.append([])
        btg_n_list.append([])
#
#--- tg_r gives the value along the dispersion axis direction
#
    for k in range(0, len(tg_m)):
            try:
                val = float(tg_r[k])
            except:
                continue
            if str(val) == 'nan':
                continue

            aval = abs(val)
            pos = int(aval / step)
            if pos < 5:
                continue
            if pos > cut:
                continue
            if inst == 'hrc-i':
                pos = int(0.5 * pos)
#
            if val >= 0:
                try:
                    bsamp_p_list[pos].append(samp[k])
                    bpi_p_list[pos].append(pi[k])
                    btg_p_list[pos].append(tg_d[k])
                except:
                    pass

            else:
                try:
                    bsamp_n_list[pos].append(samp[k])
                    bpi_n_list[pos].append(pi[k])
                    btg_n_list[pos].append(tg_d[k])
                except:
                    pass
#
#--- get the heights of each bin 
#
    bheight_p_list = []
    barea_p_list   = []
    for k in range(0, len(btg_p_list)):
        try:
            dmin = min(btg_p_list[k])
            dmax = max(btg_p_list[k])
            h    = dmax - dmin
        except:
            h    = 0.0
        bheight_p_list.append(h)
        barea_p_list.append(factor*h*step)
    bheight_n_list = []
    barea_n_list   = []
    for k in range(0, len(tg_n_list)):
        try:
            dmin = min(btg_n_list[k])
            dmax = max(btg_n_list[k])
            h    = dmax - dmin
        except:
            h    = 0.0
        bheight_n_list.append(h)
        barea_n_list.append(factor*h*step)

    out = [samp_p_list, samp_n_list, pi_p_list, pi_n_list]
    out = out + [bsamp_p_list, bsamp_n_list, bpi_p_list, bpi_n_list]
    out = out + [area_p_list, area_n_list, barea_p_list, barea_n_list]
    return  out

#-----------------------------------------------------------------------------------------
#-- save_ind_dist: save distribution data for each bin of given obsid                   --
#-----------------------------------------------------------------------------------------

def save_ind_dist(sdir, samp_p_list, samp_n_list, pi_p_list, pi_n_list):
    """
    save distribution data for each bin of given obsid
    input:  sdir    --- the directroy path to save the data 
            samp_p_list --- positive side samp data
            samp_n_list --- negative side samp data
            pi_p_list   --- positive side pi data
            pi_n_list   --- negative side pi data
    output: <sdir>/samp_p_list<#> etc
    """
    for k in range(0, 18):
        try:
            odata  = samp_p_list[k]
            odata2 = pi_p_list[k]

            oname  = sdir + '/samp_p_list_' + str(k)
            oname2 = sdir + '/pi_p_list_'   + str(k)

            print_out_dist(odata,  oname)
            print_out_dist(odata2, oname2)
        except:
            pass

        try:
            odata  = samp_n_list[k]
            odata2 = pi_n_list[k]

            oname = sdir + '/samp_n_list_' + str(k)
            oname2= sdir + '/pi_n_list_'   + str(k)

            print_out_dist(odata,  oname)
            print_out_dist(odata2, oname2)
        except:
            pass

#-----------------------------------------------------------------------------------------
#-- save_ind_dist_center: save distribution data for each bin of given obsid           ---
#-----------------------------------------------------------------------------------------

def save_ind_dist_center(sdir, pi_list, samp_list):
    """
    save distribution data for each bin of given obsid
    input:  sdir    --- the directroy path to save the data 
            samp_p_list --- positive side samp data
            samp_n_list --- negative side samp data
            pi_p_list   --- positive side pi data
            pi_n_list   --- negative side pi data
    output: <sdir>/samp_p_list<#> etc
    """
    try:
        oname = sdir + '/samp_center_list'
        print_out_dist(samp_list, oname)
    except:
        pass

    try:
        oname = sdir + '/pi_center_list'
        print_out_dist(pi_list, oname)
    except:
        pass

#-----------------------------------------------------------------------------------------
#-- print_out_dist: print out distribution data                                       ----
#-----------------------------------------------------------------------------------------

def print_out_dist(odata, oname):
    """
    print out distribution data
    input:  odata   --- data (one dimension)
            oname   --- output file name
    """
    if len(odata) < 10:
        return False

    with open(oname, 'w') as fo:
        for ent in odata:
            fo.write(str(ent)+ '\n')

    return True

#-----------------------------------------------------------------------------------------
#-- create_stat_table: create a table of stat results for a given data sets             --
#-----------------------------------------------------------------------------------------

def create_stat_table(dlist):
    """
    create a table of stat results for a given data sets
    input:  dlist   --- a list of lists of data
    output: a list of lists of:
            avg     --- mean of the data
            std     --- standard error of the data
            med     --- median of the data
            cnt     --- the number of the data used
            if the data does not create proper stats, return  -999
    """
    avg = []
    std = []
    med = []
    cnt = []
    dlen = len(dlist)
    for k in range(0, dlen):
        stats = get_stat(dlist[k])
#
#--- if the statistics do not give proper values, use -999
#
        if stats == False:
            avg.append(-999)
            std.append(-999)
            med.append(-999)
            cnt.append(-999)
        else:
            avg.append(stats[0])
            std.append(stats[1])
            med.append(stats[2])
            cnt.append(stats[3])

    return [avg, std, med, cnt]

#-----------------------------------------------------------------------------------------
#-- get_stat: find a basic stat of the data                                            ---
#-----------------------------------------------------------------------------------------

def get_stat(data):
    """
    find a basic stat of the data
    input:  data    --- a list of data
    output: avg     --- mean of the data
            std     --- standard error of the data
            med     --- median of the data
            cnt     --- the number of the data used
            if the data does not create proper stats, return False
    """
    if len(data) > 0:
        avg = numpy.mean(data)
        std = numpy.std(data)
        med = numpy.median(data)
        cnt = len(data)
        return [avg, std, med, cnt]
    else:
        return False
    

#-----------------------------------------------------------------------------------------
#-- compute_samp: compute scaled samp                                                  ---
#-----------------------------------------------------------------------------------------

def compute_samp(sumamps, amp_sf, inst):
    """
    compute scaled samp
    input:  sumamp  --- a list of sumamp
            amp_sf  --- a list of amp_sf
            inst    --- instrument; either hrc-i or hrc-s
    output: samp    --- a list of scaled samp
    """ 
    
    if inst.lower() == 'hrc-i':
        const = 128
    else:
        const = 148

    samp = []
    for k in range(0, len(sumamps)):
        val = sumamps[k] * 2**(amp_sf[k] - 1) / const
        samp.append(val)

    return samp

#-----------------------------------------------------------------------------------------
#-- get_data: extract data from a fits file                                            ---
#-----------------------------------------------------------------------------------------

def get_data(fits, col_list):
    """
    extract data from a fits file
    input:  fits        --- fits file name
            col_list    --- a list of column names which you want to extract data
    output: save        --- a list of lists of data
    """

    hd    = pyfits.open(fits)
    tdata = hd[1].data

    save = []
    for col in col_list:
        #exec("save.append(tdata['%s'])" % (col))
        save.append(tdata[col])

    return save

#-----------------------------------------------------------------------------------------
#-- add_coldata_to_fits: adding data from the second fits file to the first fits file   --
#-----------------------------------------------------------------------------------------

def add_coldata_to_fits(ofits1, ofits2, col_names, outfile):
    """
    adding data from the second fits file to the first fits file
    input:  ofits1      --- the first fits file
            ofits2      --- the second fits file
            col_names   --- a list of column names which you want to copy from the 
                            second fits file
            outfile     --- output fits file name
    output: outfile     --- resulted fits file
    """
    otable = pyfits.open(ofits1)[1].data
    ocols  = otable.columns
    atable = pyfits.open(ofits2)[1].data
    acols  = atable.columns
#
#--- extract informaiton of the columns and the column data which are not overlapped with
#--- the data to be added
#
    out        = ocols.names
    ocol_names = numpy.setdiff1d(out, col_names)
    save = []
    for col in ocol_names:
        cent = ocols[col]
        data = otable[col]
        ndat = pyfits.Column(name=cent.name, format=cent.format, unit=cent.unit, array=data)
        save.append(ndat)
#
#--- set the column data
#
    oldcols  = pyfits.ColDefs(save)
#
#--- extreact information of the columns and the column data to be added from the
#--- second fits file
#
    save = []
    for col in col_names:
        cent = acols[col]
        data = atable[col]
        ndat = pyfits.Column(name=cent.name, format=cent.format, unit=cent.unit, array=data)
        save.append(ndat)
#
#--- define the column data
#
    newcols  = pyfits.ColDefs(save)
#
#--- combine the data
#
    #hdu = pyfits.BinTableHDU.from_columns(ocols + newcols)
    hdu = pyfits.BinTableHDU.from_columns(oldcols + newcols)
#
#--- create the new fits file
#
    mcf.rm_files(outfile)
    hdu.writeto(outfile)

#-----------------------------------------------------------------------------------------
#-- select_letg_arm: extract area defined by letg arm area either include or exclude     -
#-----------------------------------------------------------------------------------------

def select_letg_arm(ofits, outfile, include=1):
    """
    extract area defined by letg arm area either include or exclude
    input:  ofits   --- fits file name
            outfile --- output fits file name
            include --- if 1: include the area, 0: exclude the area
    output: outfile --- resulted fits file

    note: polygons are from:
        $CALDB/data/chandra/hrc/tgmask2/letgD1999-07-22regN0002.fits"[ROWID=SOURCE]"
    """
#
#--- letg polygon areas
#
    if include == 1:
        area1 = '(-0.0057803997770,  0.000531000027 '
        area1 = area1 + ',-0.29866999387741, 0.000531000027       '
        area1 = area1 + ',-1.1548999548,     0.00211400003172     '
        area1 = area1 + ',-1.1548999548,    -0.00211400003172     '
        area1 = area1 + ',-0.29866999387741,-0.000531000027       '
        area1 = area1 + ',-0.0057803997770, -0.000531000027      '
        area1 = area1 + ',-0.0057803997770,  0.000531000027)      '
    
        area2 = '(0.0057803997770,   0.000531000027'
        area2 = area2 + ',0.29866999387741,  0.000531000027       '
        area2 = area2 + ',    1.1548999548,  0.00211400003172     '
        area2 = area2 + ',    1.1548999548, -0.00211400003172     '
        area2 = area2 + ',0.29866999387741, -0.000531000027       '
        area2 = area2 + ',0.0057803997770,  -0.000531000027      '
        area2 = area2 + ',0.0057803997770,   0.000531000027)'
#
#--- letg polygon background area
#
    else:
        area1 = '(-0.0057803997770,   0.00730999978259  '
        area1 = area1 + ',-0.29866999387741,   0.00730999978259  '
        area1 = area1 + ',-1.1548999548,       0.02314000017941  '
        area1 = area1 + ',-1.1548999548,       0.00200000009499  '
        area1 = area1 + ',-0.29866999387741,   0.00200000009499  '
        area1 = area1 + ',-0.0057803997770,    0.00200000009499  '
        area1 = area1 + ',-0.0057803997770,    0.00730999978259) '
     
        area2 = '(0.0057803997770,   -0.00730999978259          '
        area2 = area2 + ',0.29866999387741,  -0.00730999978259  '
        area2 = area2 + ',1.1548999548,      -0.02314000017941  '
        area2 = area2 + ',1.1548999548,      -0.00200000009499  '
        area2 = area2 + ',0.29866999387741,  -0.00200000009499  '
        area2 = area2 + ',0.0057803997770,   -0.00200000009499  '
        area2 = area2 + ',0.0057803997770,   -0.00730999978259) '
#
#--- create dmcopy command
#
    cmd = 'dmcopy "' + ofits + '[(tg_r,tg_d)='
    cmd = cmd + ' polygon' + area1
    cmd = cmd + '+polygon' + area2
    cmd = cmd + ']" outfile=' +  outfile + ' clobber="yes"'

    os.system(cmd)

#-----------------------------------------------------------------------------------------
#-- get_bkg_subtracted_stats: create background adjusted statistics                     --
#-----------------------------------------------------------------------------------------

def get_bkg_subtracted_stats(data, bkg, tspan, area, barea, cent=0):
    """
    create background adjusted statistics
    input:  data    --- main data set
            bkg     --- backgroun data set
            tspan   --- time span in seconds
            area    --- main data area size in pixs
            barea   --- background data area size in pixs
    output: avg     --- average
            sig     --- sigma 
            med     --- median
            crate   --- data count rate
            brate   --- background count rate
            tspan   --- time span (in sec)
    """
#
#--- initialize the bins
#
    dbin = [0] * 300 
    bbin = [0] * 300 
    dtot = 0
    btot = 0
#
#--- center data histogram
#
    for k in range(0, len(data)):
        pos = int(data[k])
        if pos >= 300:
            continue

        dbin[pos] += 1
        dtot += 1
#
#--- background histogram
#
    for k in range(0, len(bkg)):
        pos = int(bkg[k])
        if pos >= 300:
            continue

        bbin[pos] += 1
        btot += 1
#
#--- compute weighted avg
#
    add  = 0
    wadd = 0 
    test = 0
    ddd  = 0
    for k in range(0, 300):
        weight = math.sqrt(dbin[k] + bbin[k]/(barea/area)**2)
        if weight == 0.0:
            continue

        add  += float(k) * weight
        wadd += weight
        test += 1./weight

    avg = add / wadd
#
#--- compute sig
#
    var = 0
    for k in range(0, len(data)):
        diff = data[k] - avg
        var += diff * diff

    sig = math.sqrt(var/len(data))
    ser = sig /math.sqrt(len(data))
#
#--- find med
#
    middle = 0.5 * dtot
    asum   = 0
    med    = avg
    for k in range(0, 300):
        asum += dbin[k]
        if asum >= middle:
            med = k
            break
#
#--- count rate per sec
#
    if cent == 1:
        crate = dtot / tspan 
        brate = btot / tspan 
        crate = crate - brate
        cerr  = math.sqrt(dtot) / tspan
        berr  = math.sqrt(btot) / tspan
    else:
        if area == 0:
            crate = dtot / tspan 
            cerr  = math.sqrt(dtot) / tspan
        else:
            crate = dtot / tspan /  area
            cerr  = math.sqrt(dtot) / tspan / area
        if barea == 0:
            brate = btot / tspan 
            berr  = math.sqrt(btot) /tspan
        else:
            brate = btot / tspan / barea
            berr  = math.sqrt(btot)  / tspan / barea

    return [avg, ser, med, crate, brate, tspan, cerr, berr]

#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    obsid = sys.argv[1]
    evt1  = sys.argv[2]
    out = extract_hz43_stat(obsid, evt1)
    if out:
        print("Data OK")
    else:
        print("Data FAILED")
