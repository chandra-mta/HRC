#!/usr/bin/env /data/mta/Script/Python3.6/envs/ska3/bin/python

#####################################################################################
#                                                                                   #
#   extract_arlac_stat.py: find out ArLac observations and update the data tables   #
#                                                                                   #
#           author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                   #
#           Last Update: Jan 12, 2021                                               #
#                                                                                   #
#####################################################################################

import sys
import os
import string
import re
import math
import time
import random
import astropy.io.fits  as pyfits
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts3/house_keeping/dir_list'

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
sys.path.append(hrc_dir)

import mta_common_functions as mcf
import filter_evt_file      as fef
import hrc_common_functions as hcf
#
#--- temp writing file name
#
rtail  = int(time.time() * random.random())
zspace = '/tmp/zspace' + str(rtail)

NULL    = 'NULL'
outlist = ['samp_p_list', 'samp_n_list', 'pi_p_list',  'pi_n_list']
pi      = math.pi

#-----------------------------------------------------------------------------------------
#-- run_process: find out arlac observations and update the data tables                  --
#-----------------------------------------------------------------------------------------

def run_process(out):
    """
    find out ArLac observations and update the data tables
    input:  out --- a list of lists of obsids of hrc i and hrc s
    output: <data_dir>/<pi/samp>_<p/n>_list_<i/s>_<#>, pi_center_list_<i/s>, samp_center_list_<i/s>
    """
    if (len(out[0]) > 0) or (len(out[1])  > 0):

        pdir = 'Past_data/'
        chk  = data_dir + pdir
        if not os.path.isdir(chk):
            cmd = 'mkdir -p ' + chk
            os.system(cmd)
#
#--- remove the older data
#
        cmd = 'rm -rf ' + data_dir + pdir + '/*'
        os.system(cmd)

        cmd = 'cp -f ' + data_dir + 'samp* ' + data_dir +  pdir + '.'
        os.system(cmd)
        cmd = 'cp -f ' + data_dir + 'pi* '   + data_dir +  pdir + '.'
        os.system(cmd)

    data   = out[0]         #--- a list of obsids of hrc i
    data2  = out[1]         #--- a list of obsids of hrc s
#
#--- analyze hrc i obsids
#
    for obsid in data:
        extract_arlac_stat(obsid, 'i')
#
#--- analyze hrc s obsids
#
    for obsid in data2:
        extract_arlac_stat(obsid, 's')

#-----------------------------------------------------------------------------------------
#-- extract_arlac_stat: extract arlac data, analyze and append the results to data table  --
#-----------------------------------------------------------------------------------------

def extract_arlac_stat(obsid, inst=''):
    """
    extract arlac data, analyze and append the results to data table
    input:  obsid   --- obsid
    output: see pinrt_out_<1/2>
    """
#
#--- create the saving directory if it does not exist
#
    sdir = data_dir + 'Fittings/' + str(obsid)
    if not os.path.isdir(sdir):
        cmd = 'mkdir -p  ' + sdir
        os.system(cmd)

    print(str(obsid))

    out = extract_data(obsid, inst)
    if out == False:
        print("Analysis Failed")
        cmd = 'rm *.fits*'
        os.system(cmd)
        ffile = house_keeping + 'skip_list'
        with open(ffile, 'a') as fo:
            fo.write(str(obsid) + '\n')

        return False
#
#---    out list below contains:
#---    hdata: [inst, date, tstart, tstop, exposure, ra_nom2, dec_nom2, 
#---            ra_targ2, dec_targ2, ra_est, dec_est,ra_nom, dec_nom, 
#---            roll_nom, ra_targ, dec_targ, ra_nom3, dec_nom3, ra_targ3, 
#---            dec_targ3, carea, barea]
#---    others: [x, y, pi, amp_sf, sumamps, samp] all of them are lists
#

    [hdata, cent_data, cent_bkg] = out
#
    tspan = hdata[4]                               #--- exposure time
#
#--- write down the center part results: [avg, sig, med, crate, brate, tspan]
#
    if (len(cent_data) > 0) and (len(cent_data[1]) >=50):

        pi_stats = get_bkg_subtracted_stats(cent_data[2], cent_bkg[2], tspan, hdata[-2], hdata[-1])
        sp_stats = get_bkg_subtracted_stats(cent_data[5], cent_bkg[5], tspan, hdata[-2], hdata[-1])

        save_ind_dist_center(sdir, cent_data[2], cent_data[5])

        oname = 'pi_list'
        print_out(obsid, hdata, pi_stats, oname)
    
        oname = 'samp_list'
        print_out(obsid, hdata, sp_stats, oname)
    else:
        print("Obs did not meet the criteria.")
        sfile = house_keeping + 'skip_list'
        with open(sfile, 'a') as fo:
            fo.write(str(obsid) + '\n')
#
#--- clean up 
#
    cmd = 'rm *.fits*'
    os.system(cmd)

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
    [inst, date, tstart, tstop, exposure, ra_nom_d, dec_nom_d, ra_targ_d,\
     dec_targ_d, ra_est, dec_est, ra_nom, dec_nom, roll_nom, ra_targ, dec_targ,\
     detx_nom, dety_nom, detx_targ, dety_targ, carea, barea] = head

    [avg, sig, med, crate, brate, tspan] = stats

    [yoffset, zoffset] = compute_y_z_offset(ra_nom, dec_nom, roll_nom, ra_targ, dec_targ)

    line = str(tstart)                                 + '\t'       #0
    line = line + format(int(float(obsid)), '5d')      + '\t'       #1
    line = line + str(date)                            + '\t'       #2
    line = line + format(int(exposure), '5d')          + '\t'       #3
    line = line + format(round(ra_nom_d, 2),   '.2f')  + '\t'       #4
    line = line + format(round(dec_nom_d, 2),  '.2f')  + '\t'       #5
    line = line + format(round(ra_targ_d, 2),  '.2f')  + '\t'       #6
    line = line + format(round(dec_targ_d, 2), '.2f')  + '\t'       #7

    line = line + format(round(ra_est,  2),    '.2f')  + '\t'       #8
    line = line + format(round(dec_est, 2),    '.2f')  + '\t'       #9

    line = line + format(round(yoffset, 8),    '.8f')  + '\t'       #10
    line = line + format(round(zoffset, 8),    '.8f')  + '\t'       #11

    line = line + format(round(avg, 3),        '.3f')  + '\t'       #12
    line = line + format(round(sig, 3),        '.3f')  + '\t'       #13
    line = line + str(int(med))                        + '\t'       #14
    line = line + format(round(crate, 8),      '.8f')  + '\t'       #15
    line = line + format(round(brate, 8),      '.8f')  + '\t'       #16

    line = line + format(round(ra_nom,   10), '.10f')  + '\t'       #17
    line = line + format(round(dec_nom,  10), '.10f')  + '\t'       #18
    line = line + format(round(roll_nom, 10), '.10f')  + '\t'       #19
    line = line + format(round(ra_targ,  10), '.10f')  + '\t'       #20
    line = line + format(round(dec_targ, 10), '.10f')  + '\t'       #21
    line = line + format(round(detx_nom, 2),   '.2f')  + '\t'       #22
    line = line + format(round(dety_nom, 2),   '.2f')  + '\t'       #23
    line = line + format(round(detx_targ, 2),  '.2f')  + '\t'       #24
    line = line + format(round(dety_targ, 2),  '.2f')  + '\t'       #25
    line = line + format(round(carea, 2),      '.2f')  + '\t'       #26
    line = line + format(round(barea, 2),      '.2f')  + '\n'       #27

    if inst.lower() == 'hrc-i':
        outfile = data_dir + oname + '_i'
    else:
        outfile = data_dir + oname + '_s'

    with open(outfile, 'a') as fo:
        fo.write(line)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def convert_to_coordinates(ra, dec, evt, ctype):
    """
    convert ra dec to sky, det, or chip coordinates
    input:  ra      --- ra in degree
            dec     --- dec in degree
            evt     --- evt1 file 
            ctype   --- sky/det/chip
    output: [newx, newy]
    """
    cmd = 'dmcoords ' + evt + ' opt=cel ra=' + str(ra) + ' dec=' 
    cmd = cmd + str(dec) + ' verbose=1 > ' + zspace
    hcf.run_ascds(cmd)

    newx= ''
    newy= ''
#
#--- extract sky coordindats
#
    info = mcf.read_data_file(zspace, remove=1)
    info.reverse()
    for ent in info:
        if ctype == 'sky':
            mc = re.search('SKY\(X,Y\)', ent)
        elif ctype == 'det':
            mc = re.search('DETX,DETY',  ent)
        elif ctype == 'chip':
            mc = re.search('CHIP',       ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            newx = float(atemp[1])
            newy = float(atemp[2])

    return [newx, newy]

#-----------------------------------------------------------------------------------------
#-- extract_data: extract needed information for a given observation                    --
#-----------------------------------------------------------------------------------------

def extract_data(obsid, inst):
    """
    extract needed information for a given observation 
    input:  obsid   --- obsid
    output: a list of lists:
                    [inst, date, start,  crsv, crsu, amp_sf, pha, pi, sumamps, samp, chip_id, 
                     tg_m, tg_mlam, tg_srcid, tg_part, tg_smap]
            note: inst,date, start are strings, but all others are lists.
    """
#
#--- extract evt1 data
#
    evt1 = hcf.run_arc5gl(0, 0, obsid = obsid, operation='retrieve', level ='1', filetype='evt1')
    if evt1 == "":
        return False
#
#--- find which instrument
#
    hdr      = pyfits.getheader(evt1, 1)
    inst     = hdr['detnam']
    date     = hdr['date-obs']
    start    = hdr['tstart']
    stop     = hdr['tstop']
    try:
        expo = hdr['exposure']
    except:
        return False

    ra_nom   = float(hdr['ra_nom'])
    dec_nom  = float(hdr['dec_nom'])
    roll_nom = float(hdr['roll_nom'])
    ra_targ  = float(hdr['ra_targ'])
    dec_targ = float(hdr['dec_targ'])

    [ra_nom2,  dec_nom2]  = convert_to_coordinates(ra_nom,  dec_nom,  evt1, 'sky')
    [ra_targ2, dec_targ2] = convert_to_coordinates(ra_targ, dec_targ, evt1, 'sky')

    [ra_nom3,  dec_nom3]  = convert_to_coordinates(ra_nom,  dec_nom,  evt1, 'det')
    [ra_targ3, dec_targ3] = convert_to_coordinates(ra_targ, dec_targ, evt1, 'det')

    hdata  = [inst, date, start, stop, expo, ra_nom2, dec_nom2, ra_targ2, dec_targ2]
    hdata2 = [ra_nom, dec_nom, roll_nom, ra_targ, dec_targ]
    hdata3 = [ra_nom3, dec_nom3, ra_targ3, dec_targ3]

    [yoffset, zoffset] = compute_y_z_offset(ra_nom, dec_nom, roll_nom, ra_targ, dec_targ)
#
#--- clean the data
#
    clean1 = 'evt1_clean.fits'
    try:
        fef.filter_evt_file(evt1, clean1, obsid)
    except:
        return False
#
#--- get the target area
#
    carea = 'clean_targ.fits'
    barea = 'clean_targ_bkg.fits'
    chk   = get_target_data(obsid, evt1, clean1, carea, barea, ra_targ,\
                            dec_targ, yoffset, zoffset, inst)

    if chk:
        hdata.append(chk[0])
        hdata.append(chk[1])
        tarea = chk[2]              #--- center area in pix
        sarea = chk[3]              #--- background area in pix
        cent_data = get_data(carea, ['x', 'y', 'pi','amp_sf', 'sumamps'])
        samp      = compute_samp(cent_data[4], cent_data[3], inst)
        cent_data.append(samp)

        cent_bkg  = get_data(barea, ['x', 'y', 'pi','amp_sf', 'sumamps'])
        samp      = compute_samp(cent_bkg[4], cent_bkg[3], inst)
        cent_bkg.append(samp)
    else:
        hdata.append(0.0)
        hdata.append(0.0)
        cent_data = []
        cent_bkg  = []
        tarea = 0
        sarea = 0

    hdata = hdata + hdata2 + hdata3
    hdata.append(tarea)
    hdata.append(sarea)

    return [hdata, cent_data, cent_bkg]

#-----------------------------------------------------------------------------------------
#-- get_target_data: extract center part of the data and the background                 --
#-----------------------------------------------------------------------------------------

def get_target_data(obsid, evt1, clean, carea, barea, ra_targ, dec_targ, yoffset, zoffset, inst):
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
#--- set the radius around the source and the background area
#
    [radius, annul1, annul2, tarea, sarea, cx, cy] = set_radius(yoffset, zoffset, inst)
#
#--- check whether manually determined sky coordinates are available 
#
    try:
        #cmd = ' celldetect ' + evt1 + ' outfile=zinfo.fits clobber=yes'
        cmd = ' tgdetect ' + evt1 + ' none outfile=zinfo.fits clobber=yes'
        hcf.run_ascds(cmd)
    
        t     = pyfits.open('zinfo.fits')
        tdata = t[1].data
    except:
        return  False

    try:
        data  = tdata['net_counts'][0]
    except:
        return False

    psave = []
    for ent in ('x','y', 'net_counts', 'net_counts_err', 'bkg_counts', 'bkg_counts_err',\
                'net_rate', 'net_rate_err', 'bkg_rate', 'bkg_rate_err'):
        try:
            data  = tdata[ent][0]
            psave.append(data)
        except:
            pass
    t.close()

    x = psave[0]
    y = psave[1]

    mcf.rm_files('zinfo.fits')
#
#--- center part
#
    cmd   = 'dmcopy "' + clean + '[(x,y)=circle(' + str(x) + ',' + str(y) 
    cmd   = cmd + ',' + str(radius) + ')]" outfile=' + carea + ' clobber=yes'
#
#--- sometime ascds does not work with dmcopy. if so, use ciao
#
    hcf.run_ciao(cmd)
    if not os.path.isfile(carea):
        hcf.run_ascds(cmd)

#
#--- background area; annulus around the source position
#
    if cx == 0:     #--- hrc s case and hrc i offset < 13.6 arcmin
        px = x
        py = y
    else:           #--- hrc i with offset > 13.6 arcmin
        px = cx
        py = cy

    cmd   = 'dmcopy "' + clean + '[(x,y)=annulus(' + str(px) + ',' + str(py)  
    cmd   = cmd + ',' +  str(annul1) + ',' + str(annul2)
    cmd   = cmd + ')]" outfile=' + barea + ' clobber=yes'

    hcf.run_ciao(cmd)
    if not os.path.isfile(barea):
        hcf.run_ascds(cmd)
    
    line  = str(psave[0])
    for  k in range(1, len(psave)):
        line  = line + '\t'  + str(psave[k])

    line = line + '\n'

    ofile = data_dir + 'Fittings/' + str(obsid) + '/center_info'
    with open(ofile, 'w') as fo:
        fo.write(line)

    return [float(x), float(y), tarea, sarea]

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
    output: <sidr>/samp_p_list<#> etc
    """
    try:
        oname = sdir + '/samp_list'
        print_out_dist(samp_list, oname)
    except:
        pass

    try:
        oname = sdir + '/pi_list'
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

    with open(oname, 'w') as fo:
        for ent in odata:
            fo.write(str(ent) + '\n')

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
        save.append(tdata[col])

    hd.close()

    return save

#-----------------------------------------------------------------------------------------
#-- get_bkg_subtracted_stats: create background adjusted statistics                     --
#-----------------------------------------------------------------------------------------

def get_bkg_subtracted_stats(data, bkg, tspan, carea, barea):
    """
    create background adjusted statistics
    input:  data    --- main data set
            bkg     --- backgroun data set
            tspan   --- time span in seconds
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
        weight = math.sqrt(dbin[k] + bbin[k] * (carea/barea)**2)
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
    crate = dtot / tspan 
    brate = btot / tspan  

    return [avg, ser, med, crate, brate, tspan]

#-----------------------------------------------------------------------------------------
#-- add_coldata_to_fits: adding data from the second fits file to the first fits file   --
#-----------------------------------------------------------------------------------------

def add_coldata_to_fits(ofits1, ofits2, col_names, outfile):
    """
    adding data from the second fits file to the first fits file
    input:  ofits1  --- the first fits file
    ofits2  --- the second fits file
    col_names   --- a list of column names which you want to copy from the 
    second fits file
    outfile --- output fits file name
    output: outfie  --- resulted fits file
    """
    otable = pyfits.open(ofits1)[1].data
    ocols  = otable.columns
    atable = pyfits.open(ofits2)[1].data
    acols  = atable.columns
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
    hdu = pyfits.BinTableHDU.from_columns(ocols + newcols)
#
#--- create the new fits file
#
    mcf.rm_files(outfile)
    hdu.writeto(outfile)

#-----------------------------------------------------------------------------------------
#-- extract_offset_info: extract yoffset and zoffset information                      ----
#-----------------------------------------------------------------------------------------

def extract_offset_info(obsid, inst): 
    """
    extract yoffset and zoffset information
    input:  obsid   --- obsid
            inst    --- instrument either 'hrc-i' or 'hrc-s'
    output: [yoffset, zoffset]
    """

    if inst.lower() == 'hrc-i':
        ifile = '/data/hrc/i/datahrcidx.fits'
    else:
        ifile = '/data/hrc/s/datahrcsidx.fits'

    t     = pyfits.open(ifile)
    tdata = t[1].data
    t.close()

    obsids = tdata['obsid'][0]
    y_set  = tdata['yoffset'][0]
    z_set  = tdata['zoffset'][0]

    yoffset = 0
    zoffset = 0
    obsid = int(float(obsid))

    for k in range(0, len(obsids)):
        val = int(float(obsids[k]))
        if val == obsid:
            yoffset = float(y_set[k])
            zoffset = float(z_set[k]) 
            break

    return [yoffset, zoffset]

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def compute_y_z_offset(ra_nom, dec_nom, roll_nom, ra_targ, dec_targ):

    const1  = pi / 180.0
    const2  = 3600.0/0.13175
    const3  = 0.13175/60.0

    roll    = -(roll_nom - 180.0)   * const1
    dx      = -(ra_targ - ra_nom)   * math.cos(dec_nom * const1) * const2
    dy      =  (dec_targ - dec_nom) * const2
    xroll   =  dx * math.cos(roll) + dy  * math.sin(roll) 
    yroll   = -dx * math.sin(roll) + dy  * math.cos(roll)
    yoffset = xroll * const3
    zoffset = yroll * const3

    return[yoffset, zoffset]

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def set_radius(yoffset, zoffset, inst):

    const3  = 0.13175/60.0

    radius  = (2.0* math.sqrt(yoffset**2+zoffset**2)) / 60. 
     
    if inst == 'i' and radius > 13.5:
        x       = 16363
        y       = 16591
        annul1  = 4780
        annul2  = 6145
        radius /= const3
    else:
        x      = 0
        y      = 0

        if radius <  0.417:
            radius = 0.417
#
#--- convert into pixel
#
        radius /= const3

        annul1 = 2.0 * radius
        annul2 = 4.0 * radius
#
#--- compute the area in pixels
#
    carea = pi * radius**2
    barea = pi * (annul2**2 - annul1**2)

    return [radius, annul1, annul2, carea, barea, x, y]


#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) == 2:
        obsid = sys.argv[1]
        obsid.strip()

        if obsid.lower() == 'all':
            run_process('all')
        else:
            run_process('new')

    elif len(sys.argv) == 3:
        obsid = sys.argv[1]
        obsid = int(float(obsid))
        inst  = sys.argv[2]
        extract_arlac_stat(obsid, inst)
        exit(1)

    else:
        run_process('new')

