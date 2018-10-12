import numpy as np
from IPython import embed
import sys
import sncosmo
import astropy
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--jla', dest='jla', action='store_true')
parser.add_argument('--cadencesim', dest='cadencesim', action='store_true')
args = parser.parse_args()

def clean_lines(lines, stringlist = [""]):
    # Start by getting rid of spaces.
    lines = [item.strip() for item in lines]

    # Check for strings to exclude.
    lines = [item for item in lines if stringlist.count(item) == 0]

    # Get rid of comments
    lines = [item for item in lines if item[0] != "#"]
    return lines

def read_param(flnm, param, default = None, ind = 1):
    # This one is commonly used for photometry parameter files.
    fp = open(flnm)
    lines = fp.read()
    fp.close()

    lines = lines.split('\n')
    lines = clean_lines(lines)

    for line in lines:
        parsed = line.split(None)
        if parsed[0] == param:
            # print "Reading " + param + " from " + flnm

            try:
                # Yeah, I know eval is bad. But it works with all types!
                return eval(parsed[ind])
            except:
                return parsed[ind]
    print
    print "Couldn't find ", param, flnm
    print "Returning default ", default
    print

    return default

    #################################################################################
if args.jla:
    lcfile = np.genfromtxt('./JLA/jla_light_curves/jla_lc.txt', dtype=None)
if args.cadencesim:
    lcfile = np.genfromtxt('./cadence_sim/lc/sim_lc.txt', dtype=None)
    #################################################################################


names = ['SN', 'survey', 'z', 'chisq', 'dof', 'DayMax', 'RestFrameMag_0_B', 'X0', 'X1', 'Color', 'Skew', 'DayMax_coverr', 'DayMax_err', 'RestFrameMag_0_B_coverr',
    'RestFrameMag_0_B_err', 'X0_coverr','X0_err', 'X1_coverr', 'X1_err', 'Skew_coverr','Skew_err', 'Color_coverr',
    'Color_err','CovDayMaxX0', 'CovDayMaxRestFrameMag_0_B', 'CovDayMaxX1', 'CovDayMaxSkew', 'CovDayMaxColor', 'CovX0X1', 'CovX0Skew',
    'CovX0Color','CovRestFrameMag_0_BX1', 'CovRestFrameMag_0_BSkew', 'CovRestFrameMag_0_BColor', 'CovX1Skew', 'CovX1Color', 'CovSkewColor',
    'DayMax_p2', 'DayMax_p15', 'DayMax_p50', 'DayMax_p84', 'DayMax_p97',
    'x0_p2', 'x0_p15', 'x0_p50', 'x0_p84', 'x0_p97',
    'mB_p2', 'mB_p15', 'mB_p50', 'mB_p84', 'mB_p97',
    'x1_p2', 'x1_p15', 'x1_p50', 'x1_p84', 'x1_p97',
    's_p2', 's_p15', 's_p50', 's_p84', 's_p97',
    'Color_p2', 'Color_p15', 'Color_p50', 'Color_p84', 'Color_p97']
dtype = ('S20', 'S10', np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64,
    np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64,
    np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64,
    np.float64, np.float64, np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64,
    np.float64,np.float64,np.float64,np.float64,np.float64)
outdata = astropy.table.Table(masked=True, names=names, dtype=dtype)
for i, lc in enumerate(lcfile):
    #################################################################################
    if args.jla:
        data = sncosmo.read_lc(lc, format='salt2')
        nickname = data.meta['SN']
    if args.cadencesim:
        data = sncosmo.read_lc(lc, format='ascii')
        nickname = lc.split('.')[0].split('/')[-1]
    #################################################################################

    try:
        nickname = str(int(nickname))
    except:
        pass
    #################################################################################
    if args.jla:
        fitfile = './fit_results/emcee/JLA/%s/%s.dat' %(nickname, nickname)
    if args.cadencesim:
        fitfile = './fit_results/emcee/cadencesim/%s/%s.dat' %(nickname, nickname)
    #################################################################################
    if not os.path.isfile(fitfile):
        print 'output file does not exist! %s' %fitfile
        continue
    z = read_param(fitfile, 'Redshift')
    tmax = read_param(fitfile, 'DayMax')
    mB = read_param(fitfile, 'RestFrameMag_0_B')
    x1 = read_param(fitfile, 'X1')
    x0 = read_param(fitfile, 'X0')
    c = read_param(fitfile, 'Color')
    s = read_param(fitfile, 'Skew')
    tm_err = np.sqrt(read_param(fitfile, 'CovDayMaxDayMax'))
    tm_err2 = np.sqrt(read_param(fitfile, 'CovDayMaxDayMax2'))
    nan_test = read_param(fitfile, 'CovRestFrameMag_0_BRestFrameMag_0_B')
    nan_test2 = read_param(fitfile, 'CovRestFrameMag_0_BRestFrameMag_0_B2')
    if (nan_test == 'nan') | (nan_test2 == 'nan'):
        mB_err = nan_test
        mB_err2 = nan_test
    else:
        mB_err = np.sqrt(read_param(fitfile, 'CovRestFrameMag_0_BRestFrameMag_0_B'))
        mB_err2 = np.sqrt(read_param(fitfile, 'CovRestFrameMag_0_BRestFrameMag_0_B2'))
    chisq = read_param(fitfile, 'CHISQ')
    dof = read_param(fitfile, 'DoF')
    survey = read_param(fitfile, 'survey')
    x1_err = np.sqrt(np.float(read_param(fitfile, 'CovX1X1')))
    x1_err2 = np.sqrt(read_param(fitfile, 'CovX1X12'))
    x0_err = np.sqrt(read_param(fitfile, 'CovX0X0'))
    x0_err2 = np.sqrt(read_param(fitfile, 'CovX0X02'))
    s_err = np.sqrt(np.float(read_param(fitfile, 'CovSkewSkew')))
    s_err2 = np.sqrt(read_param(fitfile, 'CovSkewSkew2'))
    c_err = np.sqrt(read_param(fitfile, 'CovColorColor'))
    c_err2 = np.sqrt(read_param(fitfile, 'CovColorColor2'))
    t0_x0 = read_param(fitfile, 'CovDayMaxX0')
    t0_mB = read_param(fitfile, 'CovDayMaxRestFrameMag_0_B')
    t0_x1 = read_param(fitfile, 'CovDayMaxX1')
    t0_s  = read_param(fitfile, 'CovDayMaxSkew')
    t0_c  = read_param(fitfile, 'CovDayMaxColor')
    x0_x1 = read_param(fitfile, 'CovX0X1')
    x0_s  = read_param(fitfile, 'CovX0Skew')
    x0_c  = read_param(fitfile, 'CovX0Color')
    mB_x1 = read_param(fitfile, 'CovRestFrameMag_0_BX1')
    mB_s  = read_param(fitfile, 'CovRestFrameMag_0_BSkew')
    mB_c  = read_param(fitfile, 'CovRestFrameMag_0_BColor')
    x1_s  = read_param(fitfile, 'CovX1Skew')
    x1_c  = read_param(fitfile, 'CovX1Color')
    s_c  = read_param(fitfile, 'CovSkewColor')
    t0_p2  = read_param(fitfile, 'x0_p2')
    t0_p15 = read_param(fitfile, 'x0_p15')
    t0_p50 = read_param(fitfile, 'x0_p50')
    t0_p84 = read_param(fitfile, 'x0_p84')
    t0_p97 = read_param(fitfile, 'x0_p97')
    x0_p2  = read_param(fitfile, 'x0_p2')
    x0_p15 = read_param(fitfile, 'x0_p15')
    x0_p50 = read_param(fitfile, 'x0_p50')
    x0_p84 = read_param(fitfile, 'x0_p84')
    x0_p97 = read_param(fitfile, 'x0_p97')
    mB_p2  = read_param(fitfile, 'mB_p2')
    mB_p15 = read_param(fitfile, 'mB_p15')
    mB_p50 = read_param(fitfile, 'mB_p50')
    mB_p84 = read_param(fitfile, 'mB_p84')
    mB_p97 = read_param(fitfile, 'mB_p97')
    x1_p2  = read_param(fitfile, 'x1_p2')
    x1_p15 = read_param(fitfile, 'x1_p15')
    x1_p50 = read_param(fitfile, 'x1_p50')
    x1_p84 = read_param(fitfile, 'x1_p84')
    x1_p97 = read_param(fitfile, 'x1_p97')
    s_p2  = read_param(fitfile, 's_p2')
    s_p15 = read_param(fitfile, 's_p15')
    s_p50 = read_param(fitfile, 's_p50')
    s_p84 = read_param(fitfile, 's_p84')
    s_p97 = read_param(fitfile, 's_p97')
    c_p2  = read_param(fitfile, 'Color_p2')
    c_p15 = read_param(fitfile, 'Color_p15')
    c_p50 = read_param(fitfile, 'Color_p50')
    c_p84 = read_param(fitfile, 'Color_p84')
    c_p97 = read_param(fitfile, 'Color_p97')


    outdata.add_row([str(nickname),survey,z, chisq, dof, tmax, mB, x0, x1,c, s, tm_err, tm_err2, mB_err, mB_err2, x0_err, x0_err2, x1_err, x1_err2, s_err, s_err2, c_err, c_err2, t0_x0, t0_mB, t0_x1, t0_s, t0_c, x0_x1, x0_s, x0_c, mB_x1, mB_s, mB_c, x1_s, x1_c, s_c,
    t0_p2, t0_p15, t0_p50, t0_p84, t0_p97,
    x0_p2, x0_p15, x0_p50, x0_p84, x0_p97,
    mB_p2, mB_p15, mB_p50, mB_p84, mB_p97,
    x1_p2, x1_p15, x1_p50, x1_p84, x1_p97,
    s_p2, s_p15, s_p50, s_p84, s_p97,
    c_p2, c_p15, c_p50, c_p84, c_p97])

outdata.sort('survey')

if args.jla:
    outdata.write('./fit_backups/jla_fits.dat', format='ascii.fixed_width', delimiter=' ')
if args.cadencesim:
    outdata.write('./fit_backups/cadencesim_fits.dat', format='ascii.fixed_width', delimiter=' ')
print outdata