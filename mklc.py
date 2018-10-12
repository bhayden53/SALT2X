import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo
from Salt2X import *
import helpers
import astropy
import sncosmo
import sys
from IPython import embed
import progressbar as pb
from astropy.io import fits
import commands

ar = 0.115
af = 0.03
beta = 3.1
MB = -19.1
sigint = 0.1

mean = [0,0,0]
cov = [[ 1   ,0.74,0],
        [0.74,1   ,2e-2],
        [0   ,2e-2,5e-3]]

def save_img(dat, imname):

    commands.getoutput("rm -f " + imname)
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()

lcfile = './JLA/jla_light_curves/jla_lc.txt'
lc = np.genfromtxt(lcfile, dtype=None)
lcfits = Table.read('./fit_backups/emcee_jla_2s_skewfit_100716.dat', format='ascii')
modeldir = os.environ['SNCOSMO_MODELDIR']
helpers.get_JLA_bandpasses()
helpers.register_JLA_magsys()
standard_zps = {'STANDARD::U':9.724, 'STANDARD::B':9.907, 'STANDARD::V':9.464, 'STANDARD::R':9.166, 'STANDARD::I':8.846}
FourShooter_zps = {'4SHOOTER2::Us':9.724, '4SHOOTER2::B':9.8744, '4SHOOTER2::V':9.4789, '4SHOOTER2::R':9.1554, '4SHOOTER2::I':8.8506}
keplercam_zps = {'KEPLERCAM::Us':9.6922, 'KEPLERCAM::B':9.8803, 'KEPLERCAM::V':9.4722, 'KEPLERCAM::r':9.3524, 'KEPLERCAM::i':9.2542}
swope_zps = {'SWOPE2::u':10.514, 'SWOPE2::g':9.64406, 'SWOPE2::r':9.3516, 'SWOPE2::i':9.2500, 'SWOPE2::B':9.876433, 'swope2::v_lc3009':9.471276, 'swope2::v_lc3014':9.476626, 'swope2::v_lc9844':9.477482}
sdss_zps = {'SDSS::u':0.06791, 'SDSS::g':-0.02028, 'SDSS::r':-0.00493, 'SDSS::i':-0.01780, 'SDSS::z':-0.01015}

def do_stuff(ctr):
    p = pb.ProgressBar(maxval=740, widgets = [pb.Percentage(),pb.Bar(),pb.ETA()]).start()
    pbctr = 0
    for i,l in enumerate(lc):
        p.update(pbctr)
        pbctr += 1

        restcut = (3000,7000)
        data = sncosmo.read_lc(l, format='salt2')
        try:
            z = data.meta['Redshift']
        except:
            pass
        try:
            z = data.meta['Z_CMB']
        except:
            pass
        try:
            survey = data.meta['SURVEY']
        except:
            pass
        try:
            survey = data.meta['SET']
        except:
            pass
        nickname = data.meta['SN']
        try:
            nickname = str(int(nickname))
        except:
            pass
        mwebv = data.meta['MWEBV']
        dust = sncosmo.CCM89Dust()
        data = astropy.table.Table(data, masked=True)
        #rename columns so that my fitter can handle things
        data.rename_column('Filter', 'tmp')
        data.rename_column('Date', 'time')
        data.rename_column('Flux', 'flux')
        data.rename_column('Fluxerr', 'fluxerr')
        data.rename_column('MagSys', 'zpsys')
        data.rename_column('ZP', 'zp')

        if survey == 'SNLS':
            sn_nickname = l.split('/')[-1].split('.')[0].split('-')[-1]
            band = []
            for j, bp in enumerate(data['tmp']):
                band.append( '%s-%s' %(sn_nickname, bp) )
            band = astropy.table.Column(band, name='band')
            data.add_column(band)
            data.remove_column('tmp')
        else:
            data.rename_column('tmp', 'band')

        # deal with swope filters
        mask = (data['band'] == 'SWOPE2::V')
        nswopev = len(mask.nonzero()[0])
        if nswopev > 0:
            band = []
            for j, bp in enumerate(data['band']):
                if (bp == 'SWOPE2::V'):
                    if (data['time'][j] < 53749):
                        band.append('swope2::v_lc3014')
                    elif (data['time'][j] < 53760):
                        band.append('swope2::v_lc3009')
                    else:
                        band.append('swope2::v_lc9844')
                else:
                    band.append(bp)
            data.remove_column('band')
            band = astropy.table.Column(band, name='band')
            data.add_column(band)

            ind = np.where( (data['band'] == 'SWOPE2::V') & (data['time']>53749.) & ((data['time']<=53760.)) )
            data['band'][ind] = 'swope2::v_lc3009'
            ind = np.where( (data['band'] == 'SWOPE2::V') & (data['time']>53760.) )
            data['band'][ind] = 'swope2::v_lc9844'
            # print ind

        #deal with filter coverage
        #also deal with STANDARD filter zeropoints
        unique_bands = np.unique(data['band'])
        fit_bands = []
        nofit_bands = []
        # print unique_bands
        tmperr = np.copy(data['fluxerr'])
        for ub in unique_bands:
            # print ub
            bp = sncosmo.get_bandpass(ub)
            rest = bp.wave_eff / (1.0+z)
            if (rest >= restcut[0]) & (rest <= restcut[1]):
                fit_bands.append(ub)
            else:
                nofit_bands.append(ub)
            if 'STANDARD' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(standard_zps[ub])
                errcor = 10**(-0.4*standard_zps[ub])
                data['fluxerr'][ind] *= errcor
            if '4SHOOTER2' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(FourShooter_zps[ub])
                errcor = 10**(-0.4*FourShooter_zps[ub])
                data['fluxerr'][ind] *= errcor
            if 'KEPLERCAM' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(keplercam_zps[ub])
                errcor = 10**(-0.4*keplercam_zps[ub])
                data['fluxerr'][ind] *= errcor
                # print ub
                # print data['zp'][ind]
            if 'swope' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(swope_zps[ub])
                errcor = 10**(-0.4*swope_zps[ub])
                data['fluxerr'][ind] *= errcor
            if 'sdss' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(sdss_zps[ub])
                errcor = 10**(-0.4*sdss_zps[ub])
                data['fluxerr'][ind] *= errcor

        for nfb in nofit_bands:
            mask = data['band'] == nfb
            for c in data.colnames:
                data[c].mask = (data[c].mask | mask)

        mwebv = data.meta['MWEBV']

        mask = data['band'].mask.nonzero()[0]
        data.remove_rows(mask)

        ind = np.where(lcfits['SN'] == nickname)
        t0 = lcfits['DayMax'][ind][0]

        x1r, x1f, c = np.random.multivariate_normal(mean,cov,size=1)[0]
        mu = cosmo.distmod(z).value
        absmag = MB - ar*x1r - af*x1f + beta*c + np.random.normal(scale=sigint, size=1)[0]
        mB = mu + absmag

        source = Salt2XSource(version='2.4', modeldir=modeldir)
        model  = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])

        model.set(z=z, x1=x1f, s=x1r, c=c, t0=t0, mwebv=mwebv)
        model.set_source_peakabsmag(absmag, 'bessellb', 'ab')
        flux = model.bandflux(data['band'], data['time'], zp=data['zp'], zpsys=data['zpsys'])
        whocares, saltcov = model.bandfluxcov(data['band'], data['time'], data['zp'],data['zpsys'])

        # handle model cov blowups
        saltcov = np.copy(saltcov)
        diag = np.copy(saltcov.diagonal())
        model_snr = whocares/np.sqrt(diag)

        ind = np.where((np.abs(model_snr) < 1) & ~np.isnan(model_snr))

        diag[ind] = diag[ind] * np.abs(model_snr[ind])**2
        np.fill_diagonal(saltcov, diag)


        diagerr = np.diag(data['fluxerr']**2)
        fullcov = saltcov + diagerr
        try:
            np.linalg.cholesky(fullcov)
        except:
            print 'Cholesky failed... exiting'
            sys.exit()
            
        noise = np.random.multivariate_normal(np.zeros(len(diagerr)), fullcov, size=1)[0]

        #lower zp, lower flux
        data['flux'] = flux + noise

        data.meta['x1r'] = x1r
        data.meta['x1f'] = x1f
        data.meta['c'] = c
        data.meta['alpha_r'] = ar
        data.meta['alpha_f'] = af
        data.meta['beta'] = beta
        data.meta['MB'] = MB
        data.meta['mB'] = mu+MB
        data.meta['DayMax'] = t0
        data.meta['cosmology'] = 'Planck15'

        data.rename_column('band', 'Filter')
        data.rename_column('time', 'Date')
        data.rename_column('flux', 'Flux')
        data.rename_column('fluxerr', 'Fluxerr')
        data.rename_column('zpsys', 'MagSys')
        data.rename_column('zp', 'ZP')

        if survey == 'SNLS':
            for row in data:
                tmp = row['Filter']
                ind = tmp.find('MEGACAM')
                row['Filter'] = row['Filter'][ind:]

        sncosmo.write_lc(data,'./cadence_sim/lc/%s_%s.list' %(nickname, ctr), format='salt2')
    p.finish()

do_stuff('A')
do_stuff('B')
do_stuff('C')
do_stuff('D')

