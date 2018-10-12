import numpy as np
import sncosmo
import matplotlib as mpl
from astropy.table import Table, Column
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
from Salt2X import *
import helpers
import fitters
import emcee
import sys
from iminuit import describe
import argparse
import os
import traceback
import astropy
from IPython import embed
from sncosmo import photdata
import time

# change
import sfdmap

parser = argparse.ArgumentParser()
parser.add_argument('--emcee',  dest='emcee',  action='store_true')
parser.add_argument('--jla',    dest='jla',    action='store_true')
parser.add_argument('--cadencesim', dest='cadencesim',    action='store_true')
parser.add_argument('--noskew', dest='noskew', action='store_true')
parser.add_argument('--nsamp',  dest='nsamp',  default=5000, type=float)
parser.add_argument('--specific', '-s', dest='specific', default=None, type=str)
args = parser.parse_args()

modeldir = os.environ['SNCOSMO_MODELDIR']
scratch = os.environ['SCRATCH']
scratch = '.'

def emcee_chain_maxlike( chain, key ):
    maxlike = np.argmax( chain['lnprob'] )
    return chain[key][maxlike]

if args.jla:
    # zeropoints taken from the JLA magnitude system files
    standard_zps = {'STANDARD::U':9.724, 'STANDARD::B':9.907, 'STANDARD::V':9.464, 'STANDARD::R':9.166, 'STANDARD::I':8.846}
    FourShooter_zps = {'4SHOOTER2::Us':9.724, '4SHOOTER2::B':9.8744, '4SHOOTER2::V':9.4789, '4SHOOTER2::R':9.1554, '4SHOOTER2::I':8.8506}
    keplercam_zps = {'KEPLERCAM::Us':9.6922, 'KEPLERCAM::B':9.8803, 'KEPLERCAM::V':9.4722, 'KEPLERCAM::r':9.3524, 'KEPLERCAM::i':9.2542}
    swope_zps = {'SWOPE2::u':10.514, 'SWOPE2::g':9.64406, 'SWOPE2::r':9.3516, 'SWOPE2::i':9.2500, 'SWOPE2::B':9.876433, 'swope2::v_lc3009':9.471276, 'swope2::v_lc3014':9.476626, 'swope2::v_lc9844':9.477482}
    sdss_zps = {'SDSS::u':0.06791, 'SDSS::g':-0.02028, 'SDSS::r':-0.00493, 'SDSS::i':-0.01780, 'SDSS::z':-0.01015}
    
    # registering the JLA bandpasses with sncosmo
    helpers.get_JLA_bandpasses()
    helpers.register_JLA_magsys()
    
    if not os.path.exists('./fit_results'):
        os.makedirs('./fit_results')

    lcfile = './JLA/jla_light_curves/jla_lc.txt'
    lc = np.genfromtxt(lcfile, dtype=None)
    restcut = (3000,7000)

    for i, sn in enumerate(lc):
        if args.specific is not None:
            if sn != args.specific:
                continue
        print '*' * 60
        print sn, i
        
        # get the light curve and if it exists, the covmat
        data = sncosmo.read_lc(sn, format='salt2')
        covmat = sn.replace('lc-', 'covmat_lc-')
        covmat = covmat.replace('.list', '.dat')
        
        if os.path.isfile(covmat):
            covmat = np.loadtxt(covmat, skiprows=1)
            c = Column(covmat,'Fluxcov')
            data.add_column(c)
        else:
            covmat = np.diag(data['Fluxerr']**2)
            c = Column(covmat,'Fluxcov')
            data.add_column(c)

        # deal with the different metadata keywords
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

        #rename columns 
        data.rename_column('Filter', 'tmp')
        data.rename_column('Date', 'time')
        data.rename_column('Flux', 'flux')
        data.rename_column('Fluxerr', 'fluxerr')
        data.rename_column('Fluxcov', 'cov')
        data.rename_column('MagSys', 'zpsys')
        data.rename_column('ZP', 'zp')

        # SNLS need special bandpasses, and we have to make a new column to deal with dtype issues
        if survey == 'SNLS':
            sn_nickname = sn.split('/')[-1].split('.')[0].split('-')[-1]
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
        for ub in unique_bands:
            print ub
            bp = sncosmo.get_bandpass(ub)
            rest = bp.wave_eff / (1.0+z)
            if (rest >= restcut[0]) & (rest <= restcut[1]):
                fit_bands.append(ub)
            else:
                nofit_bands.append(ub)
            if 'STANDARD' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(standard_zps[ub])
            if '4SHOOTER2' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(FourShooter_zps[ub])
            if 'KEPLERCAM' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(keplercam_zps[ub])
            if 'swope' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(swope_zps[ub])
            if 'sdss' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(sdss_zps[ub])

        for nfb in nofit_bands:
            mask = np.array(data['band'] != nfb)
            data = sncosmo.select_data(data,mask)

        # build the normal salt model and the salt2x model
        mwebv = data.meta['MWEBV']
        dust = sncosmo.CCM89Dust()

        if args.emcee:
            if not os.path.exists('./plots/emcee/jla'):
                os.makedirs('./plots/emcee/jla')
            if not os.path.exists('./plots/emcee/jla/triangle'):
                os.makedirs('./plots/emcee/jla/triangle')
            if not os.path.exists('./plots/emcee/jla/salt'):
                os.makedirs('./plots/emcee/jla/salt')

            # make the 2stretch source, apply dust, set it to the right z
            source = Salt2XSource(version='2.4', modeldir=modeldir)
            model  = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])

            SaltSource = sncosmo.SALT2Source(version='2.4', modeldir=modeldir)
            SaltModel = sncosmo.Model(source=SaltSource, effects=[dust], effect_names=['mw'], effect_frames=['obs'])

            model.set(z=z, mwebv=float(mwebv), mwr_v=3.1)
            SaltModel.set(z=z, mwebv=float(mwebv), mwr_v=3.1)

            if args.noskew:
                emfit = fitters.emcee_salt_fit_noskew(data,model,SaltModel)
            else:
                emfit = fitters.emcee_salt_fit(data,model,SaltModel)

            try:
                try:
                    cov, res = emfit.normal_salt_fit(nickname)
                except:
                    traceback.print_exc()
                    print 'normal_salt_fit failed exception'
                    continue

                nsamp = args.nsamp
                t0 = time.time()
                print 'sampling %s times...' %(args.nsamp*emfit.nwalkers)
                fit = emfit.run(nsamples=nsamp)
                print 'time %s samples: %s' %(nsamp*emfit.nwalkers, time.time()-t0)
                
                # make the chains into a dictionary keyed by param name
                chains = emfit.chain_dict(fit)
                
                # occasionally samples go into x0 < 0. Rather than restricting
                # at the likelihood level, we sample longer until we have 
                # enough valid chains
                # in this final implementation of the code I don't think
                # the while loop is ever actually entered for JLA
                ind = np.where(~np.isnan(chains['mB']))
                print 'good samples: %s' %len(chains['mB'][ind])
                if len(chains['mB'][ind]) >= args.nsamp*emfit.nwalkers:
                    keep_going = False
                else:
                    keep_going = True
                ctr = 0
                while keep_going:
                    ctr += 1
                    print 'sampling %i00 times...' %args.nsamp
                    fit = emfit.keep_going(fit, float(mwebv), nsamples=nsamp)
                    print 'time %s samples: %s' %(nsamp, time.time()-t0)
                    chains = emfit.chain_dict(fit)
                    ind = np.where(~np.isnan(chains['mB']))
                    print 'good samples: %s' %len(chains['mB'][ind])
                    if len(chains['mB'][ind]) >= args.nsamp*100:
                        keep_going = False
                    if ctr >= 5:
                        print 'Sampled 5 times. Just stopping'
                        print 'good samples: %s' %len(chains['mB'][ind])
                        keep_going = False
                print 'finally done sampling!'


                # these are just crude errors to be printed in the LC plots
                err = {}
                for k in chains.keys():
                    e = np.percentile(chains[k], [16,84])
                    err[k] = 0.5*(e[1]-e[0])

                # using the maximum likelihood sample just for plotting purposes
                maxlike = np.argmax(chains['lnprob'])
                model.set(t0=chains['t0'][maxlike], x1=chains['x1'][maxlike], s=chains['s'][maxlike], c=chains['c'][maxlike], x0=chains['x0'][maxlike], z=z, mwebv=float(mwebv))

                # the errors passed in here are errors in the measured parameters. Best-fit taken from maximum likelihood sample
                sncosmo.plot_lc(emfit.data, model=model, errors=err, fname='./plots/emcee/jla/%s.pdf' %(nickname),color='black')
                if args.noskew:
                    triangle_keys = ['mB', 'c', 't0', 'x1']
                else:
                    triangle_keys = ['mB', 'c','s', 't0', 'x1']

                helpers.save(chains, './chains/%s.chains' %(nickname))

                # triangle plots
                emfit.plots(chains, nickname, triangle_keys, outdir='./plots/emcee/jla/triangle')
            except Exception as e:
                # as of final release of the code no JLA SNe fail the try except:
                # occasionally a simulated LC will fail, usually due to poor S/N
                traceback.print_exc()
                continue
                
            # output lcfit file
            outdir = os.path.abspath('./fit_results/emcee/JLA/%s' %(nickname))
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                
            if args.noskew:
                params = [chains['x0'][maxlike],chains['x1'][maxlike],chains['c'][maxlike],chains['t0'][maxlike]]
            else:
                params = [chains['x0'][maxlike],chains['x1'][maxlike],chains['s'][maxlike],chains['c'][maxlike],chains['t0'][maxlike]]
                
            chisq = emfit.chi2(params)
            dof = len(emfit.data['flux'])
            print chisq, dof
            #calculate mB

            # dumps statistics from the chains to a file, though the standardization code uses the chain files directly
            helpers.dump_emcee_results(chains, outdir, nickname, z, chisq, dof, survey)

    outfile.close()

if args.cadencesim:
    # The zeropoints in the sim files are already corrected for the magnitude systems registered to sncosmo
    standard_zps = {'STANDARD::U':0, 'STANDARD::B':0, 'STANDARD::V':0, 'STANDARD::R':0, 'STANDARD::I':0}
    FourShooter_zps = {'4SHOOTER2::Us':0, '4SHOOTER2::B':0, '4SHOOTER2::V':0, '4SHOOTER2::R':0, '4SHOOTER2::I':0}
    keplercam_zps = {'KEPLERCAM::Us':0, 'KEPLERCAM::B':0, 'KEPLERCAM::V':0, 'KEPLERCAM::r':0, 'KEPLERCAM::i':0}
    swope_zps = {'SWOPE2::u':0, 'SWOPE2::g':0, 'SWOPE2::r':0, 'SWOPE2::i':0, 'SWOPE2::B':0, 'swope2::v_lc3009':0, 'swope2::v_lc3014':0, 'swope2::v_lc9844':0}
    sdss_zps = {'SDSS::u':0, 'SDSS::g':0, 'SDSS::r':0, 'SDSS::i':0, 'SDSS::z':0}
    
    # registering the JLA bandpasses with sncosmo
    helpers.get_JLA_bandpasses()
    helpers.register_JLA_magsys()
    
    if not os.path.exists('./fit_results'):
        os.makedirs('./fit_results')
    lcfile = './cadence_sim/lc/sim_lc.txt'
    lc = np.genfromtxt(lcfile, dtype=None)
    restcut = (3000,7000)

    for i, sn in enumerate(lc):
        if args.specific is not None:
            if sn != args.specific:
                continue
        print '*' * 60
        print sn, i
        
        # get the light curve
        data = sncosmo.read_lc(sn, format='salt2')

        covmat = np.diag(data['FluxPsferr']**2)
        c = Column(covmat,'Fluxcov')
        data.add_column(c)

        # deal with the different metadata keywords
        try:
            z = data.meta['REDSHIFT']
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
        band_nickname = data.meta['SN']
        nickname = sn.split('.')[0].split('/')[-1]
        try:
            nickname = str(int(nickname))
        except:
            pass
            
        #rename some columns
        data.rename_column('Filter', 'tmp')
        data.rename_column('Date', 'time')
        data.rename_column('FluxPsf', 'flux')
        data.rename_column('FluxPsferr', 'fluxerr')
        data.rename_column('Fluxcov', 'cov')
        data.rename_column('MagSys', 'zpsys')
        data.rename_column('ZP', 'zp')


        # SNLS need special bandpasses, and we have to make a new column to deal with dtype issues
        if survey == 'SNLS':
            sn_nickname = sn.split('/')[-1].split('.')[0].split('-')[-1]
            band = []
            for j, bp in enumerate(data['tmp']):
                band.append( '%s-%s' %(band_nickname, bp) )
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
        for ub in unique_bands:
            bp = sncosmo.get_bandpass(ub)
            rest = bp.wave_eff / (1.0+z)
            if (rest >= restcut[0]) & (rest <= restcut[1]):
                fit_bands.append(ub)
            else:
                nofit_bands.append(ub)
            if 'STANDARD' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(standard_zps[ub])
            if '4SHOOTER2' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(FourShooter_zps[ub])
            if 'KEPLERCAM' in ub:
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(keplercam_zps[ub])
            if 'swope' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(swope_zps[ub])
            if 'sdss' in ub.lower():
                ind = np.where(data['band'] == ub)
                data['zp'][ind] = data['zp'][ind] - float(sdss_zps[ub])

            print ub, bp.wave_eff, rest
            
        for nfb in nofit_bands:
            mask = np.array(data['band'] != nfb)
            data = sncosmo.select_data(data,mask)
            
        mwebv = data.meta['MWEBV']
        dust = sncosmo.CCM89Dust()

        if args.emcee:
            if not os.path.exists('./plots/emcee/cadencesim'):
                os.makedirs('./plots/emcee/cadencesim')
            if not os.path.exists('./plots/emcee/cadencesim/triangle'):
                os.makedirs('./plots/emcee/cadencesim/triangle')
            if not os.path.exists('./plots/emcee/cadencesim/salt'):
                os.makedirs('./plots/emcee/cadencesim/salt')

            # make the 2stretch source, apply dust, set it to the right z
            source = Salt2XSource(version='2.4', modeldir=modeldir)
            model  = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])

            SaltSource = sncosmo.SALT2Source(version='2.4', modeldir=modeldir)
            SaltModel = sncosmo.Model(source=SaltSource, effects=[dust], effect_names=['mw'], effect_frames=['obs'])

            model.set(z=z, mwebv=float(mwebv), mwr_v=3.1)
            SaltModel.set(z=z, mwebv=float(mwebv), mwr_v=3.1)

            if args.noskew:
                emfit = fitters.emcee_salt_fit_noskew(data,model,SaltModel)
            else:
                emfit = fitters.emcee_salt_fit(data,model,SaltModel)

            try:
                try:
                    cov,res = emfit.normal_salt_fit(nickname)
                    data = emfit.data
                    invcov = emfit.invcov

                except:
                    traceback.print_exc()
                    print 'normal_salt_fit failed exception'
                    continue

                nsamp = args.nsamp
                t0 = time.time()
                print 'sampling %i times...' %(nsamp*emfit.nwalkers)
                fit = emfit.run(nsamples=nsamp)
                print 'time %s samples: %s' %(nsamp*emfit.nwalkers, time.time()-t0)
                
                # make the chains into a dictionary keyed by param name
                chains = emfit.chain_dict_x0(fit)
                
                # occasionally samples go into x0 < 0. Rather than restricting
                # at the likelihood level, we sample longer until we have 
                # enough valid chains
                # in this final implementation of the code I don't think
                # the while loop is ever actually entered for JLA
                ind = np.where(~np.isnan(chains['mB']))
                print 'good samples: %s' %len(chains['mB'][ind])
                if len(chains['mB'][ind]) >= args.nsamp*emfit.nwalkers:
                    keep_going = False
                else:
                    keep_going = True
                ctr = 0
                while keep_going:
                    ctr += 1
                    print 'sampling %s times...' %args.nsamp*emfit.nwalkers
                    fit = emfit.keep_going(fit, float(mwebv), nsamples=nsamp)
                    print 'time %s samples: %s' %(nsamp*emfit.nwalkers, time.time()-t0)
                    chains = emfit.chain_dict_x0(fit)
                    ind = np.where(~np.isnan(chains['mB']))
                    print 'good samples: %s' %len(chains['mB'][ind])
                    if len(chains['mB'][ind]) >= args.nsamp*emfit.nwalkers:
                        keep_going = False
                    if ctr >= 5:
                        print 'Sampled 5 times. Just stopping'
                        print 'good samples: %s' %len(chains['mB'][ind])
                        keep_going = False
                print 'finally done sampling!'

                # these are just crude errors to be printed in the LC plots
                err = {}
                for k in chains.keys():
                    e = np.percentile(chains[k], [16,84])
                    err[k] = 0.5*(e[1]-e[0])

                # using the maximum likelihood sample just for plotting purposes
                maxlike = np.argmax(chains['lnprob'])
                model.set(t0=chains['t0'][maxlike], x1=chains['x1'][maxlike], s=chains['s'][maxlike], c=chains['c'][maxlike], x0=chains['x0'][maxlike], z=z, mwebv=float(mwebv))

                # the errors passed in here are errors in the measured parameters. Best-fit taken from maximum likelihood sample
                sncosmo.plot_lc(emfit.data, model=model, errors=err, fname='./plots/emcee/cadencesim/%s.pdf' %(nickname),color='black')
                triangle_keys = ['mB', 'c', 't0', 'x1']

                helpers.save(chains, './chains/%s.chains' %(nickname))

                # triangle plots
                emfit.plots(chains, nickname, triangle_keys, outdir='./plots/emcee/cadencesim/triangle')
            except Exception as e:
                # as of final release of the code no JLA SNe fail the try except:
                # occasionally a simulated LC will fail, usually due to poor S/N
                traceback.print_exc()
                continue
                
            # output lcfit file
            outdir = os.path.abspath('./fit_results/emcee/cadencesim/%s' %(nickname))
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            if args.noskew:
                params = [chains['x0'][maxlike],chains['x1'][maxlike],chains['c'][maxlike],chains['t0'][maxlike]]
            else:
                params = [chains['x0'][maxlike],chains['x1'][maxlike],chains['s'][maxlike],chains['c'][maxlike],chains['t0'][maxlike]]
                
            chisq = emfit.chi2(params)
            dof = len(emfit.data['flux'])
            print chisq, dof
            #calculate mB

            # dumps statistics from the chains to a file, though the standardization code uses the chain files directly
            helpers.dump_emcee_results(chains, outdir, nickname, z, chisq, dof, survey)

    outfile.close()

