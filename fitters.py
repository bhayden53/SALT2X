#!/usr/bin/env python

# from emcee.utils import MPIPool
import numpy as np
import corner as triangle
import emcee
from iminuit import Minuit
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from sncosmo import photdata
import copy
import sncosmo
from IPython import embed
import sys
import copy
import pickle
import helpers
from copy import deepcopy
import os

scratch = os.environ['SCRATCH']

#note: the skew/noskew settings require some hand tuning. it doesn't work all by itself

class emcee_salt_fit(object):

    def __init__(self,data,model,SaltModel):
        self.data = data

        test_data = photdata.photometric_data(deepcopy(data))

        self.model = model
        self.z = self.model.get('z')
        self.bands = np.unique(self.data['band'])
        self.x1_prior = (0,np.sqrt(2))
        self.s_prior = (0,np.sqrt(2))
        self.c_prior = (0,0.15)
        self.SaltModel = SaltModel
        self.ndim = 5
        self.nwalkers = 100
        self.tmax_guess, self.x0_start = sncosmo.fitting.guess_t0_and_amplitude(test_data, self.model, minsnr=3)
        self.tmax_bounds = sncosmo.fitting.t0_bounds(test_data, self.model)

        # try guess tmax ourselves
        t = np.arange(self.tmax_bounds[0], self.tmax_bounds[1]+1, 1)
        test_chi = 1e20
        test_time = 0
        test_x0 = 0
        x0arr = []


        for time in t:
            params = [self.x0_start, 0,0,0, time]
            self.SaltModel.set(t0=time,x0=self.x0_start,x1=0,c=0, z=self.z)
            res, fitted_model = sncosmo.fit_lc(deepcopy(self.data), self.SaltModel, ['x0'], guess_amplitude=False, guess_t0=False, modelcov=False)
            c = res.chisq

            x0arr.append(res.parameters[2])
            if np.isnan(c) or np.isinf(c):
                continue
            if c < test_chi:
                test_chi = c
                test_time = time
                test_x0 = res.parameters[2]

        self.tmax_guess = test_time
        self.x0_start = test_x0
        self.nest_bounds = {'t0':self.tmax_bounds, 'x1':[-4,4], 'c':[-1,1], 'x0':[np.min(x0arr), np.max(x0arr)]}

    def normal_salt_fit(self, nickname):

        bounds = {'t0':self.tmax_bounds}
        exception_bounds = {'t0':self.tmax_bounds, 'x1':[-4,4], 'c':[-1,1]}

        self.SaltModel.set(t0=self.tmax_guess,x0=self.x0_start)
        res, fitted_model = sncosmo.fit_lc(deepcopy(self.data), self.SaltModel, ['t0','x0','x1','c'], bounds=exception_bounds, guess_amplitude=False, guess_t0=False, modelcov=True)

        sncosmo.plot_lc(self.data, model=fitted_model, fname='%s/plots/emcee/cadencesim/salt/%s.pdf' %(scratch,nickname),color='black')

        x1 = fitted_model.get('x1')
        c = fitted_model.get('c')


        if (x1<exception_bounds['x1'][0]) | (x1>exception_bounds['x1'][1]) | (c<exception_bounds['c'][0]) | (c>exception_bounds['c'][1]):
            print 'in exception bounds if for first data phase cut iteration'
            self.SaltModel.set(t0=self.tmax_guess,x0=self.x0_start)
            res, fitted_model = sncosmo.fit_lc(deepcopy(self.data), self.SaltModel, ['t0','x0','x1','c'], bounds=exception_bounds, guess_amplitude=False, guess_t0=False, modelcov=True)


        tmax              = fitted_model.get('t0')
        phase             = (self.data['time']-tmax)/(1+self.z)
        phase_mask = np.array(((phase>-15) & (phase<45)))

        self.data = sncosmo.select_data(self.data, phase_mask)

        self.SaltModel.set(t0=self.tmax_guess,x0=self.x0_start)
        res, fitted_model = sncosmo.fit_lc(deepcopy(self.data), self.SaltModel, ['t0','x0','x1','c'], bounds=exception_bounds, guess_amplitude=False, guess_t0=False, modelcov=True)

        x1 = fitted_model.get('x1')
        c = fitted_model.get('c')

        if (x1<exception_bounds['x1'][0]) | (x1>exception_bounds['x1'][1]) | (c<exception_bounds['c'][0]) | (c>exception_bounds['c'][1]):
            raise 

        tmax              = fitted_model.get('t0')
        # tmax              = tmax
        phase             = (self.data['time']-tmax)/(1+self.z)
        phase_mask        = np.array(((phase>-15) & (phase<45)))
        self.data = sncosmo.select_data(self.data, phase_mask)

        # pull out the cov
        self.whocares, self.SaltCov = fitted_model.bandfluxcov(self.data['band'], self.data['time'], self.data['zp'], self.data['zpsys'])
        sncosmo.plot_lc(self.data, model=fitted_model, fname='%s/plots/emcee/cadencesim/salt/%s.pdf' %(scratch,nickname),color='black')

        self.invcov    = np.linalg.inv(self.SaltCov + self.data['cov'])

        self.x0_salt   = fitted_model.get('x0')
        self.tmax_salt = fitted_model.get('t0')
        self.c_salt    = fitted_model.get('c')

        return self.SaltCov, res

    def mcmc_salt_fit(self, nickname):
        # do 2 SALT fits to iterate on the data to use for the fit
        bounds = {'t0':self.tmax_bounds}
        try:
            res, fitted_model = sncosmo.mcmc_lc(self.data, self.SaltModel, ['t0','x0','x1','c'], minsnr=3)
        except:
            # embed()
            print 'mcmc_salt_fit 1st iteration failed'
            self.SaltModel.set(t0=self.tmax_guess,x0=self.x0_start)
            res, fitted_model = sncosmo.mcmc_lc(self.data, self.SaltModel, ['t0','x0','x1','c'], bounds=bounds, guess_amplitude=False, guess_t0=False, modelcov=True)
        tmax              = fitted_model.get('t0')
        phase             = (self.data['time']-tmax)/(1+self.z)
        phase_mask = np.array(((phase>-15) & (phase<45)))
        self.data         = sncosmo.select_data(self.data, phase_mask)

        # fit 2
        try:
            res, fitted_model = sncosmo.mcmc_lc(self.data, self.SaltModel, ['t0','x0','x1','c'], minsnr=3)
        except:
            print 'mcmc_salt_fit 2nd iteration failed'
            raise #nothing succeeded!
            self.SaltModel.set(t0=self.tmax_guess,x0=self.x0_start)
            res, fitted_model = sncosmo.mcmc_lc(self.data, self.SaltModel, ['t0','x0','x1','c'], bounds=bounds, guess_amplitude=False, guess_t0=False, modelcov=True)
        tmax              = fitted_model.get('t0')
        phase             = (self.data['time']-tmax)/(1+self.z)
        phase_mask = np.array(((phase>-15) & (phase<45)))
        self.data         = sncosmo.select_data(self.data, phase_mask)


        # pull out the cov
        self.whocares, self.SaltCov = fitted_model.bandfluxcov(self.data['band'], self.data['time'], self.data['zp'], self.data['zpsys'])
        # plot
        sncosmo.plot_lc(self.data, model=fitted_model, fname='%s/plots/emcee/cadencesim/salt/%s.pdf' %(scratch,nickname),color='black')

        # pull out the cov and some other parameters
        self.invcov    = np.linalg.inv(self.SaltCov + self.data['fluxcov'])
        self.x0_salt   = fitted_model.get('x0')
        self.tmax_salt = fitted_model.get('t0')
        self.c_salt    = fitted_model.get('c')
        chisq = sncosmo.chisq(self.data,fitted_model)
        res.chisq = chisq
        res.ndof = len(self.data) - 4

        return self.SaltCov, res

    def chi2(self,params):
        x0 = params[0]
        x1 = params[1]
        s  = params[2]
        c  = params[3]
        t0 = params[4]

        chis=0

        self.model.set(x0=x0,x1=x1,s=s,c=c, t0=t0, z=self.z)
        model_flux = self.model.bandflux(self.data['band'], self.data['time'], self.data['zp'], self.data['zpsys'])
        residuals  = (model_flux - self.data['flux'])

        chis += np.dot(residuals, np.dot(self.invcov, residuals))

        if np.isnan(chis) or np.isinf(chis):
            chis = np.inf

        return chis

    def lnprior(self, params):
        x0 = params[0]
        x1 = params[1]
        s  = params[2]
        c  = params[3]
        t0 = params[4]

        if (t0 < self.tmax_bounds[0]) | (t0 > self.tmax_bounds[1]):
            return -np.inf
        else:
            return 0


    def loglike(self,params):
        # print params
        lp = self.lnprior(params)
        if not np.isfinite(lp):
            return -np.inf

        ll = -0.5*self.chi2(params)

        if np.isnan(ll) or np.isinf(ll):
            print 'll nan'
            print params
            return -np.inf
        else:
            # print params, ll
            return ll

    def param_names(self):
        return ['x0','x1','s','c','t0']

    def chain_dict(self, sampler):
        p = self.param_names()
        d = {}
        for i,v in enumerate(p):
            d[v] = sampler.flatchain[:,i]
        d['lnprob'] = sampler.flatlnprobability
        # calculate mB
        d['mB'] = []

        for i, x0 in enumerate(d['x0']):
            #negative x0 is no bueno but should be very rare
            if x0 < 0:
                d['mB'].append(30)
            else:
                self.model.set(x0=x0, s=d['s'][i], x1=d['x1'][i], c=d['c'][i], t0=0, z=0, mwebv=0)
                # B = self.model.bandmag('bessellb', 'vega2', 0) + 9.907
                B = self.model.source.peakmag('bessellb', 'vega2', sampling=1.0) + 9.907
                d['mB'].append(B)
        d['mB'] = np.array(d['mB'])
        ind = np.where(~np.isnan(d['mB']))
        for par in p:
            d[par] = d[par][ind]
        d['lnprob'] = d['lnprob'][ind]
        d['mB'] = d['mB'][ind]
        return d

    def chain_dict_x0(self, sampler):
        p = self.param_names()
        d = {}
        for i,v in enumerate(p):
            d[v] = sampler.flatchain[:,i]
        d['lnprob'] = sampler.flatlnprobability
        # calculate mB
        d['mB'] = -2.5*np.log10(d['x0'])
        ind = np.where(d['x0'] < 0)
        d['mB'][ind] = 30
        ind = np.where(~np.isnan(d['mB']))
        for par in p:
            d[par] = d[par][ind]
        d['lnprob'] = d['lnprob'][ind]
        d['mB'] = d['mB'][ind]

        return d

    def run(self, nsamples=5000):

        nwalkers = 100

        # skew
        guess = [self.x0_salt,0,0,self.c_salt,self.tmax_salt ]
        guess = [self.x0_salt,0,0,0,self.tmax_salt ]
        print guess
        steps = [ 0.5*self.x0_start, 0.5, 0.5, 0.1, 2 ]
        randarr = np.random.rand(self.ndim * nwalkers).reshape((nwalkers, self.ndim))
        start = np.array( guess ) + np.array( steps ) * ( randarr - 0.5 )
        sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.loglike, threads=1)

        # burn
        pos, prob, state = sampler.run_mcmc(start, nsamples/2)
        sampler.reset()

        pos,prob,state = sampler.run_mcmc(pos, nsamples)

        self.pos = pos
        self.sampler = sampler
        return sampler



    def keep_going(self, sampler, mwebv, nsamples=500):
        self.model.set(z=self.z, mwebv=mwebv)
        pos,prob,state = sampler.run_mcmc(self.pos, nsamples)
        return sampler


    def plots(self, chains, cid, keys, outdir='./plots/emcee/triangle'):
        #sampler.chain.shape = (walkers, samples, ndim)
        #after reshape, shape = (walkers*sample, ndim)
        tmp = []
        for k in keys:
            tmp.append(chains[k])
        tmp = np.array(tmp)
        tmp = tmp.T
        # embed()
        # samples = sampler.chain[:,:,:].reshape((-1,self.ndim))
        fig = triangle.corner(tmp,labels=keys)
        plt.savefig('%s/%s_tri.pdf' %(outdir,cid ))
        # plt.show()

class emcee_salt_fit_noskew(emcee_salt_fit):
    def __init__(self,data,model,SaltModel):
        super( emcee_salt_fit_noskew, self).__init__(data,model,SaltModel)
        self.ndim = 4

    def chi2(self, params):
        x0 = params[0]
        x1 = params[1]
        s  = x1
        c  = params[2]
        t0 = params[3]

        chis=0
        self.model.set(x0=x0,x1=x1,s=x1,c=c, t0=t0)

        model_flux = self.model.bandflux(self.data['band'], self.data['time'], self.data['zp'], self.data['zpsys'])
        residuals = (model_flux - self.data['flux'])
        chis += np.dot(residuals, np.dot(self.invcov, residuals))

        if np.isnan(chis) or np.isinf(chis):
            chis = np.inf
        return chis

    def lnprior(self, params):
        x0 = params[0]
        x1 = params[1]
        s  = x1
        c  = params[2]
        t0 = params[3]

        if (t0 < self.tmax_bounds[0]) | (t0 > self.tmax_bounds[1]):
            return -np.inf
        else:
            return 0

    def loglike(self,params):
        # print params
        lp = self.lnprior(params)
        if not np.isfinite(lp):
            return -np.inf

        ll = -0.5*self.chi2(params)

        if np.isnan(ll) or np.isinf(ll):
            print 'll nan'
            print params
            return -np.inf
        else:
            # print params, ll
            return ll

    def param_names(self):
        return ['x0','x1','c','t0']

    def chain_dict(self, sampler):
        p = self.param_names()
        d = {}
        for i,v in enumerate(p):
            d[v] = sampler.flatchain[:,i]
        d['lnprob'] = sampler.flatlnprobability
        # calculate mB
        d['mB'] = []
        d['s'] = d['x1']
        for i, x0 in enumerate(d['x0']):
            #negative x0 is no bueno but should be very rare
            if x0 < 0:
                d['mB'].append(30)
            else:
                self.model.set(x0=x0, s=d['s'][i], x1=d['x1'][i], c=d['c'][i], t0=0, z=0, mwebv=0)
                # B = self.model.bandmag('bessellb', 'vega2', 0) + 9.907
                B = self.model.source.peakmag('bessellb', 'vega2', sampling=1.0) + 9.907
                d['mB'].append(B)
        d['mB'] = np.array(d['mB'])
        ind = np.where(~np.isnan(d['mB']))
        for par in p:
            d[par] = d[par][ind]
        d['s'] = d['s'][ind]
        d['lnprob'] = d['lnprob'][ind]
        d['mB'] = d['mB'][ind]
        return d

    def chain_dict_x0(self, sampler):
        p = self.param_names()
        d = {}
        for i,v in enumerate(p):
            d[v] = sampler.flatchain[:,i]
        d['s'] = d['x1']
        d['lnprob'] = sampler.flatlnprobability
        # calculate mB
        d['mB'] = -2.5*np.log10(d['x0'])
        ind = np.where(d['x0'] < 0)
        d['mB'][ind] = 30
        ind = np.where(~np.isnan(d['mB']))
        for par in p:
            d[par] = d[par][ind]
        d['lnprob'] = d['lnprob'][ind]
        d['mB'] = d['mB'][ind]

        return d

    def run(self, nsamples=5000):

        nwalkers = 100

        #noskew
        guess = [self.x0_salt,0,self.c_salt,self.tmax_salt ]
        print guess
        steps = [ 0.5*self.x0_start, 0.5, 0.1, 2 ]
        randarr = np.random.rand(self.ndim * nwalkers).reshape((nwalkers, self.ndim))
        start = np.array( guess ) + np.array( steps ) * ( randarr - 0.5 )
        sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.loglike, threads=1)

        # burn
        pos, prob, state = sampler.run_mcmc(start, nsamples/2)
        sampler.reset()

        pos,prob,state = sampler.run_mcmc(pos, nsamples)

        self.pos = pos
        self.sampler = sampler
        return sampler


    def keep_going(self, sampler, mwebv, nsamples=500):
        self.model.set(z=self.z, mwebv=mwebv)
        pos,prob,state = sampler.run_mcmc(self.pos, nsamples)
        return sampler
