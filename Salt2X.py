import numpy as np
from scipy.interpolate import (InterpolatedUnivariateSpline as Spline1d,
                               RectBivariateSpline as Spline2d,
                               splmake, spleval)
# from astropy.utils import OrderedDict as odict
from collections import OrderedDict as odict
from astropy.utils.misc import isiterable
from astropy import (cosmology, units as u, constants as const)
from astropy.extern import six

import sncosmo
from sncosmo.io import read_griddata_ascii
from sncosmo import registry
from sncosmo.bandpasses import get_bandpass, Bandpass
from sncosmo.magsystems import get_magsystem

import abc
import os
from copy import copy as cp
from textwrap import dedent
from math import ceil

class Salt2XSource(sncosmo.SALT2Source):
    _param_names = ['x0','x1','s','c']
    param_names_latex = ['x_0','x_1','s','c']

    def __init__(self, modeldir=None,
                 m0file='salt2_template_0.dat',
                 m1file='salt2_template_1.dat',
                 clfile='salt2_color_correction.dat',
                 cdfile='salt2_color_dispersion.dat',
                 errscalefile='salt2_lc_dispersion_scaling.dat',
                 lcrv00file='salt2_lc_relative_variance_0.dat',
                 lcrv11file='salt2_lc_relative_variance_1.dat',
                 lcrv01file='salt2_lc_relative_covariance_01.dat',
                 name=None, version=None):
        super(Salt2XSource, self).__init__(modeldir=modeldir,
                                            m0file=m0file,
                                            m1file=m1file,
                                            clfile=clfile,
                                            cdfile=cdfile,
                                            errscalefile=errscalefile,
                                            lcrv00file=lcrv00file,
                                            lcrv11file=lcrv11file,
                                            lcrv01file=lcrv01file,
                                            name=name,version=version)

        self._parameters = np.array([1.,0.,0.,0.])

    def _flux(self,phase,wave):
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        sr = self._parameters[2]
        sf = self._parameters[1]
        #transitioning x1 smoothly from -3 to 3 days
        ind = np.where((phase>=-3) & (phase<=3))
        transition = self._smooth_transition(phase[ind], sr, sf)
        x1m1 = np.copy(m1)
        m1_rise = sr*x1m1[phase<=0]
        m1_fall = sf*x1m1[phase>0]
        m1_trans = (x1m1[ind].T*transition).T
        x1m1[phase<=0] = m1_rise
        x1m1[phase>0] = m1_fall
        x1m1[ind] = m1_trans

        return (self._parameters[0] * (m0 + x1m1) *
                10. ** (-0.4 * self._colorlaw(wave) * self._parameters[3]))

    def _smooth_transition(self,phase,sr,sf):
        # b = (sr+sf)/2.0
        # m = (sf-b)/3.0
        # s = m*phase + b
        s = sr + 0.5*(np.sin(phase*np.pi/6.)+1)*(sf-sr)
        return s
        