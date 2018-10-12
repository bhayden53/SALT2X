import sncosmo
from astropy.table import Table
import numpy as np
import glob
from IPython import embed
import sys

sne = glob.glob('./lc/*.list')

names = ['SN','tmax', 'color', 'x1r', 'x1f', 'mB', 'mwebv']
dtype = ['S50',np.float64,np.float64,np.float64,np.float64,np.float64, np.float64]
t = Table(names=names, dtype=dtype)

for i, sn in enumerate(sne):
    data = sncosmo.read_lc(sn, format='salt2')
    daymax = data.meta['DAYMAX']
    x1r = data.meta['X1R']
    x1f = data.meta['X1F']
    c = data.meta['C']
    mB = data.meta['MB']
    name = sn.split('.')[-2].split('/')[-1]
    mwebv = data.meta['MWEBV']
    t.add_row([name,daymax,c,x1r,x1f,mB,mwebv])

t.write('./cadence_sim_inputs.txt', format='ascii.fixed_width', delimiter=' ', overwrite=True)
