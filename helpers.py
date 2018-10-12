import numpy as np
import sncosmo
import glob
from astropy import units as u
import os
import gzip
import cPickle

def dump_emcee_results(chains, outfile, sn_name, z, chisq, dof, survey, maxll=True):
    res = {}
    for k in chains.keys():
        res[k] = np.percentile(chains[k], [16,50,84])
    if maxll:
        maxlike = np.argmax(chains['lnprob'])
        for k in chains.keys():
            # print k, res[k][1]
            res[k][1] = chains[k][maxlike]
            # print k, res[k][1]

    cov_t0_t0 = np.cov(chains['t0'], chains['t0'])
    cov_x0_x0 = np.cov(chains['x0'], chains['x0'])
    cov_mB_mB = np.cov(chains['mB'], chains['mB'])
    cov_x1_x1 = np.cov(chains['x1'], chains['x1'])
    cov_s_s   = np.cov(chains['s'], chains['s'])
    cov_c_c   = np.cov(chains['c'], chains['c'])
    cov_t0_x0 = np.cov(chains['t0'], chains['x0'])
    cov_t0_mB = np.cov(chains['t0'], chains['mB'])
    cov_t0_x1 = np.cov(chains['t0'], chains['x1'])
    cov_t0_s = np.cov(chains['t0'], chains['s'])
    cov_t0_c = np.cov(chains['t0'], chains['c'])
    cov_x0_x1 = np.cov(chains['x0'], chains['x1'])
    cov_x0_s = np.cov(chains['x0'], chains['s'])
    cov_x0_c = np.cov(chains['x0'], chains['c'])
    cov_mB_x1 = np.cov(chains['mB'], chains['x1'])
    cov_mB_s = np.cov(chains['mB'], chains['s'])
    cov_mB_c = np.cov(chains['mB'], chains['c'])
    cov_x1_s = np.cov(chains['x1'], chains['s'])
    cov_x1_c = np.cov(chains['x1'], chains['c'])
    cov_s_c = np.cov(chains['s'], chains['c'])


    file = open('%s/%s.dat' %(outfile, sn_name), 'w')
    file.write('Salt2Model\n')
    file.write('BEGIN_OF_FITPARAMS Salt2Model\n')
    file.write('CHISQ %s -1\n' %chisq)
    file.write('DoF %s -1\n' %dof)
    file.write('survey %s -1\n' %survey)
    file.write('DayMax %s %s\n' %(res['t0'][1], 0.5*(res['t0'][2]-res['t0'][0])))
    file.write('Redshift %s 0 F\n' %(z))
    file.write('Color %s %s\n' %(res['c'][1], 0.5*(res['c'][2]-res['c'][0])))
    file.write('X0 %s %s\n' %(res['x0'][1], 0.5*(res['x0'][2]-res['x0'][0])))
    file.write('RestFrameMag_0_B %s %s\n' %(res['mB'][1], 0.5*(res['mB'][2]-res['mB'][0])))
    file.write('X1 %s %s\n' %(res['x1'][1], 0.5*(res['x1'][2]-res['x1'][0])))
    file.write('Skew %s %s\n' %(res['s'][1], 0.5*(res['s'][2]-res['s'][0])))
    file.write('CovDayMaxDayMax %s -1\n' %(cov_t0_t0[0,0]))
    file.write('CovDayMaxDayMax2 %s -1\n' %((0.5*(res['t0'][2]-res['t0'][0]))**2))
    file.write('CovX0X0 %s -1\n' %(cov_x0_x0[0,0]))
    file.write('CovX0X02 %s -1\n' %((0.5*(res['x0'][2]-res['x0'][0]))**2))
    file.write('CovX1X1 %s -1\n' %(cov_x1_x1[0,0]))
    file.write('CovX1X12 %s -1\n' %((0.5*(res['x1'][2]-res['x1'][0]))**2))
    file.write('CovSkewSkew %s -1\n' %(cov_s_s[0,0]))
    file.write('CovSkewSkew2 %s -1\n' %((0.5*(res['s'][2]-res['s'][0]))**2))
    file.write('CovColorColor %s -1\n' %(cov_c_c[0,0]))
    file.write('CovColorColor2 %s -1\n' %((0.5*(res['c'][2]-res['c'][0]))**2))
    file.write('CovRestFrameMag_0_BRestFrameMag_0_B %s -1\n' %(cov_mB_mB[0,0]))
    file.write('CovRestFrameMag_0_BRestFrameMag_0_B2 %s -1\n' %(0.5*(res['mB'][2]-res['mB'][0]))**2)
    file.write('CovDayMaxX0 %s -1\n' %(cov_t0_x0[0,1]))
    file.write('CovDayMaxRestFrameMag_0_B %s -1\n' %(cov_t0_mB[0,1]))
    file.write('CovDayMaxX1 %s -1\n' %(cov_t0_x1[0,1]))
    file.write('CovDayMaxSkew %s -1\n' %(cov_t0_s[0,1]))
    file.write('CovDayMaxColor %s -1\n' %(cov_t0_c[0,1]))
    file.write('CovX0X1 %s -1\n' %(cov_x0_x1[0,1]))
    file.write('CovX0Skew %s -1\n' %(cov_x0_s[0,1]))
    file.write('CovX0Color %s -1\n' %(cov_x0_c[0,1]))
    file.write('CovRestFrameMag_0_BX1 %s -1\n' %(cov_mB_x1[0,1]))
    file.write('CovRestFrameMag_0_BSkew %s -1\n' %(cov_mB_s[0,1]))
    file.write('CovRestFrameMag_0_BColor %s -1\n' %(cov_mB_c[0,1]))
    file.write('CovX1Skew %s -1\n' %(cov_x1_s[0,1]))
    file.write('CovX1Color %s -1\n' %(cov_x1_c[0,1]))
    file.write('CovSkewColor %s -1\n' %(cov_s_c[0,1]))
    # percentiles for evaluating gaussianity of chains
    file.write('DayMax_p2 %s -1\n' %(np.percentile(chains['t0'], 2.27501)))
    file.write('DayMax_p15 %s -1\n' %(np.percentile(chains['t0'], 15.8655)))
    file.write('DayMax_p50 %s -1\n' %(np.percentile(chains['t0'], 50)))
    file.write('DayMax_p84 %s -1\n' %(np.percentile(chains['t0'], 84.1345)))
    file.write('DayMax_p97 %s -1\n' %(np.percentile(chains['t0'], 97.725)))
    file.write('x0_p2 %s -1\n'  %(np.percentile(chains['x0'], 2.27501)))
    file.write('x0_p15 %s -1\n' %(np.percentile(chains['x0'], 15.8655)))
    file.write('x0_p50 %s -1\n' %(np.percentile(chains['x0'], 50)))
    file.write('x0_p84 %s -1\n' %(np.percentile(chains['x0'], 84.1345)))
    file.write('x0_p97 %s -1\n' %(np.percentile(chains['x0'], 97.725)))
    file.write('mB_p2 %s -1\n'  %(np.percentile(chains['mB'], 2.27501)))
    file.write('mB_p15 %s -1\n' %(np.percentile(chains['mB'], 15.8655)))
    file.write('mB_p50 %s -1\n' %(np.percentile(chains['mB'], 50)))
    file.write('mB_p84 %s -1\n' %(np.percentile(chains['mB'], 84.1345)))
    file.write('mB_p97 %s -1\n' %(np.percentile(chains['mB'], 97.725)))
    file.write('x1_p2 %s -1\n'  %(np.percentile(chains['x1'], 2.27501)))
    file.write('x1_p15 %s -1\n' %(np.percentile(chains['x1'], 15.8655)))
    file.write('x1_p50 %s -1\n' %(np.percentile(chains['x1'], 50)))
    file.write('x1_p84 %s -1\n' %(np.percentile(chains['x1'], 84.1345)))
    file.write('x1_p97 %s -1\n' %(np.percentile(chains['x1'], 97.725)))
    file.write('s_p2 %s -1\n'  %(np.percentile(chains['s'], 2.27501)))
    file.write('s_p15 %s -1\n' %(np.percentile(chains['s'], 15.8655)))
    file.write('s_p50 %s -1\n' %(np.percentile(chains['s'], 50)))
    file.write('s_p84 %s -1\n' %(np.percentile(chains['s'], 84.1345)))
    file.write('s_p97 %s -1\n' %(np.percentile(chains['s'], 97.725)))
    file.write('Color_p2 %s -1\n'  %(np.percentile(chains['c'], 2.27501)))
    file.write('Color_p15 %s -1\n' %(np.percentile(chains['c'], 15.8655)))
    file.write('Color_p50 %s -1\n' %(np.percentile(chains['c'], 50)))
    file.write('Color_p84 %s -1\n' %(np.percentile(chains['c'], 84.1345)))
    file.write('Color_p97 %s -1\n' %(np.percentile(chains['c'], 97.725)))

    file.close()
    # print chains.keys()

def get_JLA_bandpasses():
    bd17 = sncosmo.get_magsystem('bd17')
    vega = sncosmo.get_magsystem('vega')

    sndata_root = '%s/filters/JLA' %(os.environ['SNDATA_ROOT'])
    filter_dir = os.environ['JLA_FILTERDIR']

    # 4shooter2
    files = glob.glob('%s/4Shooter2/*.txt' %(filter_dir))
    for f in files:
        fspl = f.split('/')[-1].split('.')[0].split('_')
        # print fspl
        # if fspl[0] == 'U':
            # fspl[0] = 'us'
        fstr = '%s::%s' %(fspl[1].lower(),fspl[0])
        d = np.genfromtxt(f, dtype=None, skip_header=1)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)

    #HST
    instmap = {'F110W':'NICMOS2', 'F160W':'NICMOS2', 'F606W':'ACSWF', 'F625W':'ACSWF',
                'F775W':'ACSWF', 'F814W':'ACSWF', 'F850LP':'ACSWF'}

    #NICMOS
    files = glob.glob('%s/NICMOS2/*.dat' %(filter_dir))
    for f in files:
        fspl = f.split('/')[-1].split('.')[0]
        # print fspl
        fstr = '%s::%s' %(instmap[fspl], fspl)
        d = np.genfromtxt(f,dtype=None)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)

    #ACS
    files = glob.glob('%s/ACSWF/*.dat' %(filter_dir))
    for f in files:
        fspl = f.split('/')[-1].split('.')[0]
        # print fspl
        fstr = '%s::%s' %(instmap[fspl], fspl)
        d = np.genfromtxt(f,dtype=None)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)

    #Keplercam
    files = glob.glob('%s/Keplercam/*.txt*' %(filter_dir))
    for f in files:
        fspl = f.split('/')[-1].split('.')[0].split('_')
        if fspl[0] == 'U':
            fstr = '%s::%s' %(fspl[1].lower(),'Us' )
            d = np.genfromtxt(f,dtype=None)
            band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
            sncosmo.registry.register(band, fstr, force=True)
        else:
            fstr = '%s::%s' %(fspl[1].lower(),fspl[0] )
            d = np.genfromtxt(f,dtype=None, skip_header=1)
            band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
            sncosmo.registry.register(band, fstr, force=True)

    #SNLS Megacam
    files = glob.glob('%s/JLA-Megacam/filters*' %(sndata_root))
    dtype = np.dtype([('f0', '<f8'), ('f1', '<f8')])
    for f in files:
        # print '%s/JLA-Megacam/%s/*.dat' %(filter_dir,f)
        files2 = glob.glob('%s/*.dat' %(f))
        for f2 in files2:
            fspl = f2.split('/')
            sn = fspl[-2].split('-')[1]
            fil = fspl[-1].split('.')[0].split('-')
            # filname = fspl[-1].split('.')[0].split('-')[-1]
            fstr = '%s-%s::%s' %(sn, fil[-2], fil[-1])
            d = np.genfromtxt(f2, dtype=dtype, skip_header=5)
            band = sncosmo.Bandpass(d['f0'], d['f1'], name=fstr)
            sncosmo.registry.register(band, fstr, force=True)
            # print d.dtype
            # print f2
            # print d['f0']

    #SDSS
    files = glob.glob('%s/SDSS/?.dat' %filter_dir)
    for f in files:
        fspl = f.split('/')[-1].split('.')[0]
        fstr = 'SDSS::%s' %(fspl)
        d = np.genfromtxt(f, dtype=None, skip_header=1)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)

    #landolt
    files = glob.glob('%s/SNLS3-Landolt-model/*.dat' %(filter_dir))
    for f in files:
        if 'shifted' not in f:
            continue
        fspl = f.split('/')[-1].split('.')[0].split('-')[0][1:2]
        fstr = 'STANDARD::%s' %(fspl.upper())
        d = np.genfromtxt(f, dtype=None)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)
        # print fspl, bd17.zpbandflux(fstr), vega.zpbandflux(fstr)


    #swope
    files = glob.glob('%s/Swope/*.txt' %(os.environ['DAVID_FILTERDIR']))
    for f in files:
        fspl = f.split('/')[-1].split('_texas')
        fstr = 'SWOPE2::%s' %(fspl[0])
        d = np.genfromtxt(f, dtype=None)
        band = sncosmo.Bandpass(d[:,0], d[:,1], name=fstr)
        sncosmo.registry.register(band, fstr, force=True)
        # print fspl, bd17.zpbandflux(fstr)


        # print instmap[fspl]

    # embed()

def register_JLA_magsys():
    #AB_B12
    spec = np.genfromtxt('%s/MagSys/ab-spec.dat' %(os.environ['SNFIT_DATA']))
    wave = spec[:,0]
    flux = spec[:,1]
    specsnc = sncosmo.Spectrum(wave,flux,unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    magsys = sncosmo.magsystems.SpectralMagSystem(specsnc)
    sncosmo.registry.register(magsys, 'AB_B12', data_class=sncosmo.MagSystem, force=True)
    test = sncosmo.get_magsystem('ab_b12')
    # embed()
    print test

    #vega2
    spec = np.genfromtxt('%s/MagSys/bd_17d4708_stisnic_003.ascii' %(os.environ['SNFIT_DATA']))
    wave = spec[:,0]
    flux = spec[:,1]
    specsnc = sncosmo.Spectrum(wave,flux,unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    magsys = sncosmo.magsystems.SpectralMagSystem(specsnc)
    sncosmo.registry.register(magsys, 'VEGA2', data_class=sncosmo.MagSystem, force=True)

    #vegahst
    spec = np.genfromtxt('%s/MagSys/alpha_lyr_stis_005.ascii' %(os.environ['SNFIT_DATA']))
    wave = spec[:,0]
    flux = spec[:,1]
    specsnc = sncosmo.Spectrum(wave,flux,unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    magsys = sncosmo.magsystems.SpectralMagSystem(specsnc)
    sncosmo.registry.register(magsys, 'VEGAHST', data_class=sncosmo.MagSystem, force=True)
    
def save(object, filename, protocol=-1):
    """Saves a compressed object to disk
	"""
    file = gzip.GzipFile(filename, 'wb')
    cPickle.dump(object, file, protocol)
    file.close()

def load( filename ):
    """Loads a compressed object from disk
    """
    file = gzip.GzipFile(filename, 'rb')
    object = cPickle.load( file )
    file.close()
    return object

if __name__ == '__main__':
    get_JLA_bandpasses()
    register_JLA_magsys()




