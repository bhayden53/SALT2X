#! /bin/bash
# total samples = NSAMP*nwalkers, which is 100 in the code.
# at a minimum nsamp*nwalkers should really be 100,000 (in the paper we use 500,000) 
# this becomes computationally expensive, so salt_2s.py has a
# --specific keyword that takes the full path to a lightcurve file.
# This allows each SN to be run in it's own thread, i.e. submitting
# to a batch job system on a typical compute cluster

NSAMP=10

cd cadence_sim/lc
find `pwd` -name "*.list" > sim_lc.txt
cd ../..
# the --noskew flag performs a normal salt fit. Remove the flag for a SALT2X run
python salt_2s.py --emcee --cadencesim --noskew --nsamp $NSAMP
python combine_emcee_fits.py --cadencesim