#!/bin/bash
# will generate a new set of simulated light curves. Wipes out any existing set!
rm cadence_sim/lc/*
python mklc.py
cd cadence_sim/lc
find `pwd` -name "*.list" > sim_lc.txt
cd ../
python scrape_input_vals.py