# SALT2X
The SALT2X light curve fitting method implemented in SNCOSMO

The JLA light curves should be placed in a directory called "JLA/jla_light_curves".
In that directory should be a list of full paths to each light curve called 
jla_lc.txt. You can create this file by calling:
`` find `pwd` -name "*.list" > jla_lc.txt ``
from the jla_light_curves directory.

You will need SNDATA_ROOT from SNANA installed and in your environment. SALT2X
gets the SNLS SN-specific filter functions from SNANA.
http://snana.uchicago.edu/

There are shell scripts for producing the simulated data and running SALT2X. 
The scripts to run SALT2X are examples only; in practice, each SN is run 
individually with 500,000 samples (after 250,000 burn in) using the --specific 
command line option and is submitted as a single job in a batch submission 
system. 