from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
from astroML.plotting import hist as aml_hist
from scipy.optimize import curve_fit

# lcfits = Table.read('../emcee_cadencesim_skewfit_021417.dat', format='ascii')
# lcfits = Table.read('../fit_backups/emcee_skewfit_cadencesim_PAPERSIM_3_fits_040918.dat', format='ascii')
lcfits = Table.read('../test.dat', format='ascii')
inputs = Table.read('./cadence_sim_inputs.txt',format='ascii')

joined = join(lcfits, inputs, join_type='inner')

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)


# Get the fitted curve
ind = np.where((joined['Skew_coverr'] < 0.7) & (joined['X1_coverr'] < 0.7) & (joined['Color_coverr'] < 0.3))
A = len(joined['color'][ind])
print A
p0 = [A, 0., 1.]

### Color ###
plt.clf()
plt.errorbar(joined['color'][ind], joined['Color'][ind], yerr=joined['Color_coverr'][ind], fmt='o', alpha=0.5)
plt.plot([-0.3,0.3], [-0.3,0.3], color='k', lw=2.5)
plt.xlim(-0.32,0.32)
plt.ylim(-0.32,0.32)
plt.xlabel('input color')
plt.ylabel('measured color (SkewFit)')
plt.savefig('./plots/color_v_color_skewfit.pdf')
# plt.show()

plt.clf()
plt.hist(joined['Color'][ind] - joined['color'][ind], weights=joined['Color_coverr'][ind], bins=25)
plt.xlabel('Color(skewfit) - Color(input)')
plt.ylabel('weighted histogram')
plt.savefig('./plots/color_v_color_skewfit_hist.pdf')
# plt.show()

# color pulls
plt.clf()
# pulls
pulls = (joined['Color'][ind] - joined['color'][ind])/joined['Color_coverr'][ind]
#pull hist
hist, bin_edges, dummy = aml_hist(pulls, bins='scott', alpha=0.25)
#gaussian fit and plot
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
x = np.linspace(np.min(pulls), np.max(pulls), 150)
y = gauss(x, *coeff)
plt.plot(x,y, label='Fitted data')
#other plotting stuff
plt.annotate('mu: %s'%np.round(coeff[1],4), xy=(0.6,0.8), xycoords='axes fraction')
plt.annotate('sigma: %s'%np.round(coeff[2],4), xy=(0.6,0.7), xycoords='axes fraction')
plt.xlabel('Color pulls')
plt.ylabel('histogram')
plt.xlim(-5,5)
plt.savefig('./plots/color_v_color_skewfit_pulls.pdf')
# plt.show()

mad  = np.median( np.abs(pulls) )
nmad = 1.4826*mad
print "color: %s +- %s nmad %s" %(np.mean(pulls), np.std(pulls), nmad)
# print 'Fitted A Color = ', coeff[0]
# print 'Fitted mean Color = ', coeff[1]
# print 'Fitted standard deviation Color = ', coeff[2]


### X1r ###
plt.clf()
plt.errorbar(joined['x1r'][ind], joined['Skew'][ind], yerr=joined['Skew_coverr'][ind], fmt='o', alpha=0.6)
plt.plot([-3,3], [-3,3], color='k', lw=2.5)
plt.xlabel('input x1r')
plt.ylabel('measured x1r (SkewFit)')
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.savefig('./plots/x1r_v_x1r_skewfit.pdf')
# plt.show()

plt.clf()
plt.hist(joined['Skew'][ind] - joined['x1r'][ind], weights=joined['Skew_coverr'][ind], bins=25)
plt.xlabel('X1r(skewfit) - X1r(input)')
plt.ylabel('weighted histogram')
plt.xlim(-2,2)
plt.savefig('./plots/x1r_v_x1r_skewfit_hist.pdf')
# plt.show()

# x1r pulls
plt.clf()
# pulls
pulls = (joined['Skew'][ind] - joined['x1r'][ind])/joined['Skew_coverr'][ind]
#pull hist
hist, bin_edges, dummy = aml_hist(pulls, bins='scott', alpha=0.25)
#gaussian fit and plot
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
x = np.linspace(np.min(pulls), np.max(pulls), 150)
y = gauss(x, *coeff)
plt.plot(x,y, label='Fitted data')
#other plotting stuff
plt.annotate('mu: %s'%np.round(coeff[1],4), xy=(0.6,0.8), xycoords='axes fraction')
plt.annotate('sigma: %s'%np.round(coeff[2],4), xy=(0.6,0.7), xycoords='axes fraction')
plt.xlabel('X1r pulls')
plt.ylabel('histogram')
plt.xlim(-5,5)
plt.savefig('./plots/x1r_v_x1r_skewfit_pulls.pdf')
# plt.show()

mad  = np.median( np.abs(pulls) )
nmad = 1.4826*mad
print "X1r: %s +- %s nmad %s" %(np.mean(pulls), np.std(pulls), nmad)
# print 'Fitted A x1r = ', coeff[0]
# print 'Fitted mean x1r = ', coeff[1]
# print 'Fitted standard deviation x1r = ', coeff[2]

### X1f ###
plt.clf()
plt.errorbar(joined['x1f'][ind], joined['X1'][ind], yerr=joined['X1_coverr'][ind], fmt='o', alpha=0.6)
plt.plot([-3,3], [-3,3], color='k', lw=2.5)
plt.xlabel('input x1f')
plt.ylabel('measured x1f (SkewFit)')
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.savefig('./plots/x1f_v_x1f_skewfit.pdf')
# plt.show()

pulls = (joined['X1'][ind] - joined['x1f'][ind])/joined['X1_coverr'][ind]
plt.clf()
plt.plot(pulls, joined['X1_coverr'][ind],'o', alpha=0.6)
# plt.plot([-3,3], [-3,3], color='k', lw=2.5)
plt.xlabel('pulls')
plt.ylabel('x1f_err')
# plt.xlim(-3,3)
# plt.ylim(-3,3)
plt.savefig('./plots/x1fpull_v_x1ferr_skewfit.pdf')
# plt.show()

plt.clf()
plt.hist(joined['X1'][ind] - joined['x1f'][ind], weights=joined['X1_coverr'][ind], bins=25)
plt.xlabel('X1f(skewfit) - X1f(input)')
plt.ylabel('weighted histogram')
plt.xlim(-2,2)
plt.savefig('./plots/x1f_v_x1f_skewfit_hist.pdf')
# plt.show()

# x1r pulls
plt.clf()
pulls = (joined['X1'][ind] - joined['x1f'][ind])/joined['X1_coverr'][ind]
# pulls
#pull hist
hist, bin_edges, dummy = aml_hist(pulls, bins='scott', alpha=0.25)
#gaussian fit and plot
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
x = np.linspace(np.min(pulls), np.max(pulls), 150)
y = gauss(x, *coeff)
plt.plot(x,y, label='Fitted data')
#other plotting stuff
plt.annotate('mu: %s'%np.round(coeff[1],4), xy=(0.6,0.8), xycoords='axes fraction')
plt.annotate('sigma: %s'%np.round(coeff[2],4), xy=(0.6,0.7), xycoords='axes fraction')
plt.xlabel('X1f pulls')
plt.ylabel('histogram')
plt.xlim(-5,5)
plt.savefig('./plots/x1f_v_x1f_skewfit_pulls.pdf')

mad  = np.median( np.abs(pulls) )
nmad = 1.4826*mad
print "x1f: %s +- %s nmad %s" %(np.mean(pulls), np.std(pulls), nmad)

# plt.clf()
# plt.errorbar(joined['tmax'][ind], joined['DayMax'][ind], yerr=joined['DayMax_coverr'][ind], fmt='o')
# plt.xlabel('input DayMax')
# plt.ylabel('measured DayMax')
# plt.show()

### DayMax ###
plt.clf()
plt.hist(joined['DayMax'][ind] - joined['tmax'][ind], weights=joined['DayMax_coverr'][ind], bins=25)
plt.xlabel('DayMax(skewfit) - DayMax(input)')
plt.ylabel('weighted histogram')
plt.savefig('./plots/DayMax_v_DayMax_skewfit_hist.pdf')
# plt.show()

plt.clf()
pulls = (joined['DayMax'][ind] - joined['tmax'][ind])/joined['DayMax_coverr'][ind]
# pulls
#pull hist
hist, bin_edges, dummy = aml_hist(pulls, bins='scott', alpha=0.25)
#gaussian fit and plot
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
x = np.linspace(np.min(pulls), np.max(pulls), 150)
y = gauss(x, *coeff)

embed()
sys.exit()
plt.plot(x,y, label='Fitted data')
#other plotting stuff
plt.annotate('mu: %s'%np.round(coeff[1],4), xy=(0.6,0.8), xycoords='axes fraction')
plt.annotate('sigma: %s'%np.round(coeff[2],4), xy=(0.6,0.7), xycoords='axes fraction')
plt.xlabel('DayMax pulls')
plt.ylabel('histogram')
plt.xlim(-5,5)
plt.savefig('./plots/DayMax_v_DayMax_skewfit_pulls.pdf')

mad  = np.median(np.abs(pulls))
nmad = 1.4826*mad
print "DayMax: %s +- %s nmad %s" %(np.mean(pulls), np.std(pulls), nmad)

# plt.clf()
# plt.errorbar(joined['mB'][ind], joined['RestFrameMag_0_B'][ind], yerr=joined['RestFrameMag_0_B_coverr'][ind], fmt='o', alpha=0.6)
# plt.plot([13,26], [13,26], color='k', lw=2)
# plt.xlabel('input mB')
# plt.ylabel('measured mB (SkewFit)')
# plt.xlim(13,26)
# plt.ylim(13,26)
# plt.savefig('./plots/mB_v_mB_skewfit.pdf')
# plt.show()
#
# plt.clf()
# aml_hist(joined['RestFrameMag_0_B'][ind] - joined['mB'][ind], weights=joined['RestFrameMag_0_B_coverr'][ind], bins='scott')
# plt.xlabel('mB(skewfit) - mB(input)')
# plt.ylabel('weighted histogram')
# plt.savefig('./plots/mB_v_mB_skewfit_hist.pdf')
# plt.show()

x1r_m_x1f_skewfit = joined['Skew'][ind] - joined['X1'][ind]
x1r_m_x1f_input = joined['x1r'][ind] - joined['x1f'][ind]
c_m_c_skewfit = joined['Color'][ind] - joined['color'][ind]
x1diff_err = np.sqrt( (joined['Skew_coverr'][ind])**2 + (joined['X1_coverr'][ind])**2 )
cdiff_err = joined['Color_coverr'][ind]

plt.clf()
plt.errorbar(x1r_m_x1f_skewfit - x1r_m_x1f_input, c_m_c_skewfit, xerr=x1diff_err, yerr=cdiff_err, fmt='o', alpha=0.4)
plt.xlabel('x1r-x1f(skewfit) - x1r-x1f(input)')
plt.ylabel('Color(skewfit) - color(input)')
plt.savefig('./plots/skew_v_color_correlations.pdf')
# plt.show()

plt.clf()
plt.errorbar(joined['mwebv'][ind], c_m_c_skewfit, yerr=cdiff_err, fmt='o', alpha=0.6)
plt.xlabel('MW E(B-V)')
plt.ylabel('Color(skewfit) - color(input)')
plt.savefig('./plots/colordiff_v_mwebv_skewfit.pdf')

# plt.show()

# ind2 = np.where(c_m_c_skewfit < -0.1)
# tmp = joined[ind][ind2]
# embed()


# embed()