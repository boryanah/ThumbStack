import sys
sys.path.append('/global/homes/b/boryanah/repos/2MPZ_vel//pairwise')

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

from pixell import powspec, utils
from utils_power import bin_mat

# location
data_dir = "/global/homes/b/boryanah/repos/2MPZ_vel/pairwise/camb_data/"

# load cmb power from theory
camb_theory = powspec.read_spectrum(f"{data_dir}/camb_theory.dat", scale=True) # scaled by 2pi/l/(l+1) to get C_ell
cl_th = camb_theory[0, 0, :10000]
ell_th = np.arange(cl_th.size)
ell_ksz, cl_ksz = np.loadtxt(f"{data_dir}/cl_ksz_bat.dat", unpack=True)
cl_ksz /= (ell_ksz*(ell_ksz+1)/(2.*np.pi))

#plt.plot(ell_ksz, cl_ksz)
#plt.xscale('log')
#plt.show()

# load power
#survey = "Planck"
survey = "ACT"
if survey == "Planck":
    cl_data = np.load(f"{data_dir}/Planck_power.npy")
    ell_data = np.load(f"{data_dir}/Planck_ell.npy")
elif survey == "ACT":
    #cl_data = np.load("../pairwise/camb_data/ACT_DR5_f090_65uK_binned_power.npy")
    #ell_data = np.load("../pairwise/camb_data/ACT_DR5_f090_65uK_centers.npy")
    cl_data = np.load(f"{data_dir}/ACT_D56_binned_power.npy")
    ell_data = np.load(f"{data_dir}/ACT_D56_centers.npy")

# bin the power spectrum
bins = np.linspace(100, 4000, 300)
ell_data_binned, cl_data_binned = bin_mat(ell_data, cl_data, bins)

# predict the noise power spectrum
if survey == "Planck":
    Delta_T = 47. # uK arcmin
    theta_FWHM = 5. # arcmin
elif survey == "ACT":
    Delta_T = 10. # uK arcmin
    theta_FWHM = 1.6 # arcmin
Delta_T *= utils.arcmin
theta_FWHM *= utils.arcmin
    
# noise power is Delta_T^2 b(ell)^-2
nl_th = Delta_T**2.*np.exp(theta_FWHM**2*ell_th**2/(8.*np.log(2.)))
bl_th = np.exp(-0.5*theta_FWHM**2*ell_th**2/(8.*np.log(2.)))
bl_data_binned = np.exp(-0.5*theta_FWHM**2*ell_data_binned**2/(8.*np.log(2.)))

# get the filtering function
if survey == "Planck":
    ell_min, ell_max = 100, 3000
elif survey == "ACT":
    ell_min, ell_max = 100, 10000
fl_data_binned = np.interp(ell_data_binned, ell_ksz, cl_ksz)/cl_data_binned
fl_data_binned[(ell_data_binned < ell_min) | (ell_data_binned > ell_max)] = 0.
fl_data = np.interp(ell_data, ell_ksz, cl_ksz)/cl_data
fl_data[(ell_data < ell_min) | (ell_data > ell_max)] = 0.    
fl_th = np.interp(ell_th, ell_ksz, cl_ksz)/((cl_th+nl_th)*bl_th**2)
fl_th[(ell_th < ell_min) | (ell_th > ell_max)] = 0.
fl_data_th = np.interp(ell_th, ell_data_binned, fl_data_binned)

# wait so why does og save data and not theory?
np.save(f"data/{survey}_filter_kSZ.npy", fl_th)
#np.save(f"data/{survey}_filter_kSZ.npy", fl_data_th) # og
np.save(f"data/{survey}_ell_kSZ.npy", ell_th)

# taper
a = 0.1
x = 4.60512/a # solution of 0.5 (1+/-tanh(ax)) = 0.0001 is x = -/+4.60512/a 
ell_min = ell_min + x
ell_max = ell_max - x
taper_th = 0.5*(1. + np.tanh(a*(ell_th-ell_min))) + 0.5*(1. - np.tanh(a*(ell_th-ell_max))) - 1.
fl_data_th_taper = fl_data_th*taper_th
fl_th_taper = fl_th*taper_th
#np.save(f"data/{survey}_filter_taper_kSZ.npy", fl_data_th_taper) # og
np.save(f"data/{survey}_filter_taper_kSZ.npy", fl_th_taper)

plt.plot(ell_th, taper_th, label='Taper')
plt.savefig(f"figs/{survey}_taper.png")
plt.show()

# for plotting
power = 2.

plt.figure(1)
plt.title("CMB power")
plt.plot(ell_data, cl_data*ell_data**power)
plt.plot(ell_data_binned, cl_data_binned*ell_data_binned**power, label="binned")
plt.plot(ell_th, cl_th*ell_th**power, lw=1, color='k')
plt.plot(ell_th, (cl_th+nl_th)*bl_th**2*ell_th**power, lw=2, color='k')
plt.plot(ell_th, nl_th*ell_th**power, lw=1, color='k')
plt.legend()
plt.xlim([0, 3000])
plt.ylim([0, 30000])
plt.savefig(f"figs/{survey}_power.png")
#plt.show()

plt.figure(2)
plt.title("filter = kSZ/Cl_data")
plt.plot(ell_data_binned, fl_data_binned*bl_data_binned/np.max(fl_data_binned*bl_data_binned), label=r'$f(\ell) = C_\ell^{kSZ} b_\ell/(C_\ell^{TT,data})$')
#plt.plot(ell_data, fl_data/np.max(fl_data))
plt.plot(ell_th, fl_th*bl_th/np.max(fl_th*bl_th), label=r'$f(\ell) = C_\ell^{kSZ} b_\ell/[(C_\ell^{TT,theo}+N_\ell^{TT})/b_\ell^2]$')
plt.plot(ell_th, fl_data_th_taper*bl_th/np.max(fl_data_th_taper*bl_th), label='Tapered')
plt.legend()
plt.xlim([0, 10000])
plt.savefig(f"figs/{survey}_filter.png")

plt.figure(3)
plt.title("kSZ file David sent (normed)")
plt.plot(ell_ksz, cl_ksz*ell_ksz**2)
plt.xlim([0, 10000])
plt.savefig(f"figs/{survey}_ksz.png")
plt.show()
