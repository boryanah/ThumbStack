import gc, sys

from pixell import enmap, enplot, reproject, curvedsky as cs, utils
import numpy as np

import pandas as pd
from astropy.io import fits

import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def eshow(x, fn, **kwargs):
    ''' Define a function to help us plot the maps neatly '''
    plots = enplot.get_plots(x, **kwargs)
    #enplot.show(plots, method = "python")
    enplot.write(fn, plots)

typ = sys.argv[1]

"""
cd /global/homes/b/boryanah/repos/ThumbStack/desi/screening
python filter_map.py dr5_f090
python filter_map.py dr5_f150
python filter_map.py dr6
"""

# define
if typ == "dr5_f150":
    #pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_rebeamed_f090.fits' # DR5 f150
    pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_masked.fits' # DR5 f150
elif typ == "dr5_f090":
    #pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree.fits' # DR5 f090
    pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked_rebeamed_f150.fits' # DR5 f090
else:
    pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits' # DR6 ILC
    #pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/ilc_fullRes_TT.fits' # DR6 ILC
#pathMap = '/pscratch/sd/b/boryanah/websky/unlensed_map.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/lensed_map.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_4096_ph201.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_fwhm1.6_noise15.0.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_tau_8192_ph201.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tau_8192_ph201.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_fwhm1.6.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_fwhm2.1.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_fwhm1.4.fits' #
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_8192_fwhm1.4.fits' #
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_order1_fwhm1.4.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_8192_ph201_order1_brute_fwhm1.4.fits' 
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_fwhm1.6.fits' 
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_tau_8192_ph201_fwhm1.6.fits' 
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tau_8192_ph201_fwhm1.6.fits' 
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tau_8192_ph201_fwhm2.1.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tau_8192_ph201_fwhm1.4.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tau_8192_ph201_MTNG_fwhm1.4.fits'
#pathMap = '/pscratch/sd/b/boryanah/websky/abacus/lensed_map_tau_8192_ph201_MTNG_fwhm1.4.fits' 
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended.fits' # Already smoothed, I am 90% sure

"""
#ACT_map_sm = enmap.read_map(pathMap.split(".fits")[0]+"_half_small_large_small_tau_screening.fits")
ACT_map_sm = enmap.read_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits")
#eshow(ACT_map_sm, pathMap.split(".fits")[0]+"_half_small_large_small_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 20})
eshow(ACT_map_sm, pathMap.split(".fits")[0]+"_small_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 20})
quit()
"""

want_multiply = True
if want_multiply:
    ACT_map_sm = enmap.read_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits")
    ACT_map_lg = enmap.read_map(pathMap.split(".fits")[0]+"_large_tau_screening.fits")
    #enmap.write_map(pathMap.split(".fits")[0]+"_small_large_tau_screening.fits", ACT_map_sm * np.sign(ACT_map_lg))
    #quit()
    
    want_plot = False
    if want_plot:
        # plot map
        ACT_map_sm = enmap.read_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits")
        ACT_map_lg = enmap.read_map(pathMap.split(".fits")[0]+"_large_tau_screening.fits")
        eshow(ACT_map_sm, pathMap.split(".fits")[0]+"_small_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
        eshow(ACT_map_lg, pathMap.split(".fits")[0]+"_large_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
        quit()

    ACT_map = ACT_map_sm*ACT_map_lg
else:
    # read
    ACT_map = enmap.read_map(pathMap)

    
if ACT_map.shape[0] == 3:
    ACT_map = ACT_map[0]
ACT_mask = enmap.read_map(pathMask)

# mask
ACT_map *= ACT_mask 
shape, wcs = ACT_map.shape, ACT_map.wcs
del ACT_mask
gc.collect()

# convert to alms (curved sky)
ACT_alms = cs.map2alm(ACT_map, lmax=10000)
del ACT_map
gc.collect()

# load ells for filter
#L_cut = 2000. # og
#L_jump = 500. # og
#L_width = 150. # og
L_cut = 1700.
L_jump = 300.
L_width = 100.
factor = 1. # 1 is the OG; 2 is the new one (when multiply is True, should be 1; first run should have multiply False and factor = 2)
L_cut /= factor
L_jump /= factor
L_width /= factor
ell_ksz = np.arange(10001)
flg_ksz = np.zeros(len(ell_ksz))
flg_ksz[ell_ksz < L_cut] = 1.
choice = (ell_ksz >= L_cut) & (ell_ksz < L_cut+L_width)
flg_ksz[choice] = np.cos((ell_ksz[choice] - L_cut)*np.pi/(2.*L_width))
fsm_ksz = np.zeros(len(ell_ksz))
fsm_ksz[ell_ksz > L_cut+L_jump] = 1.
choice = (ell_ksz <= L_cut+L_jump) & (ell_ksz > L_cut+L_jump-L_width)
fsm_ksz[choice] = np.sin((ell_ksz[choice] - L_cut+L_jump-L_width)*np.pi/(2.*L_width))

# We can now filter the kappa alms by the filter produced above
ACT_alms_sm = hp.almxfl(ACT_alms, fsm_ksz)
if not want_multiply:
    ACT_alms_lg = hp.almxfl(ACT_alms, flg_ksz)

# We now convert the filtered alms to map space
ACT_map_sm = cs.alm2map(ACT_alms_sm, enmap.empty(shape, wcs))
if not want_multiply:
    ACT_map_lg = cs.alm2map(ACT_alms_lg, enmap.empty(shape, wcs))

# write map
if want_multiply:
    enmap.write_map(pathMap.split(".fits")[0]+"_half_small_large_small_tau_screening.fits", ACT_map_sm)
else:
    enmap.write_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits", ACT_map_sm)
    enmap.write_map(pathMap.split(".fits")[0]+"_large_tau_screening.fits", ACT_map_lg)
