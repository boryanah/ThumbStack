import gc

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

# define
pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits' # DR6 ILC
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended.fits' # Already smoothed, I am 90% sure

"""
# plot map
ACT_map_sm = enmap.read_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits")
ACT_map_lg = enmap.read_map(pathMap.split(".fits")[0]+"_large_tau_screening.fits")
eshow(ACT_map_sm, pathMap.split(".fits")[0]+"_small_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
eshow(ACT_map_lg, pathMap.split(".fits")[0]+"_large_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
quit()
"""

# read
ACT_map = enmap.read_map(pathMap)
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
ell_ksz = np.arange(10001)
flg_ksz = np.zeros(len(ell_ksz))
flg_ksz[ell_ksz < 2000.] = 1.
choice = (ell_ksz >= 2000.) & (ell_ksz < 2150.)
flg_ksz[choice] = np.cos((ell_ksz[choice] - 2000.)*np.pi/300.)
fsm_ksz = np.zeros(len(ell_ksz))
fsm_ksz[ell_ksz > 2500.] = 1.
choice = (ell_ksz <= 2500.) & (ell_ksz > 2350.)
fsm_ksz[choice] = np.sin((ell_ksz[choice] - 2350.)*np.pi/300.)

# We can now filter the kappa alms by the filter produced above
ACT_alms_sm = hp.almxfl(ACT_alms, fsm_ksz)
ACT_alms_lg = hp.almxfl(ACT_alms, flg_ksz)

# We now convert the filtered alms to map space
ACT_map_sm = cs.alm2map(ACT_alms_sm, enmap.empty(shape, wcs))
ACT_map_lg = cs.alm2map(ACT_alms_lg, enmap.empty(shape, wcs))

# write map
enmap.write_map(pathMap.split(".fits")[0]+"_small_tau_screening.fits", ACT_map_sm)
enmap.write_map(pathMap.split(".fits")[0]+"_large_tau_screening.fits", ACT_map_lg)
