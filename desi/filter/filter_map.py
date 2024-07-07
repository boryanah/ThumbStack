import gc

from pixell import enmap, utils, powspec, enplot, reproject #, pointsrcs
import numpy as np

import pandas as pd
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

# read
ACT_map = enmap.read_map(pathMap)
ACT_mask = enmap.read_map(pathMask)

# mask
ACT_map *= ACT_mask
shape, wcs = ACT_map.shape, ACT_map.wcs
del ACT_mask
gc.collect()

# fourier transform
kmap = enmap.fft(ACT_map, normalize="phys") 
#lymap, lxmap = enmap.lmap(shape, wcs) # matches kmap
modlmap = enmap.modlmap(shape, wcs)
print(kmap.shape)
print(kmap.dtype)

# load kSZ filter 
fl_ksz = np.load("data/ACT_filter_taper_kSZ.npy") # F(ell)
ell_ksz = np.load("data/ACT_ell_kSZ.npy")
fl_ksz /= np.max(fl_ksz)
fl_fun = interp1d(ell_ksz, fl_ksz, bounds_error=False, fill_value=0.)
print("coverage", ell_ksz.min(), ell_ksz.max(), modlmap.min(), modlmap.max())

# apply
fl_map = fl_fun(modlmap)
kmap *= fl_map
ACT_map_tmp = enmap.ifft(kmap, normalize="phys")
ACT_map[:, :] = np.real(ACT_map_tmp[:, :])
print(ACT_map_tmp.shape)
print(ACT_map.dtype)
del kmap, ACT_map_tmp
gc.collect()

# write map
enmap.write_map(pathMap.split(".fits")[0]+"_kSZ_filter.fits", ACT_map)

# plot map
eshow(ACT_map, pathMap.split(".fits")[0]+"_kSZ_filter", **{"colorbar":True, "ticks": 5, "downgrade": 4})

# use in usual measurement (but also produce stacked image)
