import gc
import os
from pathlib import Path
import sys
from importlib import reload
sys.path.append('/global/u1/b/boryanah/repos/ThumbStack/')
sys.path.append('/global/u1/b/boryanah/reconstruct_DESI/')

import numpy as np
import pandas as pd
import fitsio
import healpy as hp
import astropy
from astropy.table import Table, vstack
from astropy.io import fits

# ThumbStack
import catalog
reload(catalog)
from catalog import *

import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

# Reconstruct DESI
from visualize_recon import calc_velocity
from pixell import enmap, enplot, reproject, curvedsky as cs, utils

# params
desi_dir = Path("/pscratch/sd/b/boryanah/ACTxDESI/DESI/")
tracer = "BGS_BRIGHT_pz"
#logm_str = "" # og
logm_str = "_logm10.5" # TESTING
pz_bin = 1
save_fn = f"{tracer}{logm_str}.txt"

# filenames
cat_fn = f"/global/cfs/cdirs/desi/users/boryanah/reconstruction_DESI/galaxies/dr9_bgs_basic_remov_isle_nobs1_ebv0.15_stardens2500_sigmaz0.0500{logm_str}.npz"
recon_fn = f"/global/cfs/cdirs/desi/users/boryanah/reconstruction_DESI/recon/displacements_dr9_bgs_basic_randoms-1-0-4_remov_isle_nobs1_ebv0.15_stardens2500_sigmaz0.0500{logm_str}_R12.50_nmesh1024_recsym_MG.npz"

# read CMB map
pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits'
cmb_map = enmap.read_map(pathMap.split('.fits')[0]+"_large_tau_screening.fits")

# load catalog file
data = np.load(cat_fn)
cat = {}
cat['RA'] = data['RA']
cat['DEC'] = data['DEC']
cat['Z'] = data['Z_PHOT_MEDIAN']
LOGM = data['LOGM']
cat = Table(cat)
del data; gc.collect()

# Converting to Pandas and rename the columns
print(cat.keys())
df = cat.to_pandas()
df2 = df.loc[:, ('RA', 'DEC', 'Z')]

# New Columns with Mstellar given from above
colnames = ['coordX', 'coordY', 'coordZ', 'dX', 'dY',
            'dZ', 'dXKaiser', 'dYKaiser', 'dZKaiser',
            'vX', 'vY', 'vZ', 'vR', 'vTheta', 'vPhi']

# set to 0
for col in colnames:
    df2[col] = 0

# set the stellar mass
df2['MStellar'] = 10.**LOGM #1.e11 
del LOGM; gc.collect()
massConversion = MassConversionKravtsov14()

# load reconstruction file
data = np.load(recon_fn)
# recon scalars
unit_los = data['unit_los']
displacements = data['displacements']
assert unit_los.shape[0] == len(cat)
mean_z = data['mean_z']
growth_factor = data['growth_factor']
Hubble_z = data['Hubble_z']
h = data['h']
del data; gc.collect()
    
# calculate reconstructed velocities
Velocity, Velocity_sphere = calc_velocity(displacements, unit_los, 1./(1.+mean_z), growth_factor, Hubble_z, h, want_rsd=False)
df2['vX'] = Velocity[:, 0]
df2['vY'] = Velocity[:, 1]
df2['vZ'] = Velocity[:, 2]

# vTheta and vPhi are incorrect, but it doesn't matter as long as vR is fine
df2['vR'] = Velocity_sphere[:, 2]
df2['vTheta'] = Velocity_sphere[:, 0]
df2['vPhi'] = Velocity_sphere[:, 1]

# interpolate the map to the given sky coordinates
sourcecoord = np.array([df2['DEC'], df2['RA']]) * np.pi/180.   # convert from degrees to radians

# use nearest neighbor interpolation # tuks cmb_map
df2['vZ'] = cmb_map.at(sourcecoord, prefilter=False, mask_nan=False, order=0) # note we are overwriting
print("number of galaxies in bin", sourcecoord.shape[1])
del sourcecoord; gc.collect()

# Saving the fits catalogs as txt files
df2.to_csv(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}', header=None, index=None, sep=' ', mode='w')
print(f'Saved DESI_pz{pz_bin:d}')

# TODO: do we need to update this?
# cosmological parameters
u = UnivMariana()

# Lastly we read & save these again using the built in ThumbStack catalog function, so it adds the last few columns (like mVir) for us.
Catalog(u, massConversion, name=f"DESI_pz{pz_bin:d}", nameLong=f"DESI pz bin {pz_bin:d}", pathInCatalog=str(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}'), save=True, cat_fn=save_fn)
