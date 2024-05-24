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

# params
desi_dir = Path("/pscratch/sd/b/boryanah/ACTxDESI/DESI/")
#tracer = "BGS_BRIGHT-21.5" # og
tracer = "BGS_BRIGHT" # new
pz_bin = 1
save_fn = f"{tracer}.txt"
#cat_NGC = Table.read(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_NGC.fits", format='fits') # og
#cat_SGC = Table.read(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_SGC.fits", format='fits') # og
cat_NGC = Table.read(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_NGC_logM11.0.fits", format='fits') # new
cat_SGC = Table.read(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_SGC_logM11.0.fits", format='fits') # new

cat = vstack((cat_NGC, cat_SGC))

# TESTING trying to see if it helps with first 2 bins -- not really
choice = cat['Z'] > 0.08
    cat = cat[choice]

# a bit weird
cat.remove_column('BITWEIGHTS')

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
df2['MStellar'] = 1.e11
massConversion = MassConversionKravtsov14()

# vTheta and vPhi are incorrect, but it doesn't matter as long as vR is fine
df2['vR'] = cat['VEL_LOS']
df2['vTheta'] = cat['VEL_ACROSS']
df2['vPhi'] = cat['VEL_ACROSS']

# Saving the fits catalogs as txt files
df2.to_csv(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}', header=None, index=None, sep=' ', mode='w')
print(f'Saved DESI_pz{pz_bin:d}')

# TODO: do we need to update this?
# cosmological parameters
u = UnivMariana()

# Lastly we read & save these again using the built in ThumbStack catalog function, so it adds the last few columns (like mVir) for us.
Catalog(u, massConversion, name=f"DESI_pz{pz_bin:d}", nameLong=f"DESI pz bin {pz_bin:d}", pathInCatalog=str(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}'), save=True, cat_fn=save_fn)
