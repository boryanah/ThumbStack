import time

import numpy as np
from astropy.table import Table
import pandas as pd

from cosmoprimo.fiducial import Planck2018FullFlatLCDM, AbacusSummit, DESI, TabulatedDESI

# Tracer type
#tracer = "LRG" # og
tracer = "BGS_BRIGHT-21.5" #


# Importing catalogs
main_directory = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1/unblinded/'
post_rec_directory = "desipipe/baseline_2pt/recon_recsym/"

# North Galactic Center:
pre_rec_NGC = Table.read(main_directory+f'{tracer}_NGC_clustering.dat.fits', format='fits')
#pre_rec_NGC = pre_rec_NGC.to_pandas()

post_rec_NGC = Table.read(main_directory+post_rec_directory+f'{tracer}_NGC_clustering.dat.fits', format='fits')
#post_rec_NGC = post_rec_NGC.to_pandas()

# South Galactic Center:
pre_rec_SGC = Table.read(main_directory+f'{tracer}_SGC_clustering.dat.fits', format='fits')
#pre_rec_SGC = pre_rec_SGC.to_pandas()

post_rec_SGC = Table.read(main_directory+post_rec_directory+f'{tracer}_SGC_clustering.dat.fits', format='fits')
#post_rec_SGC = post_rec_SGC.to_pandas()

# Numbers are low!
print("post pre NGC", len(post_rec_NGC), len(pre_rec_NGC))
print("post pre NGC", len(post_rec_SGC), len(pre_rec_SGC))

# All (for visualization purposes)
#pre_rec = pd.concat([pre_rec_NGC, pre_rec_SGC])

"""
# For matching to rongpu's bins
if tracer == "LRG":
    bins_directory = "/pscratch/sd/b/bried24/ACTxDESI/DESI/dr9_lrg_pzbins.fits"
if "BGS" in tracer:
    bins_directory = "/pscratch/sd/b/boryanah/kSZ_pairwise/BGS_BRIGHT_full_noveto_vac_marvin.dat.fits"
bins = Table.read(bins_directory, format='fits')
bins = bins.to_pandas()

# select only high mass objects
choice = (bins['Z'] < 0.7) & (bins['Z'] > 0.08) & (bins['LOGM'] > 11)

# choose only objects above a certain mass threshold

# Start timer
start_time = time.perf_counter()
pre_rec_NGC = pd.merge(pre_rec_NGC, bins, on=["TARGETID", "RA", "DEC", "Z"])
pre_rec_SGC = pd.merge(pre_rec_SGC, bins, on=["TARGETID", "RA", "DEC", "Z"])

# this is the new thing!!!!!!!
post_rec_NGC = pd.merge(post_rec_NGC, bins, on=["TARGETID", "RA", "DEC", "Z"])
post_rec_SGC = pd.merge(post_rec_SGC, bins, on=["TARGETID", "RA", "DEC", "Z"])

# End timer
print("pre_rec_NGC.columns", pre_rec_NGC.columns)
end_time = time.perf_counter()

# Calculate elapsed time
elapsed_time = end_time - start_time
print("Elapsed time: ", elapsed_time)
"""

# Let's now obtain the displacements and velocities
cosmo = DESI()

# Displacements:
chi_pre_NGC = cosmo.comoving_radial_distance(pre_rec_NGC['Z'])
chi_post_NGC = cosmo.comoving_radial_distance(post_rec_NGC['Z'])

chi_pre_SGC = cosmo.comoving_radial_distance(pre_rec_SGC['Z'])
chi_post_SGC = cosmo.comoving_radial_distance(post_rec_SGC['Z'])

disp_NGC_LOS = chi_pre_NGC - chi_post_NGC
disp_SGC_LOS = chi_pre_SGC - chi_post_SGC

pre_rec_NGC["DISP_LOS"] = disp_NGC_LOS
pre_rec_SGC["DISP_LOS"] = disp_SGC_LOS

pre_rec_NGC["chi_pre"] = chi_pre_NGC
pre_rec_SGC["chi_pre"] = chi_pre_SGC


# Angular difference:
delta_theta_NGC = np.arccos(np.sin(np.radians(pre_rec_NGC['DEC']))*np.sin(np.radians(post_rec_NGC['DEC']))
                            +np.cos(np.radians(pre_rec_NGC['DEC']))*np.cos(np.radians(post_rec_NGC['DEC']))
                            *np.cos(np.radians(pre_rec_NGC['RA'])-np.radians(post_rec_NGC['RA'])))

delta_theta_SGC = np.arccos(np.sin(np.radians(pre_rec_SGC['DEC']))*np.sin(np.radians(post_rec_SGC['DEC']))
                            +np.cos(np.radians(pre_rec_SGC['DEC']))*np.cos(np.radians(post_rec_SGC['DEC']))
                            *np.cos(np.radians(pre_rec_SGC['RA'])-np.radians(post_rec_SGC['RA'])))

disp_NGC_ACROSS = chi_pre_NGC*delta_theta_NGC
disp_SGC_ACROSS = chi_pre_SGC*delta_theta_SGC

# Velocities:
ff_NGC = cosmo.growth_factor(pre_rec_NGC['Z'])
H_z_NGC = cosmo.hubble_function(pre_rec_NGC['Z'])
a_NGC = 1/(1+pre_rec_NGC['Z'])

ff_SGC = cosmo.growth_factor(pre_rec_SGC['Z'])
H_z_SGC = cosmo.hubble_function(pre_rec_SGC['Z'])
a_SGC = 1/(1+pre_rec_SGC['Z'])

vel_NGC_LOS = a_NGC * H_z_NGC * ff_NGC * disp_NGC_LOS
pre_rec_NGC["VEL_LOS"] = vel_NGC_LOS

vel_SGC_LOS = a_SGC * H_z_SGC * ff_SGC * disp_SGC_LOS
pre_rec_SGC["VEL_LOS"] = vel_SGC_LOS


vel_NGC_ACROSS = a_NGC * H_z_NGC * ff_NGC * disp_NGC_ACROSS
pre_rec_NGC["VEL_ACROSS"] = vel_NGC_ACROSS

vel_SGC_ACROSS = a_SGC * H_z_SGC * ff_SGC * disp_SGC_ACROSS
pre_rec_SGC["VEL_ACROSS"] = vel_SGC_ACROSS


pre_rec_NGC.write(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_NGC.fits", overwrite=False)
pre_rec_SGC.write(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_SGC.fits", overwrite=False)

#pre_rec_NGC.write(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_NGC_z0.08-0.7_logM11.dat.fits", overwrite=False)
#pre_rec_SGC.write(f"/pscratch/sd/b/boryanah/kSZ_pairwise/{tracer}_SGC_z0.08-0.7_logM11.dat.fits", overwrite=False)
