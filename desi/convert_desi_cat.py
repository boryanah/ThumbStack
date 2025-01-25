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
from astropy.table import Table
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


"""
Oburnali sme vZ samo na extended dr9 i extended dr9 sigmaz0.05

python convert_desi_cat.py 1 1
python convert_desi_cat.py 0 0
python convert_desi_cat.py 1 0 # best?
python convert_desi_cat.py 0 1

python convert_desi_cat.py 1 0 dr9_lrg_pzbins.fits
python convert_desi_cat.py 1 0 dr9_extended_lrg_pzbins.fits
python convert_desi_cat.py 1 0 dr9_lrg_pzbins-dr10_pz.fits
python convert_desi_cat.py 1 0 dr9_extended_lrg_pzbins-dr10_pz.fits

python convert_desi_cat.py 1 0 dr9_extended_lrg_pzbins.fits f150
python convert_desi_cat.py 1 0 dr9_extended_lrg_pzbins.fits f090

Boryana, what the actual fuck. This code makes absolutely no sense....!!!!!!!!!!!!!!!!!!
"""

# choose galaxy and randoms
gal_fn = sys.argv[3]
#gal_fn = "dr9_lrg_pzbins.fits" #
#gal_fn = "dr9_extended_lrg_pzbins.fits" #
#gal_fn = "dr9_lrg_pzbins-dr10_pz.fits" #
#gal_fn = "dr9_extended_lrg_pzbins-dr10_pz.fits" #
#gal_fn = "mix_dr9_lrg_pzbins.fits" # new
#gal_fn = "mix_dr9_extended_lrg_pzbins.fits" # new

# mask params
min_nobs = 1 #2 # Rongpu says 1 it's ok (DR9 we used 2)
max_ebv = 0.15
max_stardens = 2500
max_sigmaz = 0. #0.05 #0. #0.05 #0.035
remove_islands = True
apply_lrg_mask = True
skip_tidal = True # TESTING!!!!!!!!!

# read CMB map
pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits'; freq_str = ""
if len(sys.argv) > 4:
    print(sys.argv[4])
    freq_str = f"_{sys.argv[4]}"
    if sys.argv[4] == "f090":
        pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked.fits' # DR5 f090
    elif sys.argv[4] == "f150":
        pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_masked.fits' # DR5 f150
cmb_map = enmap.read_map(pathMap.split('.fits')[0]+"_large_tau_screening.fits")

# mask strings
isle_str = "_remov_isle" if remove_islands else ""
lrg_mask_str = "_lrg_mask" if apply_lrg_mask else ""
sigmaz_str = f"_sigmaz{max_sigmaz:.4f}" if not np.isclose(max_sigmaz, 0.) else ""
mask_str = f"{isle_str}_nobs{min_nobs:d}_ebv{max_ebv:.2f}_stardens{max_stardens:d}{lrg_mask_str}{sigmaz_str}"

# only record one footprint
only_footprint = None #"DES" # "S" #"N" None
only_str = f"_{only_footprint}" if only_footprint is not None else ""

# choose footprint
#comb_footprint = False #
comb_footprint =  bool(int(sys.argv[1])) #
cat_foot_str = "_allfoot" if comb_footprint else "_joinfoot"

# recon params
recon_all_pz_bin = bool(int(sys.argv[2])) #
recon_bin_str = "_allbin" if recon_all_pz_bin else '_perbin'
if "extended" in gal_fn:
    start_str = "extended_"
else:
    start_str = "main_"
if "dr10" in gal_fn:
    want_dr10 = True
else:
    want_dr10 = False
if "mix" in gal_fn:
    from tools.match_searchsorted import match
    want_mix = True
else:
    want_mix = False
if recon_all_pz_bin:
    rand_fn = start_str+"randoms-1-0-9.fits" #
else:
    rand_fn = start_str+"randoms-1-0-2.fits" # typically
    #rand_fn = start_str+"randoms-1-0-6.fits" # TESTING dr10 no difference
if want_dr10:
    version_str = "_dr10"
    rand_fn = "dr10_"+rand_fn
else:
    version_str = ""
mix_str = "_mix" if want_mix else ""
if comb_footprint:
    nmesh = 2048 
else:
    nmesh = 1024 
sr = 12.5
convention = "recsym"
rectype = "MG"

doRandomPositions = False
if doRandomPositions:
    random_str = "_randompos"
else:
    random_str = ""

# filename to save to
desi_dir = Path("/pscratch/sd/b/boryanah/ACTxDESI/DESI/")
if "extended" in gal_fn:
    save_fn = f"extended_catalog{mix_str}{version_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{freq_str}.txt"
else:
    save_fn = f"catalog{mix_str}{version_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{freq_str}.txt"
print(save_fn)

# DESI catalog location
if want_dr10:
    #fn = "/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/"+(gal_fn.split('mix_')[-1]) # og
    fn = "/pscratch/sd/b/boryanah/kSZ_pairwise/"+(gal_fn.split('mix_')[-1]) # TESTING
else:
    #fn = '/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/data_products/catalogs/'+(gal_fn.split('mix_')[-1]) # og
    fn = "/pscratch/sd/b/boryanah/kSZ_pairwise/"+(gal_fn.split('mix_')[-1]) # TESTING
cat = Table(fitsio.read(fn)) # directory change!
print(len(cat))
print(cat.keys())

# load DES footprint
nside = 256
des_footprint = hp.fitsfunc.read_map(f"/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/data_products/misc/hp_in_des_{nside}_ring.fits.gz")
des_hp_idx = np.arange(len(des_footprint))[des_footprint]
del des_footprint; gc.collect()

# TODO: Update
# Converting a flat virial mass to Stellar Mass, this is what we're using for every object for now.
massConversion = MassConversionKravtsov14()
MStellar = massConversion.fmVirTomStar(2e13)

# We run the code for the four DESI redshift bins in a for loop.
bins = [1, 2, 3, 4]

# load the cleaned catalog used in reconstruction
recon_dir = Path("/global/cfs/cdirs/desi/users/boryanah/reconstruction_DESI/")
clean_gal_fn = f"{gal_fn.split('.fits')[0]}{mask_str}.npz"
clean_data = np.load(recon_dir / "galaxies" / clean_gal_fn)
h = clean_data['h'] # TODO: move to displacement
unit_los = clean_data['unit_los']
pz_bin_clean = clean_data['pz_bin']
del clean_data; gc.collect()

# load full galaxy catalog
cat = Table.read(fn, format='fits')

# initiate mask
mask = np.ones(len(cat['DEC']), dtype=bool)

# Remove "islands" in the NGC
if remove_islands:
    mask &= ~((cat['DEC']<-10.5) & (cat['RA']>120) & (cat['RA']<260))
    print('Remove islands', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))


if want_dr10:
    # Remove region with large fraction of bad photo-z's
    bad_exposure_ra, bad_exposure_dec = 33.64, 13.38
    mask &= ~((cat['RA'] > bad_exposure_ra-1.5) & (cat['RA'] < bad_exposure_ra+1.5) & \
              (cat['DEC'] > bad_exposure_dec-1.5) & (cat['DEC'] < bad_exposure_dec+1.5))
    print('Bad exposure region', np.sum(mask), np.sum(mask)/len(mask))

    # Remove region near LMC where data and randoms do not match
    mask &= ~((cat['RA'] > 68.) & (cat['RA'] < 87.) & (cat['DEC'] < -62.))
    print('Bad LMC region', np.sum(mask), np.sum(mask)/len(mask))

    # Identify and remove islands
    nside = 256
    pix = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    pix_to_keep = np.load(f'/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/misc/not_islands_pix_{nside}.npy')
    mask &= np.in1d(pix, pix_to_keep)
    print('Remove lonely islands', np.sum(mask), np.sum(mask)/len(mask))
    del pix; gc.collect()
    
# NOBS cut
if min_nobs > 0:
    mask &= (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
    print('NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))

# additional cuts for DR10
if want_dr10:
    # DR10 NOBS cut
    mask &= (cat['DR10_PIXEL_NOBS_G'] >= min_nobs) & (cat['DR10_PIXEL_NOBS_R'] >= min_nobs) & (cat['DR10_PIXEL_NOBS_I'] >= min_nobs) & (cat['DR10_PIXEL_NOBS_Z'] >= min_nobs)
    print('DR10 NOBS', np.sum(mask), np.sum(mask)/len(mask))
    
# Apply LRG mask
if apply_lrg_mask:
    mask &= cat['lrg_mask']==0
    print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

# Martin's EBV cut
if max_ebv < 1.:
    mask &= cat['EBV']<max_ebv
    print('EBV', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

# Martin's STARDENS cut
if max_stardens > 0:
    stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_64_ring.npy')
    stardens_nside = 64
    mask_star = stardens>=max_stardens
    bad_hp_idx = np.arange(len(stardens))[mask_star]
    cat_hp_idx = hp.pixelfunc.ang2pix(stardens_nside, cat['RA'], cat['DEC'], lonlat=True, nest=False)
    mask &= ~np.in1d(cat_hp_idx, bad_hp_idx)
    print('STARDENS', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

if not np.isclose(max_sigmaz, 0.):
    if want_dr10:
        mask &= 0.5*(cat['DR10_Z_PHOT_U68_I'] - cat['DR10_Z_PHOT_L68_I']) < max_sigmaz*(1. + cat['DR10_Z_PHOT_MEDIAN_I'])
    else:
        # more
        more_dir = Path("/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/data_products/catalogs/more/")
        more = Table(fitsio.read(more_dir / ((gal_fn.split('mix_')[-1]).split('_pzbins')[0]+'_pz.fits'), columns=['Z_PHOT_STD']))
        mask &= more['Z_PHOT_STD'] < max_sigmaz*(1. + cat['Z_PHOT_MEDIAN'])
        print('SIGMAZ', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
    
# apply all cuts
cat = cat[mask]
print(len(cat['RA']), unit_los.shape[0])
assert unit_los.shape[0] == len(cat['RA'])

# match pz_bin and Z_PHOT_MEDIAN
if want_mix:
    load_dir = Path("/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/")
    if "extended" in gal_fn:
        cat_dr10 = Table(fitsio.read(load_dir / "dr9_extended_lrg_pzbins-dr10_pz.fits"))
    else:
        cat_dr10 = Table(fitsio.read(load_dir / "dr9_lrg_pzbins-dr10_pz.fits"))
    ptr = match(np.asarray(cat['TARGETID']).astype(np.int64), np.asarray(cat_dr10['TARGETID']).astype(np.int64))
    cat['Z_PHOT_MEDIAN'][ptr > -1] = cat_dr10['DR10_Z_PHOT_MEDIAN_I'][ptr[ptr > -1]]
    cat['pz_bin'][ptr > -1] = cat_dr10['DR10_pz_bin'][ptr[ptr > -1]]
    del cat_dr10
    gc.collect()
    
# additional cut because had issues with reconstruction when using beyond the 4 bins
if want_dr10:
    mask_within_bins = (cat['DR10_pz_bin'] <= 4) & (cat['DR10_pz_bin'] >= 1)
else:
    mask_within_bins = (cat['pz_bin'] <= 4) & (cat['pz_bin'] >= 1)

# loop over the photo-z bins
for pz_bin in bins:
    
    # Bin Cut
    if want_dr10:
        mask_bin = cat['DR10_pz_bin'] == pz_bin # dont touch
    else:
        mask_bin = cat['pz_bin'] == pz_bin
    cat_pz_bin = cat[mask_bin]
    
    # Converting to Pandas and rename the columns
    df = cat_pz_bin.to_pandas()
    if want_dr10:
        df2 = df.loc[:, ('RA', 'DEC', 'DR10_Z_PHOT_MEDIAN_I')]#, 'PHOTSYS')]
        df2.rename(columns={'DR10_Z_PHOT_MEDIAN_I':'Z'}, inplace=True)    
    else:
        df2 = df.loc[:, ('RA', 'DEC', 'Z_PHOT_MEDIAN')]#, 'PHOTSYS')] 
        df2.rename(columns={'Z_PHOT_MEDIAN':'Z'}, inplace=True)    
    
    # New Columns with Mstellar given from above
    colnames = ['coordX', 'coordY', 'coordZ', 'dX', 'dY',
                'dZ', 'dXKaiser', 'dYKaiser', 'dZKaiser',
                'vX', 'vY', 'vZ', 'vR', 'vTheta', 'vPhi']

    # set to 0
    for col in colnames:
        df2[col] = 0

    # set the stellar mass
    df2['MStellar'] = 10.**cat_pz_bin['LOGM'] #MStellar 

    # mask for that bin for the cleaned catalogs (needed only for unit_los)
    mask_clean = pz_bin_clean == pz_bin
    unit_los_pz_bin = unit_los[mask_clean]
    assert len(cat_pz_bin) == np.sum(mask_clean)
    assert np.all(mask_bin == mask_clean)
    
    # initialize displacements
    if recon_all_pz_bin:
        pz_str = ""
        displacements = np.zeros_like(unit_los)
        unit_los_sanity = np.zeros_like(unit_los)
    else:
        pz_str = f"_bin_{pz_bin:d}"
        displacements = np.zeros_like(unit_los_pz_bin)
        unit_los_sanity = np.zeros_like(unit_los_pz_bin)
        mask_clean = np.ones(len(cat_pz_bin), dtype=bool)        

    # are we combining the footprints or running them one by one
    if comb_footprint:
        footprints = [None]
    else:
        if want_dr10:
            footprints = ["DES", "S"]
        else:
            footprints = ["DES", "N", "S"]
    print("cat", len(cat), len(cat_pz_bin))
    
    # loop over footprints
    for want_footprint in footprints:
        print("footprint", want_footprint)
        if want_footprint is not None:
            if recon_all_pz_bin:
                cat_hp_idx = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True, nest=False)
            else:
                cat_hp_idx = hp.pixelfunc.ang2pix(nside, cat_pz_bin['RA'], cat_pz_bin['DEC'], lonlat=True, nest=False)
        else:
            if recon_all_pz_bin:
                mask_foot = np.ones(len(cat), dtype=bool)
            else:
                mask_foot = np.ones(len(cat_pz_bin), dtype=bool)
        if want_footprint == "DES":
            if recon_all_pz_bin:
                mask_foot = np.in1d(cat_hp_idx, des_hp_idx) & (cat['PHOTSYS'] == "S")
            else:
                mask_foot = np.in1d(cat_hp_idx, des_hp_idx) & (cat_pz_bin['PHOTSYS'] == "S") 
            foot_str = "_des"
        elif want_footprint == "N":
            if recon_all_pz_bin:
                mask_foot = cat['PHOTSYS'] == "N"
            else:
                mask_foot = cat_pz_bin['PHOTSYS'] == "N"
            foot_str = "_north"
        elif want_footprint == "S":
            if recon_all_pz_bin:
                mask_foot = ~np.in1d(cat_hp_idx, des_hp_idx) & (cat['PHOTSYS'] == "S")
            else:
                mask_foot = ~np.in1d(cat_hp_idx, des_hp_idx) & (cat_pz_bin['PHOTSYS'] == "S")
            foot_str = "_south"
        else:
            foot_str = ""
        if recon_all_pz_bin:# and want_footprint is not None:
            mask_foot &= mask_within_bins
        
        
        # load the file with reconstructed displacements
        recon_fn = f"displacements_{gal_fn.split('.fits')[0]}_{rand_fn.split('.fits')[0]}{mask_str}_R{sr:.2f}_nmesh{nmesh:d}_{convention}_{rectype}{pz_str}{foot_str}.npz"
        data = np.load(recon_dir / "recon" / recon_fn)
        displacements[mask_foot] = data['displacements']
        unit_los_sanity[mask_foot] = data['unit_los'] 

        # sooooo how about if want_mix load here the displacements (and also the RA DEC) and then match (but how do we match....)
        # attach dr10
        # recon_fn = f"displacements_{gal_fn.split('.fits')[0]}_{rand_fn.split('.fits')[0]}{mask_str}_R{sr:.2f}_nmesh{nmesh:d}_{convention}_{rectype}{pz_str}{foot_str}.npz"
        
    # apply bin mask if recon_all_pz_bin or do nothing
    displacements = displacements[mask_clean]
    unit_los_sanity = unit_los_sanity[mask_clean]
    print("please 0", np.sum(np.isnan(displacements)))
    if np.sum(np.isnan(unit_los_pz_bin)) > 0:
        unit_los_pz_bin[np.isnan(unit_los_pz_bin)] = unit_los_sanity[np.isnan(unit_los_pz_bin)]
    assert np.isclose(np.sum(unit_los_sanity-unit_los_pz_bin), 0.)
    del unit_los_sanity; gc.collect()
    
    # recon scalars
    mean_z = data['mean_z']
    growth_factor = data['growth_factor']
    Hubble_z = data['Hubble_z']
    h = data['h']
    del data; gc.collect()
    
    # calculate reconstructed velocities
    Velocity, Velocity_sphere = calc_velocity(displacements, unit_los_pz_bin, 1./(1.+mean_z), growth_factor, Hubble_z, h, want_rsd=False)
    df2['vX'] = Velocity[:, 0]
    df2['vY'] = Velocity[:, 1]
    df2['vZ'] = Velocity[:, 2]

    # vTheta and vPhi are incorrect, but it doesn't matter as long as vR is fine
    df2['vR'] = Velocity_sphere[:, 2]
    df2['vTheta'] = Velocity_sphere[:, 0]
    df2['vPhi'] = Velocity_sphere[:, 1]

    # if selecting special footprint
    if only_footprint is not None:
        if only_footprint == "DES":
            cat_hp_idx = hp.pixelfunc.ang2pix(nside, df2['RA'], df2['DEC'], lonlat=True, nest=False)
            mask_foot = np.in1d(cat_hp_idx, des_hp_idx) & (cat_pz_bin['PHOTSYS'] == "S")
        elif only_footprint == "S":
            cat_hp_idx = hp.pixelfunc.ang2pix(nside, df2['RA'], df2['DEC'], lonlat=True, nest=False)
            mask_foot = ~np.in1d(cat_hp_idx, des_hp_idx) & (cat_pz_bin['PHOTSYS'] == "S")
        elif only_footprint == "N":
            mask_foot = cat_pz_bin['PHOTSYS'] == "S"
        df2 = df2[mask_foot]

    # perturb positions by +/- 5 deg
    if doRandomPositions:
        df2['DEC'] += (np.random.rand(len(df2['RA']))-0.5)*10.
        df2['DEC'][df2['DEC'] > 90.] -= 90.
        df2['DEC'][df2['DEC'] <= -90.] += 90.
        df2['RA'] += (np.random.rand(len(df2['RA']))-0.5)*10.
        df2['RA'] %= 360.
        
    # interpolate the map to the given sky coordinates
    sourcecoord = np.array([df2['DEC'], df2['RA']]) * np.pi/180.   # convert from degrees to radians

    # use nearest neighbor interpolation # cmb_map
    df2['vZ'] = cmb_map.at(sourcecoord, prefilter=False, mask_nan=False, order=0) # note we are overwriting
    print("number of galaxies in bin", sourcecoord.shape[1])
    del sourcecoord; gc.collect()
    print(df2['vZ'][df2['vZ'] != 0.][:10])
    print(np.std(df2['vZ'][df2['vZ'] != 0.]))
    
    # TESTING!!!!!!!! # I THINK TECHNICALLY SHOULD BE DONE WITHIN THE FORLOOP
    # could be load evectors 
    """
    from cosmoprimo.fiducial import Planck2018FullFlatLCDM, DESI
    from pyrecon import  utils
    cosmo = DESI()
    Position = utils.sky_to_cartesian(cosmo.comoving_radial_distance(Z), RA, DEC)

    from cosmoprimo.fiducial import Planck2018FullFlatLCDM, DESI
    cosmo = DESI()
    position = unit_los_pz_bin*cosmo.comoving_radial_distance(df2['Z']) # should be Mpc and doesn't need to use cosmoprimo
    recon_fn = f"eigenvals_eigenvecs_{gal_fn.split('.fits')[0]}_{rand_fn.split('.fits')[0]}{mask_str}_R{sr:.2f}_nmesh{nmesh:d}_{convention}_{rectype}{pz_str}{foot_str}.npz"
    data = np.load(recon_dir / "recon" / recon_fn)
    data = np.load(fn)
    evals = data['evals']
    evecs = data['evecs']
    boxcenter_cellunits = data['boxcenter_cellunits']
    cellsize = boxsize/(nmesh-1)
    position /= cellsize
    position -= boxcenter_cellunits
    position = position.astype(np.int64)
    i_argmin = np.argmin(evals, axis=3)

    @njit
    def cross(a, b):
        c = np.array([a[1]*b[2] - a[2]*b[1],
                      a[2]*b[0] - a[0]*b[2],
                      a[0]*b[1] - a[1]*b[0]])
        return c

    #e3 = np.array([0., 0., 0.], dtype=np.float32)
    for i in range(position.shape[0]):
        e3 = evecs[position[i, 0], position[i, 1], position[i, 2], :, i_argmin[position[i, 0], position[i, 1], position[i, 2]]]
    """

    if not skip_tidal:
        #recon_fn = f"tidal_field_2D_{gal_fn.split('.fits')[0]}_{rand_fn.split('.fits')[0]}{mask_str}_R0.30_nmesh512{pz_str}{foot_str}.npz"
        recon_fn = f"tidal_field_2D_{gal_fn.split('.fits')[0]}_{rand_fn.split('.fits')[0]}{mask_str}_R0.06_nmesh1024{pz_str}{foot_str}.npz"
        print(str(recon_dir / "recon" / recon_fn))
        data = np.load(recon_dir / "recon" / recon_fn)
        ca = data['ca'] # mirror? need to test
        sa = data['sa'] 
        assert len(ca) == len(sa) == len(df2['RA'])
        df2['vX'] = ca # cos alpha (e_th, e2)
        df2['vY'] = sa # sin alpha (e_th, e2) --> e2 . (R e_th) = 1, R [[ca -sa], [sa ca]]
    
    # Saving the fits catalogs as txt files
    df2.to_csv(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}', header=None, index=None, sep=' ', mode='w')
    print(f'Saved DESI_pz{pz_bin:d}')

# TODO: do we need to update this?
# cosmological parameters
u = UnivMariana()

# Lastly we read & save these again using the built in ThumbStack catalog function, so it adds the last few columns (like mVir) for us.
for pz_bin in bins:
     Catalog(u, massConversion, name=f"DESI_pz{pz_bin:d}", nameLong=f"DESI pz bin {pz_bin:d}", pathInCatalog=str(desi_dir / f'DESI_pz{pz_bin:d}/{save_fn}'), save=True, cat_fn=save_fn)
