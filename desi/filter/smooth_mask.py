from pixell import enmap, utils, powspec, enplot, reproject #, pointsrcs
import healpy as hp
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter

pscratch = '/pscratch/sd/r/rhliu/projects'

pathMap = pscratch + '/ThumbStack/ACT_DR6/ilc_SZ_yy.fits'
pathMask = pscratch + '/ThumbStack/ACT_DR6/wide_mask_GAL070_apod_1.50_deg_wExtended.fits'

ACT_map = enmap.read_map(pathMap)
ACT_mask = enmap.read_map(pathMask)

with fits.open(pathMap) as hdul:
    res = np.abs(hdul[0].header['CDELT1'])
    assert res == np.abs(hdul[0].header['CDELT2'])
res_rad = np.deg2rad(res)

pathUnionLarge = pscratch + '/ThumbStack/cluster_catalogs/union_catalog_large_20220316.csv'
pathUnionRegular = pscratch + '/ThumbStack/cluster_catalogs/union_catalog_regular_20220316.csv'

pathClusterLarge = pscratch + '/ThumbStack/cluster_catalogs/large_cluster_catalog_20220316.csv'
pathClusterRegular = pscratch + '/ThumbStack/cluster_catalogs/regular_cluster_catalog_20220316.csv'
savePath = pscratch + '/ThumbStack/cluster_catalogs/'

paths = [pathUnionLarge, pathUnionRegular, pathClusterLarge, pathClusterRegular]
radii = [10, 3, 5, 3] # In arcminutes, radius to mask for each type of source
saveNames = ['mask_union_catalog_large.fits', 'mask_union_catalog_regular.fits', 'mask_large_cluster.fits', 'mask_regular_cluster.fits']

shape = ACT_map.shape
wcs = ACT_map.wcs
def Gaussian(x, A, mu, sigma):
    return A * np.exp(-1/2*(x-mu)**2/sigma**2)
threshold = Gaussian(1,1,0,1) # 1 sigma y threshold for a gaussian

arcmin = 0.000290888 # radians
factor = np.round(arcmin / res_rad, 5) # conversion factor from arcmin radius to pixel radius, should be 2 (since resolution is 0.5 arcmin)
listofMasks = []

for i, path in enumerate(paths):
    catalog = pd.read_csv(path, skiprows = 3)
    RA = np.array(catalog['# RA (deg)'])
    DEC = np.array(catalog[' Dec (deg)'])
    coords = np.deg2rad(np.vstack((DEC, RA)))
    
    Map2 = ACT_map.copy() * 0
    
    ypix,xpix = enmap.sky2pix(shape,wcs,coords)
    ypix = np.round(ypix).astype(int)
    xpix = np.round(xpix).astype(int)
    
    Map2[ypix, xpix] = 1
    
    Map_gaussian = gaussian_filter(Map2, sigma=radii[i] * factor)
    Map_gaussian = Map_gaussian / Map_gaussian.max()
    
    ACT_filter = (Map_gaussian < threshold).astype(np.float32)
    new_ndmap = enmap.ones(shape, wcs, dtype=np.float32)
    ACT_filter2 = (new_ndmap * ACT_filter)
    print('saving map: ' + saveNames[i])
    enmap.write_map(savePath + saveNames[i], ACT_filter2)
    listofMasks.append(ACT_filter2)
    
    
# # The next part of the code is optional: We use this part to combine the masks into a single mask alongside the 70% Galaxy Mask
# final_mask = ACT_mask.copy()
# for mask in listofMasks:
#     final_mask = final_mask * mask
# enmap.write_map(pscratch + '/ThumbStack/ACT_DR6/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree.fits', final_mask)
print('Done!!')
