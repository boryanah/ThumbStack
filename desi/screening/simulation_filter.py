import gc, sys

from pixell import enmap, enplot, reproject, curvedsky as cs, utils
import numpy as np

import pandas as pd
from astropy.io import fits

import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def eshow(x, fn, **kwargs):
    ''' Define a function to help us plot the maps neatly '''
    plots = enplot.get_plots(x, **kwargs)
    #enplot.show(plots, method = "python")
    enplot.write(fn, plots)

"""
# TNG300
python simulation_filter.py 69 100 Y
python simulation_filter.py 62 100 Y
python simulation_filter.py 56 100 Y
python simulation_filter.py 52 100 Y

python simulation_filter.py 69 100 tau
python simulation_filter.py 62 100 tau
python simulation_filter.py 56 100 tau
python simulation_filter.py 52 100 tau

python simulation_filter.py 69 100 dm
python simulation_filter.py 62 100 dm
python simulation_filter.py 56 100 dm
python simulation_filter.py 52 100 dm

python simulation_filter.py 69 100 b
python simulation_filter.py 62 100 b
python simulation_filter.py 56 100 b
python simulation_filter.py 52 100 b

# Illustris
python simulation_filter.py 105 100 Y
python simulation_filter.py 98 100 Y
python simulation_filter.py 92 100 Y
python simulation_filter.py 88 100 Y

python simulation_filter.py 105 100 tau
python simulation_filter.py 98 100 tau
python simulation_filter.py 92 100 tau
python simulation_filter.py 88 100 tau

python simulation_filter.py 105 100 dm
python simulation_filter.py 98 100 dm
python simulation_filter.py 92 100 dm
python simulation_filter.py 88 100 dm

python simulation_filter.py 105 100 b
python simulation_filter.py 98 100 b
python simulation_filter.py 92 100 b
python simulation_filter.py 88 100 b
"""

# parameters
#sim_name = "MTNG"
#sim_name = "TNG300"
#sim_name = "TNG100"
sim_name = "Illustris"
if sim_name == "MTNG":
    Lbox = 500. # Mpc/h
    field_dir_fp = "/mnt/alan1/boryanah/MTNG/data_fp/"
    data_dir = "/mnt/alan1/boryanah/MTNG/data_sz/"
elif sim_name == "TNG300":
    Lbox = 205. # Mpc/h
    field_dir_fp = "/pscratch/sd/b/boryanah/TNG/" # not used
    data_dir = "/pscratch/sd/b/boryanah/TNG/"
elif sim_name == "TNG100":
    Lbox = 75. # Mpc/h
    field_dir_fp = "/mnt/gosling1/boryanah/TNG100/"
    data_dir = "/mnt/marvin1/boryanah/SZ_TNG/TNG100/"
elif sim_name == "Illustris":
    Lbox = 75. # Mpc/h
    field_dir_fp = "/mnt/gosling1/boryanah/Illustris/"
    data_dir = "/pscratch/sd/b/boryanah/Illustris/" #"/mnt/marvin1/boryanah/SZ_TNG/Illustris/"
nu = int(sys.argv[2]) # GHz
if np.isclose(nu, 150):
    beam_fwhm = 1.3
elif np.isclose(nu, 100):
    beam_fwhm = 1.6
elif np.isclose(nu, 90):
    beam_fwhm = 2.1
if np.isclose(nu, 0):
    beam_str = ""
else:
    beam_str = f"_beam{beam_fwhm:.1f}"
    
fn_base_Y = f"Y_compton{beam_str}_xy"
fn_base_tau = f"tau{beam_str}_xy"
fn_base_b = f"b{beam_str}_xy"
fn_base_dm = f"dm{beam_str}_xy"


snap = int(sys.argv[1])
if sim_name == "MTNG":
    if snap == 264:
        z = 0.5 # 0. TESTING!!!!!!!!!! shouldn't matter?
    elif snap == 237:
        z = 0.25
    elif snap == 214:
        z = 0.5
    elif snap == 179:
        z = 1.
elif sim_name == "TNG300" or sim_name == "TNG100": 
    if snap == 67:
        z = 0.5
    elif snap == 69:
        z = 0.47 
    elif snap == 62:
        z = 0.629
    elif snap == 56:
        z = 0.791
    elif snap == 52:
        z = 0.924
elif sim_name == "Illustris":
    if snap == 105:
        z = 0.47 
    elif snap == 98:
        z = 0.629
    elif snap == 92:
        z = 0.791
    elif snap == 88:
        z = 0.924
a = 1./(1.+z)
if sim_name == "Illustris":
    h = 0.704
    Om_m = 0.2726
else:
    h = 0.6774
    Om_m = 0.3089

# define cosmology
cosmo = FlatLambdaCDM(H0=h*100., Om0=Om_m, Tcmb0=2.725)

# compute angular distance
d_L = cosmo.luminosity_distance(z).to(u.Mpc).value
d_C = d_L/(1.+z) # dC = dL/(1+z) Mpc
d_A = d_L/(1.+z)**2 # dA = dL/(1+z)^2 # Mpc
print("d_A [Mpc], z", z, d_A)
d_A *= h # Mpc/h
d_C *= h # Mpc/h
print("comoving distance = ", d_C)

# get size on the sky of each pixel at given redshift
Lbox_deg = (a*Lbox)/d_A*(180./np.pi) # degrees
print("Lbox deg = ", Lbox_deg)

# cell size
if sim_name == "Illustris" or sim_name == "TNG100":
    N_cell = 2000
else:
    N_cell = 10000
cell_size = Lbox/N_cell # cMpc/h
cell_size_deg = Lbox_deg/N_cell
print("cell size arcmin", cell_size_deg*60.)

# cell size dm
if sim_name == "Illustris" or sim_name == "TNG100":
    N_cell_dm = 2000
else:
    N_cell_dm = 10000
cell_size_dm = Lbox/N_cell_dm # cMpc/h
cell_size_dm_deg = Lbox_deg/N_cell_dm
print("cell size dm arcmin", cell_size_dm_deg*60.)
sys.stdout.flush()

#map_type = "Y"
#map_type = "tau"
#map_type = "b"
#map_type = "dm"
map_type = sys.argv[3]

"""
# plot map
ACT_map_sm = enmap.read_map(f"{data_dir}/{map_type}{beam_str}_snap_{snap:d}_small_tau_screening.fits")
ACT_map_lg = enmap.read_map(f"{data_dir}/{map_type}{beam_str}_snap_{snap:d}_large_tau_screening.fits")
eshow(ACT_map_sm, f"{map_type}{beam_str}_snap_{snap:d}_small_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
eshow(ACT_map_lg, f"{map_type}{beam_str}_snap_{snap:d}_large_tau_screening", **{"colorbar":True, "ticks": 5, "downgrade": 4})
quit()
"""

if map_type == "dm":
    want_dm = True
    want_sz = False
else:
    want_dm = False
    want_sz = True


# create pixell map
if want_sz:
    box = np.array([[0., 0.],[Lbox_deg, Lbox_deg]]) * utils.degree
    shape, wcs = enmap.geometry(pos=box, res=cell_size_deg * utils.degree, proj='car')
    Y_xy_map = enmap.zeros(shape, wcs=wcs)
    tau_xy_map = enmap.zeros(shape, wcs=wcs)
    b_xy_map = enmap.zeros(shape, wcs=wcs)

# create pixell map dm
if want_dm:
    box = np.array([[0., 0.],[Lbox_deg, Lbox_deg]]) * utils.degree
    shape, wcs = enmap.geometry(pos=box, res=cell_size_dm_deg * utils.degree, proj='car')
    dm_xy_map = enmap.zeros(shape, wcs=wcs)

# load MTNG
if want_sz:
    Y_xy_map[:] = np.load(f"{data_dir}/{fn_base_Y}_snap_{snap:d}.npy")
    tau_xy_map[:] = np.load(f"{data_dir}/{fn_base_tau}_snap_{snap:d}.npy")
    b_xy_map[:] = np.load(f"{data_dir}/{fn_base_b}_snap_{snap:d}.npy")
    assert N_cell == Y_xy_map.shape[0] == Y_xy_map.shape[1]
    msk = Y_xy_map.copy()
    msk *= 0.
    msk += 1.
if want_dm:
    dm_xy_map[:] = np.load(f"{data_dir}/{fn_base_dm}_snap_{snap:d}.npy")
    msk = dm_xy_map.copy()
    msk *= 0.
    msk += 1.
    assert N_cell_dm == dm_xy_map.shape[0] == dm_xy_map.shape[1]

# read
if map_type == "Y":
    ACT_map = Y_xy_map
    ACT_mask = msk
elif map_type == "tau":
    ACT_map = tau_xy_map
    ACT_mask = msk
elif map_type == "b":
    ACT_map = b_xy_map
    ACT_mask = msk
elif map_type == "dm":
    ACT_map = dm_xy_map
    ACT_mask = msk
    
# mask
ACT_map *= ACT_mask
shape, wcs = ACT_map.shape, ACT_map.wcs
del ACT_mask
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

"""
# convert to alms (curved sky)
ACT_alms = cs.map2alm(ACT_map, lmax=10000)
del ACT_map
gc.collect()

# We can now filter the kappa alms by the filter produced above
ACT_alms_sm = hp.almxfl(ACT_alms, fsm_ksz)
ACT_alms_lg = hp.almxfl(ACT_alms, flg_ksz)

# We now convert the filtered alms to map space
ACT_map_sm = cs.alm2map(ACT_alms_sm, enmap.empty(shape, wcs))
ACT_map_lg = cs.alm2map(ACT_alms_lg, enmap.empty(shape, wcs))
"""

kmap = enmap.fft(ACT_map, normalize="phys")
#lmap = enmap.lmap(shape, wcs)
#lymap, lxmap = lmap
modlmap = enmap.modlmap(shape, wcs)

flg_map = modlmap.copy()
flg_map[:] = 0.
flg_map[modlmap < 2000.] = 1.
choice = (modlmap >= 2000.) & (modlmap < 2150.)
flg_map[choice] = np.cos((modlmap[choice] - 2000.)*np.pi/300.)

fsm_map = modlmap.copy()
fsm_map[:] = 0.
fsm_map[modlmap > 2500.] = 1.
choice = (modlmap <= 2500.) & (modlmap > 2350.)
fsm_map[choice] = np.sin((modlmap[choice] - 2350.)*np.pi/300.)

ACT_map_sm = enmap.ifft(kmap*fsm_map, normalize="phys").real
ACT_map_lg = enmap.ifft(kmap*flg_map, normalize="phys").real

print(ACT_map_sm.shape)
print(ACT_map_lg)

# write map
enmap.write_map(f"{data_dir}/{map_type}{beam_str}_snap_{snap:d}_small_tau_screening.fits", ACT_map_sm)
enmap.write_map(f"{data_dir}/{map_type}{beam_str}_snap_{snap:d}_large_tau_screening.fits", ACT_map_lg)
