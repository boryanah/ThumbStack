import sys
sys.path.append('/global/homes/b/boryanah/repos/ThumbStack')

from importlib import reload
import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

import thumbstack
reload(thumbstack)
from thumbstack import *

import cmb
reload(cmb)
from cmb import *
from cmbMap import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

##################################################################################

"""
# NOTE THE CHOICE OF MAX_SIGMAZ

python run_thumbstack.py 1 1 
python run_thumbstack.py 0 0 
python run_thumbstack.py 1 0 # best?
python run_thumbstack.py 0 1 

python run_thumbstack.py 1 0 DES
python run_thumbstack.py 1 0 N
python run_thumbstack.py 1 0 S

python run_thumbstack.py 1 0 extended dr9
python run_thumbstack.py 1 0 main dr9
python run_thumbstack.py 1 0 extended dr10
python run_thumbstack.py 1 0 main dr10

python run_thumbstack.py 1 0 bgs na f090 # not applicable
python run_thumbstack.py 1 0 bgs na f150 # not applicable
python run_thumbstack.py 1 0 bgs_lrg_elg_qso na # not applicable
python run_thumbstack.py 1 0 des_y6 na # not applicable
python run_thumbstack.py 1 0 LS na # not applicable

python run_thumbstack.py 1 0 halos lensed
python run_thumbstack.py 1 0 halos unlensed
python run_thumbstack.py 1 0 halos lensed_abacus
python run_thumbstack.py 1 0 halos lensed_abacus_fwhm1.6_noise15.0
python run_thumbstack.py 1 0 halos lensed_abacus_small_large
python run_thumbstack.py 1 0 halos lensed_tau_abacus
python run_thumbstack.py 1 0 halos unlensed_tau_abacus
python run_thumbstack.py 1 0 halos tau_abacus
python run_thumbstack.py 1 0 halos lensed_abacus_fwhm1.6
python run_thumbstack.py 1 0 halos lensed_tau_abacus_fwhm1.6
python run_thumbstack.py 1 0 halos unlensed_tau_abacus_fwhm1.6
python run_thumbstack.py 1 0 halos tau_abacus_fwhm1.6

python run_thumbstack.py 1 0 halos lensed_abacus_fwhm2.1
python run_thumbstack.py 1 0 halos unlensed_tau_abacus_fwhm2.1

python run_thumbstack.py 1 0 halos lensed_abacus_fwhm1.4
python run_thumbstack.py 1 0 halos lensed_order1_abacus_fwhm1.4
python run_thumbstack.py 1 0 halos kappa_order1_abacus
python run_thumbstack.py 1 0 halos unlensed_tau_abacus_fwhm1.4
python run_thumbstack.py 1 0 halos unlensed_abacus_fwhm1.4

python run_thumbstack.py 1 0 extended dr9 f090
python run_thumbstack.py 1 0 extended dr9 f150
python run_thumbstack.py 1 0 extended dr9 kappa
"""

nProc = 256 # this is full
#nProc = 128 # this is half
#nProc = 64 # this is quarter

#version = "dr5"
version = "dr6"
#version = "websky"
freq = "f090"
if version == "dr5":
    #pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree.fits' # old mask
    pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
    freq = sys.argv[5]
    #freq = "f150"
    #freq = "f090"
    if freq == "f090":
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked_rebeamed_f150.fits" # 2.1 arcmin
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked.fits" # 2.1 arcmin
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked_rebeamed_f150_small_tau_screening.fits" # 2.1 arcmin
        pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f090_night_map_srcfree_masked_rebeamed_f150_small_large_tau_screening.fits" # 2.1 arcmin
        pathHit = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_dr5.01_s08s18_AA_f090_night_ivar.fits"
        freq_str = "_f090"
    elif freq == "f150":
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_masked.fits" # 1.3 arcmin
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_masked_small_tau_screening.fits" # 1.3 arcmin
        pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_planck_dr5.01_s08s18_AA_f150_night_map_srcfree_masked_small_large_tau_screening.fits" # 1.3 arcmin
        pathHit = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/act_dr5.01_s08s18_AA_f150_night_ivar.fits"
        freq_str = "_f150"
    #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/cmb_night_tot_f090_uncorr_coadd_map.fits" # 2.1 arcmin
    pathHit = None # TESTING note we are doing this so the units make sense
    version_str = ""
    if freq == "f090":
        version_str_bgs = "_dr5"
    elif freq == "f150":
        version_str_bgs = "_dr5_f150"
        version_str = "_f150"
    else:
        print("sths not right"); quit()
elif version == "dr6":
    pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
    if len(sys.argv) > 5:
        if sys.argv[5] == "kappa":
            pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/maps/baseline/kappa_alm_data_act_dr6_lensing_v1_baseline_masked_bias_tau_screening.fits" # TESTING
            #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/maps/baseline/kappa_alm_data_act_dr6_lensing_v1_baseline_masked.fits" # TESTING
            #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/maps/baseline/kappa_alm_data_act_dr6_lensing_v1_baseline_masked_tau_screening.fits" # TESTING
    else:
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits" # 1.6 arcmin # OG
        pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000_half_small_large_small_tau_screening.fits" # 1.6 arcmin # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/ilc_fullRes_TT_half_small_large_small_tau_screening.fits" # 1.6 arcmin # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000_small_large_tau_screening.fits" # 1.6 arcmin # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/ilc_actplanck_ymap.fits" # 1.6 arcmin # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/evals_1_bin_1.fits" # eigenvalue # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000_kSZ_filter.fits" # 1.6 arcmin # TESTING
        #pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000_small_tau_screening.fits" # 1.6 arcmin # TESTING
    version_str = "_dr6"
    if len(sys.argv) > 5:
        if sys.argv[5] == "kappa":
            version_str += "_kappa_bias"
    version_str_bgs = ""
    freq_str = ""
    pathHit = None
elif version == "websky":
    version_str = "_websky"
    version_str_bgs = ""
    freq_str = ""
    pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
    if "fwhm" in sys.argv[4]:
        fwhm = np.float((sys.argv[4].split('_fwhm')[-1]).split('_')[0])
    if "abacus" in sys.argv[4]:
        map_dir = "/pscratch/sd/b/boryanah/websky/abacus/"
        if "tau_abacus_fwhm" in sys.argv[4] and "noise15.0" in sys.argv[4]:
            map_str = f"_tau_8192_ph201_fwhm{fwhm:.1f}_noise15.0"
        elif "tau_abacus_fwhm" in sys.argv[4]:
            map_str = f"_tau_8192_ph201_fwhm{fwhm:.1f}"
        elif "fwhm" in sys.argv[4] and "noise15.0" in sys.argv[4]:
            map_str = f"_8192_ph201_fwhm{fwhm:.1f}_noise15.0"
        elif "fwhm" in sys.argv[4]:
            map_str = f"_8192_ph201_fwhm{fwhm:.1f}"
        elif "tau" in sys.argv[4]:
            map_str = "_tau_8192_ph201"
        else:
            map_str = "_8192_ph201"
        if "order" in sys.argv[4]:
            #map_str = map_str.replace("ph201_", "ph201_order1_") # TESTING!!!!!!!!!!!!!!
            map_str = map_str.replace("ph201_", "ph201_order1_brute_") # TESTING!!!!!!!!!!!!!!
        if "unlensed" in sys.argv[4] and "tau" not in sys.argv[4]:
            map_str = map_str.replace("_ph201", "")
    else:
        map_dir = "/pscratch/sd/b/boryanah/websky/"
        map_str = ""
        
    if "unlensed" in sys.argv[4]:
        pathMap = f'{map_dir}/unlensed_map{map_str}_small_tau_screening.fits'
    elif "lensed" in sys.argv[4]:
        if "small_large" in sys.argv[4]:
            pathMap = f'{map_dir}/lensed_map{map_str}_small_large_tau_screening.fits'
            version_str += "_small_large"
        else:
            pathMap = f'{map_dir}/lensed_map{map_str}_small_tau_screening.fits'
    else:
        #pathMap = f'{map_dir}/map{map_str}_small_tau_screening.fits' # this is pure tau
        pathMap = f'{map_dir}/kappa{map_str}_tau_screening.fits'#print("Not sure what you want"); quit()
        #pathMap = f'{map_dir}/map{map_str}_tau_screening.fits'#print("Not sure what you want"); quit()
    pathHit = None
#pathHit = "/pscratch/sd/r/rhliu/projects/ThumbStack/" + "act_dr5.01_s08s18_AA_f150_daynight_ivar.fits"
#pathHit = None
print(pathMap)

if "Will" not in pathMask:
    mask_str = "_oldmask"
else:
    mask_str = ""

if "kSZ_filter" in pathMap:
    filter_str = "_kSZ_filter"
else:
    filter_str = ""

if "lensing" in pathMap:
    lensing_str = "_lensing"
else:
    lensing_str = ""

if "tau_screening" in pathMap:
    screening_str = "_tau_screening"
else:
    screening_str = ""

if "tau_screening" in pathMap:
    filterType = "meanring" # og
    #filterType = "diskring" # TESTING!!!!!!!!!!!!
    Obs = 'tau'
elif "lensing" in pathMap:
    filterType = "meanring"
    Obs = 'tsz'
else:
    filterType = "diskring"
    #Obs = 'ksz' # og
    Obs = 'ksz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta
    #Obs = 'tsz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta

# Testing
want_all = False #True # combine all galaxies in the 4 redshift bins, default is False
save = True # if yes calculates delta T anew (saves mask, delta T, maybe vel is separate? but the numbers need to match)
#save = False # False means use saved stuff 
doBootstrap = True #False #True # do bootstrap or not? (not needed if e.g. shuffling velocities)
doMBins = False # we don't have mass bins yet
doVShuffle = False #True # default measurement is False
wantMF = False # default measurement is False
doOnlyFiltering = True #True # default is False; if True, compute dT decrements and return TESTING!!!!!!!!!!!; the only problem is that.... it doesn't compute covariance only stackedMap; we use this for anisotropic
want_mix = False #True # DR9 with DR10 z where possible
vshuff_str = "_vshuffle" if doVShuffle else ""
cat_type = sys.argv[3] #"extended", "main"
if cat_type == "halos":
    assert version == "websky"
cat_dr = sys.argv[4] # "dr9", "dr10" # default so far (want dr9 if mix)
if cat_dr == "dr10":
    cat_dr_str = "_dr10"
elif cat_dr == "lensed":
    cat_dr_str = "_lensed"
elif cat_dr == "lensed_abacus":
    cat_dr_str = "_lensed_abacus"
elif cat_dr == "lensed_abacus_small_large":
    cat_dr_str = "_lensed_abacus_small_large"
elif cat_dr == "lensed_abacus_fwhm1.6":
    cat_dr_str = "_lensed_abacus_fwhm1.6"
elif cat_dr == "lensed_abacus_fwhm2.1":
    cat_dr_str = "_lensed_abacus_fwhm2.1"
elif cat_dr == "lensed_abacus_fwhm1.4":
    cat_dr_str = "_lensed_abacus_fwhm1.4"
elif cat_dr == "lensed_order1_abacus_fwhm1.4":
    cat_dr_str = "_lensed_order1_abacus_fwhm1.4"
elif cat_dr == "kappa_order1_abacus":
    cat_dr_str = "_kappa_order1_abacus"
elif cat_dr == "lensed_abacus_fwhm1.6_noise15.0":
    cat_dr_str = "_lensed_abacus_fwhm1.6_noise15.0"
elif cat_dr == "lensed_tau_abacus":
    cat_dr_str = "_lensed_tau_abacus"
elif cat_dr == "unlensed_tau_abacus":
    cat_dr_str = "_unlensed_tau_abacus"
elif cat_dr == "tau_abacus":
    cat_dr_str = "_tau_abacus"
elif cat_dr == "lensed_tau_abacus_fwhm1.6":
    cat_dr_str = "_lensed_tau_abacus_fwhm1.6"
elif cat_dr == "unlensed_tau_abacus_fwhm1.6":
    cat_dr_str = "_unlensed_tau_abacus_fwhm1.6"
elif cat_dr == "unlensed_tau_abacus_fwhm2.1":
    cat_dr_str = "_unlensed_tau_abacus_fwhm2.1"
elif cat_dr == "unlensed_tau_abacus_fwhm1.4":
    cat_dr_str = "_unlensed_tau_abacus_fwhm1.4"
elif cat_dr == "unlensed_abacus_fwhm1.4":
    cat_dr_str = "_unlensed_abacus_fwhm1.4"
elif cat_dr == "tau_abacus_fwhm1.6":
    cat_dr_str = "_tau_abacus_fwhm1.6"
elif cat_dr == "unlensed":
    cat_dr_str = "_unlensed"
else:
    cat_dr_str = ""
mix_str = "_mix" if want_mix else ""
doRandomPositions = False
smooth = 0.03 #0.06 #0.18 # deg 0.06 is def
smooth_str = f"_R{smooth:.2f}" if smooth != 0.06 else ""
if doRandomPositions:
    random_str = "_randompos"
else:
    random_str = ""

# pick type of catalog
comb_footprint =  bool(int(sys.argv[1])) 
cat_foot_str = "_allfoot" if comb_footprint else "_joinfoot"

# recon params
recon_all_pz_bin = bool(int(sys.argv[2])) 
recon_bin_str = "_allbin" if recon_all_pz_bin else '_perbin'

# apply cut in redshift accuracy (this is done elsewhere and here we are simply loading the correct file)
max_sigmaz = 0.0 #0.05 # 0.0 #0.05 #0.035 # 0 or 0.05 is used in analysis
sigmaz_str = f"_sigmaz{max_sigmaz:.4f}" if not np.isclose(max_sigmaz, 0.) else ""

# feeling lazy
"""
# just one of the 3 regions
if len(sys.argv) > 3:
    only_str = f"_{sys.argv[3]}" # DES, N, S
else:
    only_str = ""
"""
only_str = ""

# catalog name
if cat_type == "main":
    cat_fn = f"catalog{mix_str}{cat_dr_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{freq_str}.txt" 
elif cat_type == "extended":
    cat_fn = f"extended_catalog{mix_str}{cat_dr_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{freq_str}.txt" 
elif cat_type == "bgs":
    #cat_fn = f"BGS_BRIGHT-21.5.txt" # spec # i think smolest? # TESTING just cause fast
    #cat_fn = f"BGS_BRIGHT.txt" # new
    #cat_fn = f"BGS_BRIGHT_pz.txt" # newest0
    cat_fn = f"BGS_BRIGHT_pz_logm10.5{random_str}{smooth_str}.txt" # newestest DEFAULT og
    tracer = cat_fn.split(".txt")[0]
elif cat_type == "bgs_lrg_elg_qso":
    cat_fn = f"BGS_LRG_ELG_QSO.txt"
elif cat_type == "des_y6":
    perc = 10
    cat_fn = f"DESY6_ABSM{perc:d}.fits"
    cat_dr_str = "_des_y6"
elif cat_type == "LS":
    LOGM_MIN = 10.5
    PHOTO_ERR = 0.1
    cat_fn = f"sweep_LOGM{LOGM_MIN:.1f}_PHOTOZ{PHOTO_ERR:.1f}.fits"
    cat_dr_str = "_LS"
elif cat_type == "halos":
    #cat_fn = f"halo_1.e13{cat_dr_str}.txt"
    #cat_fn = f"halo_3.e12{cat_dr_str}.txt"
    cat_fn = f"halo_6.e12{cat_dr_str}.txt"
print(cat_fn)


if wantMF:
    from scipy.interpolate import CubicSpline
    data = np.load("matched_filter/data/hilc_fullRes_TT_17000.npz")
    inv_loginterp = CubicSpline(data['centers'], np.log(1./data['binned_power']))
    invPowerFunc = lambda ell: np.exp(inv_loginterp(ell))
    data = np.load("matched_filter/data/theory_profile_gaussian.npz")
    theta_arcmin = data['theta_arcmin']
    prof = data['prof']
    choice = theta_arcmin < 20.
    loginterp = CubicSpline(theta_arcmin[choice]*np.pi/(180.*60.), np.log(prof[choice]))
    filterFuncRad = lambda th_rad: np.exp(loginterp(th_rad))
    apod_pix = 20
else:
    invPowerFunc = None
    filterFuncRad = None
    apod_pix = 20

# virial mass
Mvir_Msun = 10**(13.4)/0.6774

# cosmological parameters
u = UnivMariana()

# M*-Mh relation
massConversion = MassConversionKravtsov14()
# massConversion.plot()

##################################################################################
# Galaxy Catalogs (from DESI)

print("Read galaxy catalogs")
tStart = time()

# catalog name
cat_dir = '/pscratch/sd/b/boryanah/ACTxDESI/DESI/'
if cat_type == "bgs":
    output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output{version_str_bgs}_{tracer}{filter_str}{lensing_str}{screening_str}{smooth_str}/thumbstack/"
else:
    if "extended" in cat_fn:
        output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output_extended{mix_str}{cat_dr_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{vshuff_str}{version_str}{mask_str}{filter_str}{lensing_str}{screening_str}{smooth_str}/thumbstack/"
    else:
        output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output{mix_str}{cat_dr_str}{cat_foot_str}{recon_bin_str}{sigmaz_str}{only_str}{random_str}{vshuff_str}{version_str}{mask_str}{filter_str}{lensing_str}{screening_str}{smooth_str}/thumbstack/"

if cat_type == "bgs":
    catalogs = {
        "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
    }
elif cat_type == "bgs_lrg_elg_qso":
    catalogs = {
        "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
        "DESI_pz2": Catalog(u, massConversion, name="DESI_pz2", nameLong="DESI pz bin 2", out_dir=cat_dir, save=False, cat_fn=cat_fn),
        "DESI_pz3": Catalog(u, massConversion, name="DESI_pz3", nameLong="DESI pz bin 3", out_dir=cat_dir, save=False, cat_fn=cat_fn),
        "DESI_pz4": Catalog(u, massConversion, name="DESI_pz4", nameLong="DESI pz bin 4", out_dir=cat_dir, save=False, cat_fn=cat_fn),
    }
elif cat_type == "des_y6":
    catalogs = {
        "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
    }
elif cat_type == "LS":
    catalogs = {
        "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
    }
elif cat_type == "halos":
    catalogs = {
        "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
    }
else:
    if want_all:
        catalogs = {
            "DESI_all": Catalog(u, massConversion, name=["DESI_pz1", "DESI_pz2", "DESI_pz3", "DESI_pz4"], nameLong="DESI all bins", out_dir=cat_dir, save=False, cat_fn=cat_fn),
        }
    else:
        catalogs = {
            "DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
            "DESI_pz2": Catalog(u, massConversion, name="DESI_pz2", nameLong="DESI pz bin 2", out_dir=cat_dir, save=False, cat_fn=cat_fn),
            "DESI_pz3": Catalog(u, massConversion, name="DESI_pz3", nameLong="DESI pz bin 3", out_dir=cat_dir, save=False, cat_fn=cat_fn),
            "DESI_pz4": Catalog(u, massConversion, name="DESI_pz4", nameLong="DESI pz bin 4", out_dir=cat_dir, save=False, cat_fn=cat_fn),
        }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")
print("CMB_map", pathMap)

###################################################################################
# Read CMB maps

# Path List
CMB_hitlist = [pathHit]
CMB_pathlist = [pathMap]
CMB_masklist = [pathMask]
CMB_name = ['act_dr6_f90']
CMB_namepublic = ['ACT DR6 (90GHz)']
if freq == "f090":
    CMB_nu = [90.e9]
elif freq == "f150":
    CMB_nu = [150.e9]
else:
    CMB_nu = [90.e9]
cmbMap_list = []

for i, path in enumerate(CMB_pathlist):
    cmap = cmbMap(path,
                  pathMask=CMB_masklist[i],
                  pathHit=CMB_hitlist[i],
                  nu=CMB_nu[i], unitLatex=r'y',
                  name=CMB_name[i])
    cmbMap_list.append(cmap)

catalogKeys = catalogs.keys()
print(catalogKeys)

"""
# TESTING filter
from pixell import enmap
cmb_maps = []
if version == "websky":
    kappa_dir = "/pscratch/sd/b/boryanah/websky/abacus/"
else:
    kappa_dir = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/maps/baseline/"
for i in range(9):
    if version == "websky":
        #cmb_maps.append(enmap.read_map(kappa_dir + f"kappa_8192_ph201_tau_screening_diff_binc{i:d}.fits"))
        cmb_maps.append(enmap.read_map(kappa_dir + f"kappa_8192_ph201_tau_screening_same_binc{i:d}.fits"))
    else:
        cmb_maps.append(enmap.read_map(kappa_dir + f"kappa_alm_data_act_dr6_lensing_v1_baseline_masked_diff_binc{i:d}.fits"))
"""
 
###################################################################################
ts_list = [[] for _ in range(len(CMB_masklist))]

for key in list(catalogKeys):
    catalog = catalogs[key]
    #catalog.Mstellar = np.ones_like(catalog.Mstellar) * massConversion.fmVirTomStar(Mvir_Msun)
    catalog.Mstellar = np.empty_like(catalog.RA)
    catalog.Mvir = np.empty_like(catalog.RA)
    #catalog.Mvir = np.ones_like(catalog.Mvir) * Mvir_Msun
    #catalog.addIntegratedY()
    #catalog.addIntegratedTau()
    #catalog.addIntegratedKSZ()
    catalog.integratedY = np.empty_like(catalog.RA)
    catalog.integratedKSZ = np.empty_like(catalog.RA)
    catalog.integratedTau = np.empty_like(catalog.RA)

    if Obs == 'tau':
        catalog.vR = catalog.vZ # hiding here info about T_large-scales
        
    for i, cmap in enumerate(cmbMap_list):
    
        ts = ThumbStack(u, catalog, 
                        cmap.map(),  # tuksi
                        #cmb_maps,  # TESTING!!!! filter
                        cmap.mask(), 
                        cmap.hit(), 
                        catalog.name + '_' + cmap.name,
                        nameLong=None,
                        save=save, 
                        nProc=nProc,
                        filterTypes=filterType, 
                        doMBins=doMBins, 
                        doBootstrap=doBootstrap,
                        doVShuffle=doVShuffle,
                        doOnlyFiltering=doOnlyFiltering,
                        cmbNu=cmap.nu, 
                        cmbUnitLatex=cmap.unitLatex,
                        output_dir=output_dir,
                        Obs=Obs,
                        wantMF=wantMF,
                        invPowerFunc=invPowerFunc,
                        filterFuncRad=filterFuncRad,
                        apod_pix=apod_pix)
        ts_list[i].append(ts)


###################################################################################

# Next for plotting:

# the question is -- what do we need minimally for plotting and I guess recording the 2D plot
# if I am not mistaken you can do overlap flag and then go directly to computeStackedProfile with stackedMap=True

# Parameters
factor = (180.*60./np.pi)**2
if pathHit is not None:
    est = f'{Obs}_varweight'
else:
    est = f'{Obs}_uniformweight'
path = f"./figures/StackedMaps_{est}_{cat_type}{cat_dr_str}.png"

# ###############################
fig, subplots = plt.subplots(2, 2, figsize=(8,8), sharex='col', sharey='row')
subplots = subplots.ravel()

for i, key in enumerate(list(catalogKeys)):
    ax = subplots[i]
    for j in range(len(ts_list)):
        tsj = ts_list[j][i]

        stackedMap = tsj.computeStackedProfile(filterType, est, iBootstrap=None, iVShuffle=None, tTh='', stackedMap=True) # gives you a 2D map because of the stackedMap = True setting
        
        # this is empty and is shaped like cutout
        cutoutMap = tsj.cutoutGeometry()

        # size of canvas in radians 
        size = cutoutMap.posmap()[0,:,:].max() - cutoutMap.posmap()[0,:,:].min()

        # size of pixel in radians
        dx = float(size) / (cutoutMap.shape[0]-1) 
        dy = float(size) / (cutoutMap.shape[1]-1)

        # centers of the pixels (index+0.5 times size)
        x = dx * (np.arange(cutoutMap.shape[0]+1) - 0.5)
        y = dy * (np.arange(cutoutMap.shape[1]+1) - 0.5)
        x,y = np.meshgrid(x, y, indexing='ij')

        # save stackedMap filtertype cmap.name and obs or estimator rather
        np.savez(f"{output_dir}/{catalogs[key].name}_{cmbMap_list[j].name}/stackedmap_{filterType}_{est}.npz", x_grid=x, y_grid=y, stackedMap=stackedMap)

        cp=ax.pcolormesh(x*180.*60./np.pi, y*180.*60./np.pi, stackedMap, linewidth=0, rasterized=True)
        # cp.set_cmap('viridis')
        ax.set_xlim(np.min(x)*180.*60./np.pi, np.max(x)*180.*60./np.pi)
        ax.set_ylim(np.min(y)*180.*60./np.pi, np.max(y)*180.*60./np.pi)
        ax.set_xlabel('$x$ [arcmin]')
        ax.set_ylabel('$y$ [arcmin]')

    ax.set_title(key)
    
# ax.legend(fontsize=10, labelspacing=0.1)
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
fig.savefig(path, dpi=100) # bbox_inches='tight')

print('Done!!!')
quit()

# Parameters
path = "./figures/thumbstack/comparison_plot_srcfreemask.png"

# Plotting 


# Plot
fig, subplots = plt.subplots(2, 2, figsize=(8,8), sharex='col', sharey='row')
subplots = subplots.ravel()

for i, key in enumerate(list(catalogKeys)):
    ax = subplots[i]
    for j in range(len(ts_list)):
        tsj = ts_list[j][i]
    
        ax.errorbar(tsj.RApArcmin, factor * tsj.stackedProfile[filterType+"_"+est], 
                    factor * tsj.sStackedProfile[filterType+"_"+est], 
                    label=CMB_namepublic[j], lw=2, ls='-.', capsize=6)


    # ax.plot(ts2.RApArcmin, factor * ts2.stackedProfile[filterType+"_"+est+"_theory_tsz"], ls='--', 
    #     label="theory tsz")
    ax.set_title(key)
    ax.grid()
    
ax.legend(fontsize=10, labelspacing=0.1)
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
fig.savefig(path, dpi=100) # bbox_inches='tight')

fig, subplots = plt.subplots(2, 2, figsize=(8,8), sharex='col', sharey='row')
fig = plt.figure(figsize=(8,8))
# subplots = subplots.ravel()
for i, key in enumerate(list(catalogKeys)):
    # ax = subplots[i]
    ts0 = ts_list[0][i]
    ts1 = ts_list[1][i]
    plt.plot(ts0.RApArcmin, ts0.sStackedProfile[filterType+"_"+est]/ts1.sStackedProfile[filterType+"_"+est], label=key)
    # ax.set_title(key)
plt.grid()
plt.title('no CIB/fiducial')
plt.legend()
plt.savefig("./figures/thumbstack/ratio.png")
           
plt.show()

print('Done!!!')

# Now we have the multiple subplots code (similiar to above) (plotted in the notebook)

fig, subplots = plt.subplots(2, 2, figsize=(8,8), sharex=True, sharey=True)
subplots = subplots.ravel()

for i, key in enumerate(list(catalogKeys)):
    ts5 = ts5_list[i]
    ts6 = ts6_list[i]
    ax = subplots[i]
    
    ax.errorbar(ts5.RApArcmin, factor * ts5.stackedProfile[filterType+"_"+est], 
                factor * ts5.sStackedProfile[filterType+"_"+est], 
                label='ACT DR5 (150GHz $\Delta T$ map)', c='r', lw=2, ls='-.', capsize=6)

    ax.errorbar(ts6.RApArcmin, factor * ts6.stackedProfile[filterType+"_"+est], 
                factor * ts6.sStackedProfile[filterType+"_"+est], 
                label='ACT DR6 (combined y parameter map)', c='b', lw=2, capsize=6)

    ax.plot(ts6.RApArcmin, factor * ts6.stackedProfile[filterType+"_"+est+"_theory_tsz"], ls='--', 
        label="theory tsz, "+' '+ts6.name.replace('_',' '))
    
ax.legend(fontsize=10, labelspacing=0.1)
print('Done!!!')
