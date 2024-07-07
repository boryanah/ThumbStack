import gc

import numpy as np
from scipy import optimize, special

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15

r = 0.3
sr_to_arcmin2 = (180.*60./np.pi)**2.
theta_arcmin = np.linspace(1, 6, 9)
want_other = True
cosmo = Planck15

#pairs = ["_extended_allfoot_perbin_dr6", "_extended_dr10_allfoot_perbin_dr6"]
#pairs = ["_allfoot_perbin_dr6", "_dr10_allfoot_perbin_dr6"]
pairs = ["_extended_allfoot_perbin_sigmaz0.0500_dr6", "_extended_dr10_allfoot_perbin_sigmaz0.0500_dr6"]
#pairs = ["_allfoot_perbin_sigmaz0.0500_dr6", "_dr10_allfoot_perbin_sigmaz0.0500_dr6"]

for pz_bin in range(1, 5):
    print("PZ BIN", pz_bin)

    for j in range(len(pairs)):
        print("file type", pairs[j])
        
        # load delta_T
        filtmap = np.loadtxt(f"/pscratch/sd/b/boryanah/ACTxDESI/output{pairs[j]}/thumbstack/DESI_pz{pz_bin:d}_act_dr6_f90/diskring_filtmap.txt") # uK sr
        filtmask = np.loadtxt(f"/pscratch/sd/b/boryanah/ACTxDESI/output{pairs[j]}/thumbstack/DESI_pz{pz_bin:d}_act_dr6_f90/diskring_filtmask.txt") # uK sr
        overlap = np.loadtxt(f"/pscratch/sd/b/boryanah/ACTxDESI/output{pairs[j]}/thumbstack/DESI_pz{pz_bin:d}_act_dr6_f90/overlap_flag.txt")

        # load catalog
        if "extended" in pairs[j]:
            ending_new = pairs[j].split("_extended")[-1]
            ext_str = "extended_"
        else:
            ending_new = pairs[j]
            ext_str = ""
        cat = np.loadtxt(f"/pscratch/sd/b/boryanah/ACTxDESI/DESI/DESI_pz{pz_bin:d}/{ext_str}catalog{ending_new.split('_dr6')[0]}.txt")

        # select overlapping objects that are completely unmasked
        choice = overlap == 1.
        print("all objects", len(choice))
        print("overlap = 1", np.sum(choice), np.sum(choice)/len(choice))
        choice &= (np.abs(filtmask[:, -1]) < 1.)
        print("overlap = 1 and filtmask == 1.", np.sum(choice), np.sum(choice)/len(choice))

        # adding outlier mask
        nObj = np.sum(choice)
        f = lambda nSigmas: nObj * special.erfc(nSigmas / np.sqrt(2.)) - special.erfc(5. / np.sqrt(2.))
        nSigmasCut = optimize.brentq(f , 0., 1.e2)
        sigmas = np.std(filtmap[choice, :], axis=0) # sigmas has shape  nRAp
        # shape is (nObj, nRAp)
        choice *= np.prod((np.abs(filtmap[:, :]) <= nSigmasCut * sigmas[np.newaxis, :]), axis=1).astype(bool)

        # apply all masks
        cat = cat[choice]
        filtmap = filtmap[choice]
        print("overlap = 1 and filtmask == 1. and outlier", np.sum(choice), np.sum(choice)/len(choice))

        # read out RA, DEC and velocity
        RA = cat[:, 0] # deg
        DEC = cat[:, 1] # deg
        Z = cat[:, 2]
        vr = cat[:, 15] / 3.e5 # km/s # unitless
        del cat
        gc.collect()

        if j == 0:
            # load other catalog
            if "extended" in pairs[1]:
                ending_new = pairs[1].split("_extended")[-1]
                ext_str = "extended_"
            else:
                ending_new = pairs[1]
                ext_str = ""
            cat_other = np.loadtxt(f"/pscratch/sd/b/boryanah/ACTxDESI/DESI/DESI_pz{pz_bin:d}/{ext_str}catalog{ending_new.split('_dr6')[0]}.txt")

            # read out RA, DEC and velocity
            RA_other = cat_other[:, 0] # deg
            DEC_other = cat_other[:, 1] # deg
            Z_other = cat_other[:, 2]
            vr_other = cat_other[:, 15] / 3.e5 # km/s # unitless
            del cat_other
            gc.collect()

            # match by RA, DEC
            c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, distance=cosmo.comoving_distance(Z)/(1.+Z))
            c_other = SkyCoord(ra=RA_other*u.degree, dec=DEC_other*u.degree, distance=cosmo.comoving_distance(Z_other)/(1.+Z_other))

            idx, d2d, d3d = c.match_to_catalog_sky(c_other)
            #max_sep = 1.0 * u.arcsec # worse?
            #max_sep = 10.0 * u.arcmin # better for pz bin 1 main
            #max_sep = 100.0 * u.arcmin # better for pz bin 1 main
            #idx[d2d > max_sep] = -1 # I think spinning is better
            max_sep = 10.0 * u.arcmin # better for pz bin 1 main
            idx[(d2d > max_sep) & (d3d > 100*u.Mpc)] = -1
            #print("debug", len(idx), len(RA)) # because no overlap cuts, etc. applied
            #print("debug", np.min(idx), np.max(idx), len(RA_other)) # possibly finding closest but not always right?
            print("idx is -1 out of (frac no match)", np.sum(idx == -1), len(idx), np.sum(idx == -1)/len(idx))

            # substitute values of the velocity
            vr_match = vr.copy()
            vr_match[idx > -1] = vr_other[idx[idx > -1]]

            # subtract mean
            vr_match -= np.mean(vr_match)
            vr_match_rms = np.std(vr_match)

            # compute the delta_T times v_r products
            prod_match = -(filtmap * vr_match[:, None]) / vr_match_rms
            prod_match /= r
            prod_match *= sr_to_arcmin2

            # compute profile and covariance matrix
            prof_match = np.mean(prod_match, axis=0)
            cov_match = np.cov(prod_match.T) / len(vr_match)
            print("MATCH: profile and error", prof_match, np.sqrt(np.diag(cov_match)))
            print("MATCH: chi2 null", np.dot(np.dot(prof_match, np.linalg.inv(cov_match)), prof_match))
            
        # subtract mean
        vr -= np.mean(vr)
        vr_rms = np.std(vr)

        # compute the delta_T times v_r products
        prod = -(filtmap * vr[:, None]) / vr_rms
        prod /= r
        prod *= sr_to_arcmin2

        # compute profile and covariance matrix
        prof = np.mean(prod, axis=0)
        cov = np.cov(prod.T) / len(vr)
        print(f"{j+1:d} MAIN: profile and error", prof, np.sqrt(np.diag(cov)))
        print(f"{j+1:d} MAIN: chi2 null", np.dot(np.dot(prof, np.linalg.inv(cov)), prof))

        if j == 0:
            prof_first = prof.copy()
            cov_first = cov.copy()
            prof_other = np.empty(0)
            cov_other = np.empty(0)
            
        if j == 1:
            prof_other = prof.copy()
            cov_other = cov.copy()
        
        # if you don't want the second one, then skip
        if not want_other:
            continue
    print("------------")
    np.savez(f"prof_{pairs[0]}_match_bin{pz_bin:d}.npz", prof_other=prof_other, cov_other=cov_other, prof_match=prof_match, cov_match=cov_match, prof=prof_first, cov=cov_first)
