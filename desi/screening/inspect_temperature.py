import numpy as np

fn = "/pscratch/sd/b/boryanah/ACTxDESI/DESI/DESI_pz1/extended_catalog_allfoot_perbin.txt"
data = np.genfromtxt(fn)
vZ = data[:, 14] # T_large-scales
