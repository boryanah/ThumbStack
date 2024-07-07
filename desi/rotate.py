#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, gc
import importlib

import numpy as np
import matplotlib.pyplot as plt

from pixell import enmap, enplot, reproject, curvedsky as cs, utils

#sys.path.append("../../")
import rotfuncs


# In[6]:


def eshow(x, **kwargs):
    ''' Define a function to help us plot the maps neatly '''
    plots = enplot.get_plots(x, **kwargs)
    enplot.show(plots, method = "python")
    #enplot.write(fn, plots)


# In[5]:


# define
pathMap = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits' # DR6 ILC
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended.fits' # Already smoothed, I am 90% sure

# read
ACT_map = enmap.read_map(pathMap)
#ACT_mask = enmap.read_map(pathMask)


# In[7]:


eshow(ACT_map)


# In[ ]:


plot = enplot.plot(ACT_map)
enplot.show(plot)


# In[9]:


enmap.box(ACT_map.shape,ACT_map.wcs)/utils.degree


# In[316]:


del ACT_map; gc.collect()


# In[15]:


dec_min = -15*utils.degree
dec_max = -5*utils.degree
ra_min = 20.*utils.degree
ra_max = 30.*utils.degree
box = np.array([[dec_min,ra_min],[dec_max,ra_max]]) # in radians
omap = ACT_map.submap(box)
#omap = enmap.submap(imap,box) # an alternative way
enmap.write_fits("test.fits", omap)


# In[114]:


omap = enmap.read_fits("test.fits")


# In[115]:


omap.shape, omap.wcs


# In[116]:


plt.imshow(omap)


# In[117]:


enmap.box(omap.shape, omap.wcs)/utils.degree


# In[118]:


opos = omap.posmap()


# In[119]:


opos[0].min()*180./np.pi, opos[0].max()*180./np.pi, opos[1].min()*180./np.pi, opos[1].max()*180./np.pi


# In[120]:


ra = -1.#23.
dec = -1. #-7.
sourcecoord = np.array([ra, dec])*utils.degree 


# In[121]:


importlib.reload(rotfuncs)


# In[122]:


ipos = rotfuncs.recenter(opos[::-1], [0, 0, sourcecoord[0], sourcecoord[1]])[::-1]


# In[123]:


ipos[0].min()*180./np.pi, ipos[0].max()*180./np.pi, ipos[1].min()*180./np.pi, ipos[1].max()*180./np.pi


# In[124]:


#shape, wcs = enmap.geometry(np.array([[-5., -5.],[5., 5.]])*utils.degree, res=0.5*utils.arcmin)
shape, wcs = enmap.geometry(np.array([[90-5., -5.],[90+5., 5.]])*utils.degree, res=0.5*utils.arcmin) # TESTING dec ra note that dec cannot exceed 90
stampMap = enmap.zeros(shape, wcs)


# In[125]:


#stampMap = omap.copy() # cause our true map and stamp map are the same size but typically not
stampMap[:, :] = omap.at(ipos, prefilter=True, mask_nan=False, order=1)


# In[136]:


plt.imshow(stampMap)


# In[127]:


opos_new = stampMap.posmap()
print(opos_new[0].min()*180./np.pi, opos_new[0].max()*180./np.pi, opos_new[1].min()*180./np.pi, opos_new[1].max()*180./np.pi)
ipos_new = rotfuncs.recenter(opos_new[::-1], [0.05, -0.005, 0.005, -0.005])[::-1]

stampMap_new = enmap.zeros(shape, wcs)
stampMap_new[:, :] = stampMap.at(ipos_new, prefilter=True, mask_nan=False, order=1)


# In[128]:


print(ipos_new[0].min()*180./np.pi, ipos_new[0].max()*180./np.pi, ipos_new[1].min()*180./np.pi, ipos_new[1].max()*180./np.pi)


# In[129]:


opos_newnew = stampMap_new.posmap()
print(opos_newnew[0].min()*180./np.pi, opos_newnew[0].max()*180./np.pi, opos_newnew[1].min()*180./np.pi, opos_newnew[1].max()*180./np.pi)


# In[130]:


plt.imshow(stampMap_new)


# In[131]:


def moveaxis(a, o, n):
    if o < 0: o = o+a.ndim
    if n < 0: n = n+a.ndim
    if n <= o: return np.rollaxis(a, o, n)
    else: return np.rollaxis(a, o, n+1)

def ang2rect(angs, zenith=True, axis=0):
    """Convert a set of angles [{phi,theta},...] to cartesianc
    coordinates [{x,y,z},...]. If zenith is True (the default),
    the theta angle will be taken to go from 0 to pi, and measure
    the angle from the z axis. If zenith is False, then theta
    goes from -pi/2 to pi/2, and measures the angle up from the xy plane."""
    phi, theta = moveaxis(angs, axis, 0)
    ct, st, cp, sp = np.cos(theta), np.sin(theta), np.cos(phi), np.sin(phi)
    if zenith: res = np.array([st*cp,st*sp,ct])
    else:      res = np.array([ct*cp,ct*sp,st])
    return moveaxis(res, 0, axis)

def rect2ang(rect, zenith=True, axis=0):
    """The inverse of ang2rect."""
    x,y,z = moveaxis(rect, axis, 0)
    r     = (x**2+y**2)**0.5
    phi   = np.arctan2(y,x)
    if zenith: theta = np.arctan2(r,z)
    else:      theta = np.arctan2(z,r)
    return moveaxis(np.array([phi,theta]), 0, axis)

def rotmatrix(ang, raxis, axis=0):
    """Construct a 3d rotation matrix representing a rotation of
    ang degrees around the specified rotation axis raxis, which can be "x", "y", "z"
    or 0, 1, 2. If ang is a scalar, the result will be [3,3]. Otherwise,
    it will be ang.shape + (3,3)."""
    ang  = np.asarray(ang)
    raxis = raxis.lower()
    c, s = np.cos(ang), np.sin(ang)
    R = np.zeros(ang.shape + (3,3))
    if   raxis == 0 or raxis == "x": R[...,0,0]=1;R[...,1,1]= c;R[...,1,2]=-s;R[...,2,1]= s;R[...,2,2]=c
    elif raxis == 1 or raxis == "y": R[...,0,0]=c;R[...,0,2]= s;R[...,1,1]= 1;R[...,2,0]=-s;R[...,2,2]=c
    elif raxis == 2 or raxis == "z": R[...,0,0]=c;R[...,0,1]=-s;R[...,1,0]= s;R[...,1,1]= c;R[...,2,2]=1
    else: raise ValueError("Rotation axis %s not recognized" % raxis)
    return moveaxis(R, 0, axis)


# In[140]:


coords = opos_new[::-1]
coords = np.asarray(coords)
co     = coords.reshape(2,-1)
rect   = ang2rect(co, False)

print(rect[0].min(), rect[0].max(), rect[1].min(), rect[1].max(), rect[2].min(), rect[2].max())

M = rotmatrix(np.pi/1000., "z")

rect   = np.einsum("...ij,j...->i...", M, rect)

# TESTING
rect = ang2rect(co, False)

co     = rect2ang(rect, False)
ipos_newnew = co.reshape(coords.shape)



print(rect[0].min(), rect[0].max(), rect[1].min(), rect[1].max(), rect[2].min(), rect[2].max())


# In[141]:


print(ipos_newnew[0].min()*180./np.pi, ipos_newnew[0].max()*180./np.pi, ipos_newnew[1].min()*180./np.pi, ipos_newnew[1].max()*180./np.pi)


# In[142]:


plt.imshow(stampMap.at(ipos_newnew, prefilter=True, mask_nan=False, order=1))


# In[144]:


plt.imshow(stampMap)


# In[150]:


ipos = stampMap.posmap()


# In[153]:


plt.imshow(stampMap.at(ipos, prefilter=True, mask_nan=False, order=1))


# In[157]:


plt.imshow(ipos[0]) # ra
plt.colorbar()


# In[156]:


plt.imshow(ipos[1]) # dec
plt.colorbar()


# In[291]:


ang = np.pi/8#2 #10.
R = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])


# In[292]:


x = np.linspace(0., 10, 100)
y = np.linspace(10., 20, 100)
X, Y = np.meshgrid(x, y)


# In[293]:


X, Y = ipos


# In[294]:


pix = enmap.sky2pix(shape, wcs, ipos)
X, Y = pix


# In[295]:


mean_X = np.mean(X)
mean_Y = np.mean(Y)

X_new = X-mean_X
Y_new = Y-mean_Y

XY = np.array([X_new.flatten(), Y_new.flatten()])

XY_rot = np.dot(R, XY)
X_rot = XY_rot[0].reshape(X.shape) 
Y_rot = XY_rot[1].reshape(Y.shape)

#mean_XY_rot = np.dot(R, np.array([mean_X, mean_Y]))
#X_rot += mean_XY_rot[0]
#Y_rot += mean_XY_rot[1]
X_rot += mean_X
Y_rot += mean_Y


# In[300]:


ipos_rot = np.array([X_rot, Y_rot])


# In[301]:


pix_rot = np.array([X_rot, Y_rot])


# In[307]:


from scipy.interpolate import RectBivariateSpline

fun = RectBivariateSpline(X[:, 0], Y[0, :], stampMap, kx=1, ky=1)


# In[308]:


meow = fun(pix_rot[0], pix_rot[1], grid=False)


# In[309]:


plt.imshow(meow.reshape(1200, 1200))


# In[310]:


plt.imshow(stampMap)


# In[ ]:




