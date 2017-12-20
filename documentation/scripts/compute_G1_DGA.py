#! /usr/bin/env python

from __future__ import print_function,
import numpy as np
import h5py
import matplotlib.pyplot as plt
import scipy.linalg

# open the file in the read-format
f = h5py.File('adga-12345.hdf5','r')

# extract the necessary components
siwk_dga = f['selfenergy/nonloc/dga'][()]
hk = f['input/hk'][()]
dc = f['input/dc'][:,0]
mu = f['input/mu'][()]
iwf = f['input/iwfmax_small'][()]
ndim = hk.shape[0]
nqx = f['input/nqpxyz'][0]
nqy = f['input/nqpxyz'][1]
nqz = f['input/nqpxyz'][2]

# create the matsubara axis
fmats = np.linspace(-(iwf*2-1)*np.pi/beta,(iwf*2-1)*np.pi/beta,2*iwf)

# building DGA Greens function from scratch
# according to [iw + mu - dc - H(k) - Sigma(k,iw)]**(-1)
gdgainv = np.zeros((ndim,ndim,nqx,nqy,nqz,2*iwf),dtype=np.complex128)
gdgainv += -siwk_dga - hk[...,None]
gdgainv[np.arange(ndim),np.arange(ndim),...] += 1j*fmats+mu-dc[np.arange(ndim),None,None,None,None]

gdga = np.empty_like(gdgainv, dtype=np.complex128)

for ikx in xrange(nqx):
    for iky in xrange(nqy):
        for ikz in xrange(nqz):
            for iw in xrange(2*iwf):
                gdga[:,:,ikx,iky,ikz,iw] = scipy.linalg.inv(gdgainv[:,:,ikx,iky,ikz,iw])

# plot the first band at the gamma-point
plt.plot(gdga[0,0,0,0,0,:].imag)

# show the figure now
plt.show()
