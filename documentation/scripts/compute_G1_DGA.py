#! /usr/bin/env python

# This file is part of the Abinitio Dynamical Vertex Approximation (ADGA)
# package. It is an electronic structure code which allows the inclusion of
# non-local correlations beyond DMFT and the calculation of momentum-dependent
# susceptibilities.
#
# The public repository can be found at
# https://github.com/AbinitioDGA/ADGA
#
# The arXiv publication can be found at
# https://arxiv.org/abs/1710.06651
#
# Copyright (C) <2017-2019>
# <Anna Galler*, Patrick Thunstr\"om, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
# * Corresponding author. E-mail address: galler.anna@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
beta = f['input/beta'][()]

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
