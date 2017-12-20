#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

f = h5py.File('adga-1234567890.hdf5','r')
dset = f['selfenergy/nonloc/dga'][()]

fig = plt.figure()

a = np.zeros((22,22), dtype=np.complex128)
a[:20,:20] = dset[0,0,...,0,30]
a[20,:] = a[0,:]
a[:,20] = a[:,0]

ax1 = fig.add_subplot(121)
p1 = ax1.pcolormesh(np.arange(22), np.arange(22),a[...].real, rasterized=True)
ax1.set_xlim(0,21)
ax1.set_ylim(0,21)
plt.xticks([0,10.5,21],[r'$0$', r'$\pi$', r'$2\pi$'])
plt.yticks([0,10.5,21],[r'$0$', r'$\pi$', r'$2\pi$'])
ax1.set_title(r'$\mathrm{Re}[\Sigma(\mathbf{k},i\nu_0)]$')
ax1.set_aspect('equal')
plt.colorbar(p1,shrink=0.47,format='%.3f')

ax2 = fig.add_subplot(122)
ax2.set_title(r'$\mathrm{Im}[\Sigma(\mathbf{k},i\nu_0)]$')
p2 = ax2.pcolormesh(np.arange(22), np.arange(22), a[...].imag, rasterized=True)
ax2.set_xlim(0,21)
ax2.set_ylim(0,21)
plt.xticks([0,10.5,21],[r'$0$', r'$\pi$', r'$2\pi$'])
plt.yticks([0,10.5,21],[r'$0$', r'$\pi$', r'$2\pi$'])
ax2.set_aspect('equal')
plt.colorbar(p2,shrink=0.47,format='%.3f')

plt.savefig('qgrid_plane.eps',format='eps', Rasterized=True)

plt.show()
