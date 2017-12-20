#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

f = h5py.File('adga-1234567890.hdf5','r')

dset_dens = f['susceptibility/nonloc/dens'][()]
dset_magn = f['susceptibility/nonloc/magn'][()]

fig = plt.figure()

ax1 = fig.add_subplot(211)
plt.plot(np.sum(dset_dens[:,:,:,0].real, axis=(0,1)))
ax1.set_ylabel(r'$\chi_D(\mathbf{q},\omega=0)$')
plt.xticks([0,10,20,30,40],[r'$\Gamma$','$X$','$M$','$R$',r'$\Gamma$'])

ax2 = fig.add_subplot(212)
plt.plot(np.sum(dset_magn[:,:,:,0].real, axis=(0,1)))
ax2.set_ylabel(r'$\chi_M(\mathbf{q},\omega=0)$')
plt.xticks([0,10,20,30,40],[r'$\Gamma$','$X$','$M$','$R$',r'$\Gamma$'])

plt.savefig('qpath.eps',format='eps')

plt.show()
