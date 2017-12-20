#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

f = h5py.File('adga-1234567890.hdf5','r')

dset = f['selfenergy/nonloc/dga'][()]
xaxis = np.linspace(-29*np.pi/10.0,29*np.pi/10.0,60)

fig1 = plt.figure()

ax1 = fig1.add_subplot(211)
ax1.set_ylabel(r'$\mathrm{Re}[\Sigma(\mathbf{k},i\nu)]$')
plt.plot(xaxis,dset[0,0,0,0,0,:].real, label=r'$\Gamma$')
plt.plot(xaxis,dset[0,0,0,10,0,:].real, label=r'$X$')
plt.plot(xaxis,dset[0,0,10,10,0,:].real, label=r'$M$')
plt.plot(xaxis,dset[0,0,10,10,10,:].real, label=r'$R$')
plt.legend(loc='lower right')

ax1 = fig1.add_subplot(212)
ax1.set_xlabel(r'$i\nu$')
ax1.set_ylabel(r'$\mathrm{Im}[\Sigma(\mathbf{k},i\nu)]$')
plt.plot(xaxis,dset[0,0,0,0,0,:].imag, label=r'$\Gamma$')
plt.plot(xaxis,dset[0,0,0,10,0,:].imag, label=r'$X$')
plt.plot(xaxis,dset[0,0,10,10,0,:].imag, label=r'$M$')
plt.plot(xaxis,dset[0,0,10,10,10,:].imag, label=r'$R$')
plt.legend(loc='upper right')

plt.savefig('qgrid.eps',format='eps')
