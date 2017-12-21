# This file is part of the Abinitio Dynamical Vertex Approximation (ADGA)
# package. It's an electronic structure code which allows the inclusion of
# non-local correlations beyond DMFT.
#
# The public repository can be found at
# https://github.com/AbinitioDGA/ADGA
#
# The arXiv publication can be found at
# https://arxiv.org/abs/1710.06651
#
# Copyright (C) <2017, 2018> 
# <Anna Galler, Patrick Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
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
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

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