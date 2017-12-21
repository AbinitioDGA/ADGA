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
# Copyright (C) <2017, 2018> 
# <Anna Galler*, Patrick Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
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
