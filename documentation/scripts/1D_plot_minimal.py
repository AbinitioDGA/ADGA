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

# open the file in the read-format
f = h5py.File('adga-12345.hdf5','r')

# extract the data into a numpy array
dset = f['selfenergy/nonloc/dga'][()]

# this dataset has the form of ndim,ndim,npx,npy,npz,2*iwf
print(dset.shape) # prints the shape of the array

# plot the first band at the gamma-point
plt.plot(dset[0,0,0,0,0,:].imag)

# show the figure now
plt.show()
