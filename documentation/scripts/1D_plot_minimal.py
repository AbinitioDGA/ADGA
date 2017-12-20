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
