#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

f = h5py.File('adga-12345.hdf5','r')

# extract the fermionic vertex box size
iwf = f['input/iwfmax_small'][()]

dset = f['selfenergy/nonloc/dga'][()]

# plot the first band at the first fermionic frequency
# in the kz = 0 plane
f = plt.figure()

f.add_subplot(211) # 2x1 subplots
plt.pcolormesh(dset[0,0,:,:,0,iwf].imag)

f.add_subplot(212)
plt.pcolormesh(dset[0,0,:,:,0,iwf].real)

plt.show
