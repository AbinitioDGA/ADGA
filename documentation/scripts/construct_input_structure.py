#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py

def read_1p_data(iw,siw,giw,dc,mu,beta):
    # read in your data according to your format.
    # save them in iw,siw,giw,dc,mu,beta as necessary for ADGA
    pass

def read_2p_data(group, g2):
    # read in one spin-band group of your two-particle data
    pass


read_1p_data(iw,siw,giw,dc,mu,beta)

f = h5py.File('1p-data.hdf5','w')
f['.axes/iw']=iw
f.create_group('.config')
f['.config'].attrs['general.beta']=beta
f['dmft-001/mu/value']=mu
f['dmft-001/ineq-001/dc/value']=dc
f['dmft-001/ineq-001/siw/value']=siw
f['dmft-001/ineq-001/giw/value']=giw
f.close()


g = h5py.File('2p-data.hdf5','w')

# groups holds all spin-band combinations of
# the two particle Greens function
# as strings, e.g. '00001'
for group in groups: 
    read_2p_data(group, g2)
    g['worm-001/ineq-001/g4iw-worm/'+group+'/value'] = g2
g.close()


n4iwf=g2.shape[0]//2
n4iwb=g2.shape[-1]//2

g['.axes/iwf-g4'] = np.pi/beta * (2 * np.arange(-n4iwf,n4iwf) +1)
g['.axes/iwb-g4'] = np.pi/beta * 2 * np.arange(-n4iwb,n4iwb+1)

g.close()
