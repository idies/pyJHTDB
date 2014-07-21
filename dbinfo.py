import numpy as np
import os

package_dir, package_filename = os.path.split(__file__)

isotropic1024coarse = {'name'   : 'isotropic1024coarse'}

for coord in ['x', 'y', 'z']:
    isotropic1024coarse[coord + 'nodes'] = (np.pi/512)*np.array(range(1024), dtype = np.float32)
    isotropic1024coarse['n' + coord] = 1024
    isotropic1024coarse['l' + coord] = 2*np.pi
    isotropic1024coarse['d' + coord] = np.pi/512
    isotropic1024coarse[coord + 'periodic'] = True
    isotropic1024coarse[coord + 'uniform'] = True

mhd1024 = {}
for key in isotropic1024coarse.keys():
    mhd1024[key] = isotropic1024coarse[key]
mhd1024['name'] = 'mhd1024'

for coord in ['x', 'y', 'z']:
    mhd1024[coord + 'nodes'] = (np.pi/512)*np.array(range(1024), dtype = np.float32)
    mhd1024['n' + coord] = 1024
    mhd1024['l' + coord] = 2*np.pi
    mhd1024['d' + coord] = np.pi/512
    mhd1024[coord + 'periodic'] = True
    mhd1024[coord + 'uniform'] = True
mhd1024['time'] = np.array(range(1024), dtype = np.float32) * 2.56 / 1024
mhd1024['nu']  = 1.1e-4
mhd1024['eta'] = 1.1e-4
mhd1024['diss_u'] = 1.1e-2
mhd1024['diss_b'] = 2.2e-2
mhd1024['dt']     = 2.5e-3
mhd1024['T']      = 2.56

channel = {'name'   : 'channel',
           'xnodes' : np.load(os.path.join(package_dir, 'data/channel_xgrid.npy')),
           'ynodes' : np.load(os.path.join(package_dir, 'data/channel_ygrid.npy')),
           'znodes' : np.load(os.path.join(package_dir, 'data/channel_zgrid.npy')),
           'lx'     : 8*np.pi,
           'ly'     : 2.,
           'lz'     : 3*np.pi,}

for coord in ['x', 'z']:
    channel['n' + coord] = channel[coord + 'nodes'].shape[0]
    channel[coord + 'periodic'] = True
    channel[coord + 'uniform'] = True
    channel['d' + coord] = channel['l' + coord] / channel['n' + coord]

channel['ny'] = 512
channel['dy'] = channel['ynodes'][1:] - channel['ynodes'][:channel['ynodes'].shape[0]-1]
channel['dy'] = np.append(channel['dy'], [channel['dy'][0]])
channel['yperiodic'] = False
channel['yuniform'] = False
