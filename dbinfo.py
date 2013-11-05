import numpy as np
import os

package_dir, package_filename = os.path.split(__file__)

isotropic1024coarse = {'name'   : 'isotropic1024coarse',
                       'nx'     : 1024,
                       'ny'     : 1024,
                       'nz'     : 1024,
                       'lx'     : 2*np.pi,
                       'ly'     : 2*np.pi,
                       'lz'     : 2*np.pi,
                       'dx'     : np.pi/512,
                       'dy'     : np.pi/512,
                       'dz'     : np.pi/512,
                       'xnodes' : (np.pi/512)*np.array(range(1024), dtype = np.float32),
                       'ynodes' : (np.pi/512)*np.array(range(1024), dtype = np.float32),
                       'znodes' : (np.pi/512)*np.array(range(1024), dtype = np.float32)}

channel = {'name'   : 'channel',
           'xnodes' : np.load(os.path.join(package_dir, 'data/channel_xgrid.npy')),
           'ynodes' : np.load(os.path.join(package_dir, 'data/channel_ygrid.npy')),
           'znodes' : np.load(os.path.join(package_dir, 'data/channel_zgrid.npy'))}

for coord in ['x', 'y', 'z']:
    channel['n' + coord] = channel[coord + 'nodes'].shape[0]
    channel['d' + coord] = np.zeros(channel['n' + coord], dtype = channel[coord + 'nodes'].dtype)
    channel['d' + coord][:channel['n' + coord] - 1] = (channel[coord + 'nodes'][1:]
                                                     - channel[coord + 'nodes'][:channel['n' + coord] - 1])
    channel['d' + coord][channel['n' + coord] - 1] = channel['d' + coord][0]
    channel['l' + coord] = np.sum(channel['d' + coord])
