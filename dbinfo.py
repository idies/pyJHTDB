import numpy as np
import os

package_dir, package_filename = os.path.split(__file__)

isotropic1024coarse = {'name'   : 'isotropic1024coarse'}

for coord in ['x', 'y', 'z']:
    isotropic1024coarse[coord + 'nodes'] = (np.pi/512)*np.array(range(1024), dtype = np.float32)
    isotropic1024coarse['n' + coord] = 1024,
    isotropic1024coarse['l' + coord] = 2*np.pi
    isotropic1024coarse['d' + coord] = np.pi/512
    isotropic1024coarse[coord + 'periodic'] = True
    isotropic1024coarse[coord + 'uniform'] = True

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
    channel[coord + 'periodic'] = True
    channel[coord + 'uniform'] = True

channel['yperiodic'] = False
channel['yuniform'] = False
