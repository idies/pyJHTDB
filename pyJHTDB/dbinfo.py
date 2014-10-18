########################################################################
#
#  Copyright 2014 Johns Hopkins University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Contact: turbulence@pha.jhu.edu
# Website: http://turbulence.pha.jhu.edu/
#
########################################################################

import numpy as np
import os

package_dir, package_filename = os.path.split(__file__)

isotropic1024coarse = {'name'   : 'isotropic1024coarse'}

for coord in ['x', 'y', 'z']:
    isotropic1024coarse[coord + 'nodes'] = (np.pi/512)*np.array(list(range(1024)), dtype = np.float32)
    isotropic1024coarse['n' + coord] = 1024
    isotropic1024coarse['l' + coord] = 2*np.pi
    isotropic1024coarse['d' + coord] = np.pi/512
    isotropic1024coarse[coord + 'periodic'] = True
    isotropic1024coarse[coord + 'uniform'] = True
isotropic1024coarse['time'] = np.array(list(range(1024)), dtype = np.float32) * 2.048 / 1024

mhd1024 = {}
for key in list(isotropic1024coarse.keys()):
    mhd1024[key] = isotropic1024coarse[key]
mhd1024['name'] = 'mhd1024'

for coord in ['x', 'y', 'z']:
    mhd1024[coord + 'nodes'] = (np.pi/512)*np.array(list(range(1024)), dtype = np.float32)
    mhd1024['n' + coord] = 1024
    mhd1024['l' + coord] = 2*np.pi
    mhd1024['d' + coord] = np.pi/512
    mhd1024[coord + 'periodic'] = True
    mhd1024[coord + 'uniform'] = True
mhd1024['time'] = np.array(list(range(1024)), dtype = np.float32) * 2.56 / 1024
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
