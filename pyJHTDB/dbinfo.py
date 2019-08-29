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

interpolation_code = {
        # spatial interpolation
        'NoSInt'                 :   0,
        'NoSpatialInterpolation' :   0,
        'Lag4'                   :   4,
        'Lag6'                   :   6,
        'Lag8'                   :   8,
        'Lagrangian4thOrder'     :   4,
        'Lagrangian6thOrder'     :   6,
        'Lagrangian8thOrder'     :   8,
        'FD4NoInt'               :  40,
        'FD6NoInt'               :  60,
        'FD8NoInt'               :  80,
        'FD4Lag4'                :  44,
        'M1Q4'                   : 104,
        'M1Q6'                   : 106,
        'M1Q8'                   : 108,
        'M1Q10'                  : 110,
        'M1Q12'                  : 112,
        'M1Q14'                  : 114,
        'M2Q4'                   : 204,
        'M2Q6'                   : 206,
        'M2Q8'                   : 208,
        'M2Q10'                  : 210,
        'M2Q12'                  : 212,
        'M2Q14'                  : 214,
        'M3Q4'                   : 304,
        'M3Q6'                   : 306,
        'M3Q8'                   : 308,
        'M3Q10'                  : 310,
        'M3Q12'                  : 312,
        'M3Q14'                  : 314,
        'M4Q4'                   : 404,
        'M4Q6'                   : 406,
        'M4Q8'                   : 408,
        'M4Q10'                  : 410,
        'M4Q12'                  : 412,
        'M4Q14'                  : 414,
        # temporal interpolation
        'NoTInt'                 :   0,
        'PCHIPInt'               :   1,
        'NoTemporalInterpolation':   0,
        'PCHIPInterpolation'     :   1,
        }

package_dir, package_filename = os.path.split(__file__)

isotropic1024coarse = {'name'   : 'isotropic1024coarse'}

for coord in ['x', 'y', 'z']:
    isotropic1024coarse[coord + 'nodes'] = (np.pi/512)*np.array(list(range(1024)), dtype = np.float32)
    isotropic1024coarse['n' + coord] = 1024
    isotropic1024coarse['l' + coord] = 2*np.pi
    isotropic1024coarse['d' + coord] = np.pi/512
    isotropic1024coarse[coord + 'periodic'] = True
    isotropic1024coarse[coord + 'uniform'] = True
isotropic1024coarse['time'] = np.array(list(range(5028)), dtype = np.float32) * 10.056 / 5028
isotropic1024coarse['diss'] = 0.0928
isotropic1024coarse['nu']   = 0.000185


isotropic1024fine = {'name'   : 'isotropic1024fine'}
isotropic1024fine['name'] = 'isotropic1024fine'

for coord in ['x', 'y', 'z']:
    isotropic1024fine[coord + 'nodes'] = (np.pi/512)*np.array(list(range(1024)), dtype = np.float32)
    isotropic1024fine['n' + coord] = 1024
    isotropic1024fine['l' + coord] = 2*np.pi
    isotropic1024fine['d' + coord] = np.pi/512
    isotropic1024fine[coord + 'periodic'] = True
    isotropic1024fine[coord + 'uniform'] = True
isotropic1024fine['time'] = np.array(list(range(100)), dtype = np.float32) * 0.0198/99
isotropic1024fine['diss'] = 0.0928
isotropic1024fine['nu']   = 0.000185

isotropic4096 = {'name'   : 'isotropic4096'}
isotropic4096['name'] = 'isotropic4096'

for coord in ['x', 'y', 'z']:
    isotropic4096[coord + 'nodes'] = (np.pi/2048)*np.array(list(range(4096)), dtype = np.float32)
    isotropic4096['n' + coord] = 4096
    isotropic4096['l' + coord] = 2*np.pi
    isotropic4096['d' + coord] = np.pi/2048
    isotropic4096[coord + 'periodic'] = True
    isotropic4096[coord + 'uniform'] = True
isotropic4096['diss'] = 1.4144
isotropic4096['nu']   = 1.732e-4
isotropic4096['time']   = [0]


rotstrat4096 = {'name'   : 'rotstrat4096'}
rotstrat4096['name'] = 'rotstrat4096'

for coord in ['x', 'y', 'z']:
    rotstrat4096[coord + 'nodes'] = (np.pi/2048)*np.array(list(range(4096)), dtype = np.float32)
    rotstrat4096['n' + coord] = 4096
    rotstrat4096['l' + coord] = 2*np.pi
    rotstrat4096['d' + coord] = np.pi/2048
    rotstrat4096[coord + 'periodic'] = True
    rotstrat4096[coord + 'uniform'] = True
rotstrat4096['diss'] = 0.0123
rotstrat4096['nu']   = 4e-5
rotstrat4096['time']   = range(0,5)

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
channel['time'] = np.array(list(range(4000)), dtype = np.float32) * 0.0065
channel['nu'] = 5e-5

channel5200 = {
        'name'   : 'channel5200',
        'xnodes' : np.array(range(10240))*(8*np.pi/10240),
        'ynodes' : np.load(os.path.join(package_dir, 'data/channel5200_ygrid.npy')),
        'znodes' : np.array(range(7680))*(3*np.pi/7680),
        'lx'     : 8*np.pi,
        'ly'     : 2.,
        'lz'     : 3*np.pi,}

for coord in ['x', 'z']:
    channel5200['n' + coord] = channel5200[coord + 'nodes'].shape[0]
    channel5200[coord + 'periodic'] = True
    channel5200[coord + 'uniform'] = True
    channel5200['d' + coord] = channel5200['l' + coord] / channel5200['n' + coord]

channel5200['ny'] = channel5200['ynodes'].shape[0]
channel5200['dy'] = channel5200['ynodes'][1:] - channel5200['ynodes'][:channel5200['ynodes'].shape[0]-1]
channel5200['dy'] = np.append(channel5200['dy'], [channel5200['dy'][0]])
channel5200['yperiodic'] = False
channel5200['yuniform'] = False
channel5200['time'] = np.array(list(range(10)), dtype = np.float32) # FIXME
channel5200['nu'] = 8e-6

transition_bl = {
        'name'   : 'transition_bl',
        'xnodes' : np.array(range(3320))*0.292210466,
        'ynodes' : np.load(os.path.join(package_dir, 'data/transition_bl_ygrid.npy')),
        'znodes' : np.array(range(2048))*0.117244748,
        'lx'     : 969.8465,
        'ly'     : 2.,
        'lz'     : 240.,}

for coord in ['z']:
    transition_bl['n' + coord] = transition_bl[coord + 'nodes'].shape[0]
    transition_bl[coord + 'periodic'] = True
    transition_bl[coord + 'uniform'] = True
    transition_bl['d' + coord] = channel5200['l' + coord] / channel5200['n' + coord]

transition_bl['nx'] = channel5200['xnodes'].shape[0]
transition_bl['yperiodic'] = False
transition_bl['yuniform'] = True

transition_bl['ny'] = channel5200['ynodes'].shape[0]
transition_bl['dy'] = channel5200['ynodes'][1:] - channel5200['ynodes'][:channel5200['ynodes'].shape[0]-1]
transition_bl['dy'] = np.append(channel5200['dy'], [channel5200['dy'][0]])
transition_bl['yperiodic'] = False
transition_bl['yuniform'] = False
transition_bl['time'] = np.array(list(range(4701)), dtype = np.float32) * 0.25
transition_bl['nu'] = 0.00125

def generate_temp_dbinfo(
        field):
    info = {'name': 'temp',
            'xnodes' : np.arange(0, field.shape[2], 1).astype(np.float32),
            'ynodes' : np.arange(0, field.shape[1], 1).astype(np.float32),
            'znodes' : np.arange(0, field.shape[0], 1).astype(np.float32),
            'xperiodic' : True,
            'xuniform'  : True,
            'yperiodic' : True,
            'yuniform'  : True,
            'zperiodic' : True,
            'zuniform'  : True,
            'dx' : 1.,
            'dy' : 1.,
            'dz' : 1.,
            'lx' : 1.*field.shape[0],
            'ly' : 1.*field.shape[1],
            'lz' : 1.*field.shape[2]}
    return info

