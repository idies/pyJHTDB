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
import pickle
import gzip

import pyJHTDB
import pyJHTDB.generic_splines as gs

from scipy.ndimage.filters import correlate1d

class spline_interpolator:
    def __init__(
            self,
            info,
            n = 1,
            m = 1,
            compute_fast_beta = False):
        self.n = n
        self.m = m
        self.info = info
        if os.path.exists(info['name'] + '_spline_interpolator_n{0}_m{1}.pickle.gz'.format(n, m)):
            self.spline = pickle.load(gzip.open(info['name'] + '_spline_interpolator_n{0}_m{1}.pickle.gz'.format(n, m)))
        else:
            func = []
            for coord in ['x', 'y', 'z']:
                if info[coord + 'uniform'] and info[coord + 'periodic']:
                    func.append(gs.generic_spline_1D(
                            info[coord + 'nodes'][:2],
                            max_deriv = m,
                            neighbours = n,
                            period = info['l' + coord]))
                else:
                    func.append(gs.generic_spline_1D(
                            info[coord + 'nodes'],
                            max_deriv = m,
                            neighbours = n))
                func[-1].compute_derivs()
                func[-1].compute_beta()
            self.spline = {'x': func[0],
                           'y': func[1],
                           'z': func[2]}
            pickle.dump(self.spline,
                        gzip.open(info['name'] + '_spline_interpolator_n{0}_m{1}.pickle.gz'.format(n, m),
                             'wb'))
        # either channel or periodic cube, so it's cheap to compute fast betas for x and z
        self.spline['x'].compute_fast_beta()
        self.spline['z'].compute_fast_beta()
        self.bx = self.spline['x'].fast_beta
        self.bz = self.spline['z'].fast_beta
        if compute_fast_beta:
            self.spline['y'].compute_fast_beta()
            self.by = self.spline['y'].fast_beta
        else:
            self.by = self.spline['y'].beta
        return None
    def __call__(
            self,
            time = 0.,
            points = np.zeros((1,3)),
            dorder = [(0, 0, 0)],
            lTDB = None,
            getFunction = 'getVelocityAndPressureSoap'):
        if (not len(points.shape) == 2):
            return None
        field_points = np.zeros((points.shape[0], 2*self.n+2, 2*self.n+2, 2*self.n+2, 3), dtype = np.float32)
        xgrid = np.floor(points[:, 0] / self.info['dx']).astype(np.int) % self.info['nx']
        xfrac = (points[:, 0] - self.info['xnodes'][xgrid])/self.info['dx']
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 0] = (self.info['xnodes'][xgrid[p]]
                    + np.array(list(range(-self.n, self.n+2)))*self.info['dx'])[np.newaxis, np.newaxis, :]
        if self.info['yperiodic']:
            ygrid = np.floor(points[:, 1] / self.info['dy']).astype(np.int) % self.info['ny']
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy']
            for p in range(points.shape[0]):
                field_points[p, :, :, :, 1] = (self.info['ynodes'][ygrid[p]]
                        + np.array(list(range(-self.n, self.n+2)))*self.info['dy'])[np.newaxis, :, np.newaxis]
        else:
            ygrid = np.searchsorted(self.info['ynodes'], points[:, 1]).astype(np.int) - 1
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy'][ygrid]
            for p in range(points.shape[0]):
                if ygrid[p] < 0 or ygrid[p] > self.info['ny'] - 1:
                    return None
                elif ygrid[p] < self.n:
                    field_points[p, :, :2*self.n+1, :, 1] = self.info['ynodes'][np.newaxis, :2*self.n+1, np.newaxis]
                elif ygrid[p] >= self.info['ny'] - self.n - 1:
                    field_points[p, :, :2*self.n+1, :, 1] = self.info['ynodes'][np.newaxis, self.info['ny'] - (2*self.n+1):, np.newaxis]
                else:
                    field_points[p, :, :, :, 1] = self.info['ynodes'][np.newaxis, ygrid[p]-self.n:ygrid[p]+self.n+2, np.newaxis]
        zgrid = np.floor(points[:, 2] / self.info['dz']).astype(np.int) % self.info['nz']
        zfrac = (points[:, 2] - self.info['znodes'][zgrid])/self.info['dz']
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 2] = (self.info['znodes'][zgrid[p]]
                    + np.array(list(range(-self.n, self.n+2)))*self.info['dz'])[:, np.newaxis, np.newaxis]
        print('computed points where field is needed, now getting values from DB')
        ## I could in principle call getRaw[...] for each point,
        ## but that would mean a lot of calls to the DB,
        ## and we should avoid that due to the latency.
        field_values = lTDB.getData(
                time, field_points,
                sinterp = 0, tinterp = 0,
                data_set = self.info['name'],
                getFunction = getFunction)
        print('got values from DB, now interpolating')
        result = np.zeros((len(dorder), points.shape[0], field_values.shape[-1]), dtype = np.float32)
        bxi = 0
        if self.info['yperiodic']:
            ygrid[:] = 0
        bzi = 0
        if self.info['yperiodic']:
            xb = np.zeros((len(dorder), points.shape[0], len(self.bx[0][dorder[0][0]])), field_values.dtype)
            yb = np.zeros((len(dorder), points.shape[0], len(self.bx[0][dorder[0][1]])), field_values.dtype)
            zb = np.zeros((len(dorder), points.shape[0], len(self.bx[0][dorder[0][2]])), field_values.dtype)
            for p in range(points.shape[0]):
                for o in range(len(dorder)):
                    xb[o, p] = np.array([self.bx[     bxi][dorder[o][0]][k](xfrac[p])
                                         for k in range(len(self.bx[     bxi][dorder[o][0]]))]).astype(field_values.dtype)
                    yb[o, p] = np.array([self.by[ygrid[p]][dorder[o][1]][k](yfrac[p])
                                         for k in range(len(self.by[ygrid[p]][dorder[o][1]]))]).astype(field_values.dtype)
                    zb[o, p] = np.array([self.bz[     bzi][dorder[o][2]][k](zfrac[p])
                                         for k in range(len(self.bz[     bzi][dorder[o][2]]))]).astype(field_values.dtype)
            result = np.einsum('opkjil,opi,opj,opk->opl', field_values[None, :], xb, yb, zb)
        else:
            for p in range(points.shape[0]):
                xb = np.zeros((len(dorder), len(self.bx[     bxi][dorder[0][0]])), field_values.dtype)
                yb = np.zeros((len(dorder), len(self.bx[ygrid[p]][dorder[0][1]])), field_values.dtype)
                zb = np.zeros((len(dorder), len(self.bx[     bxi][dorder[0][2]])), field_values.dtype)
                for o in range(len(dorder)):
                    xb[o] = np.array([self.bx[     bxi][dorder[o][0]][k](xfrac[p])
                                   for k in range(len(self.bx[     bxi][dorder[o][0]]))]).astype(field_values.dtype)
                    yb[o] = np.array([self.by[ygrid[p]][dorder[o][1]][k](yfrac[p])
                                   for k in range(len(self.by[ygrid[p]][dorder[o][1]]))]).astype(field_values.dtype)
                    zb[o] = np.array([self.bz[     bzi][dorder[o][2]][k](zfrac[p])
                                   for k in range(len(self.bz[     bzi][dorder[o][2]]))]).astype(field_values.dtype)
                    #print ['{0:6}'.format(yb[k]) for k in range(2*self.n + 2)]
                    #print ['{0:6}'.format(field_points[p, 0, k, 0, 1]) for k in range(2*self.n + 2)]
                    #print ['{0:6}'.format(field_values[p, 0, k, 0, 1]) for k in range(2*self.n + 2)]
                result[:, p] = np.einsum('okjil,oi,oj,ok->ol', field_values[None, p], xb, yb, zb)
        return result
    def refine_grid(
            self,
            data = None,
            i0 = None, i1 = None,
            j0 = None, j1 = None,
            k0 = None, k1 = None,
            dorder = [(0, 0, 0)],
            factor = 2):
        """
            meant to be called for regularly spaced data, otherwise results make no sense.
        """
        beta_vals = np.empty((len(dorder), 3, factor, len(self.bx[0][0])), dtype = data.dtype)
        for o in range(len(dorder)):
            for i in range(factor):
                beta_vals[o, 0, i] = np.array([self.bx[0][dorder[o][0]][k](i*1./factor)
                                               for k in range(len(self.bx[0][0]))])
                beta_vals[o, 1, i] = np.array([self.bx[0][dorder[o][1]][k](i*1./factor)
                                               for k in range(len(self.bx[0][0]))])
                beta_vals[o, 2, i] = np.array([self.bx[0][dorder[o][2]][k](i*1./factor)
                                               for k in range(len(self.bx[0][0]))])
        if len(data.shape) == 3:
            result = np.empty((len(dorder), (k1 - k0)*factor, (j1 - j0)*factor, (i1 - i0)*factor), dtype = data.dtype)
            for cx in range(factor):
                for cy in range(factor):
                    for cz in range(factor):
                        result[:, cz:result.shape[1]:factor, cy:result.shape[2]:factor, cx:result.shape[3]:factor] = sum(sum(sum(
                                data     [None, k0+kk-self.n:k1+kk-self.n, j0+jj-self.n:j1+jj-self.n, i0+ii-self.n:i1+ii-self.n]
                              * beta_vals[   :, 0,     None,        None,        None, cx, ii] for ii in range(len(self.bx[0][0])))
                              * beta_vals[   :, 1,     None,        None,        None, cy, jj] for jj in range(len(self.bx[0][0])))
                              * beta_vals[   :, 2,     None,        None,        None, cz, kk] for kk in range(len(self.bx[0][0])))
        elif len(data.shape) == 4:
            result = np.empty((len(dorder), (k1 - k0)*factor, (j1 - j0)*factor, (i1 - i0)*factor, 3), dtype = data.dtype)
            for cx in range(factor):
                for cy in range(factor):
                    for cz in range(factor):
                        for coord in range(3):
                            for o in range(len(dorder)):
                                tmp = correlate1d(data[:, :, :, coord], np.array(beta_vals[o, 0, cx, :]), axis = 2)
                                tmp = correlate1d(                 tmp, np.array(beta_vals[o, 1, cy, :]), axis = 1)
                                tmp = correlate1d(                 tmp, np.array(beta_vals[o, 2, cz, :]), axis = 0)
                                result[ o,
                                        cz:result.shape[1]:factor,
                                        cy:result.shape[2]:factor,
                                        cx:result.shape[3]:factor,
                                        coord] = tmp[self.n:result.shape[1]+self.n,
                                                     self.n:result.shape[2]+self.n,
                                                     self.n:result.shape[3]+self.n]
                        #result[:, cz:result.shape[1]:factor, cy:result.shape[2]:factor, cx:result.shape[3]:factor] = sum(sum(sum(
                        #        data     [None, k0+kk-self.n:k1+kk-self.n, j0+jj-self.n:j1+jj-self.n, i0+ii-self.n:i1+ii-self.n,            :]
                        #      * beta_vals[   :, 0,     None,        None,        None, cx, ii, None] for ii in range(len(self.bx[0][0])))
                        #      * beta_vals[   :, 1,     None,        None,        None, cy, jj, None] for jj in range(len(self.bx[0][0])))
                        #      * beta_vals[   :, 2,     None,        None,        None, cz, kk, None] for kk in range(len(self.bx[0][0])))
        return result

