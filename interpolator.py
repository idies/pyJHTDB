import numpy as np
import os
import pickle

import pyTDB
import pyTDB.generic_splines as gs


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
        if os.path.exists(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m)):
            self.spline = pickle.load(open(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m),
                                           'r'))
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
                        open(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m),
                             'w'))
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
                    + np.array(range(-self.n, self.n+2))*self.info['dx'])[np.newaxis, np.newaxis, :]
        if self.info['yperiodic']:
            ygrid = np.floor(points[:, 1] / self.info['dy']).astype(np.int) % self.info['ny']
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy']
            for p in range(points.shape[0]):
                field_points[p, :, :, :, 1] = (self.info['ynodes'][ygrid[p]]
                        + np.array(range(-self.n, self.n+2))*self.info['dy'])[np.newaxis, :, np.newaxis]
        else:
            ygrid = np.searchsorted(self.info['ynodes'], points[:, 1]).astype(np.int) - 1
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy'][ygrid]
            for p in range(points.shape[0]):
                if ygrid[p] < 0 or ygrid[p] > self.info['ny'] - 1:
                    return None
                elif ygrid[p] < self.n:
                    field_points[p, :, :2*self.n+1, :, 1] = self.info['ynodes'][np.newaxis, :2*self.n+1, np.newaxis]
                elif ygrid[p] > self.info['ny'] - self.n - 1:
                    field_points[p, :, :2*self.n+1, :, 1] = self.info['ynodes'][np.newaxis, self.info['ny'] - (2*self.n+1):, np.newaxis]
                else:
                    field_points[p, :, :, :, 1] = self.info['ynodes'][np.newaxis, ygrid[p]-self.n:ygrid[p]+self.n+2, np.newaxis]
        zgrid = np.floor(points[:, 2] / self.info['dz']).astype(np.int) % self.info['nz']
        zfrac = (points[:, 2] - self.info['znodes'][zgrid])/self.info['dz']
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 2] = (self.info['znodes'][zgrid[p]]
                    + np.array(range(-self.n, self.n+2))*self.info['dz'])[:, np.newaxis, np.newaxis]
        print 'computed points where field is needed, now getting values from DB'
        ## I could in principle call getRaw[...] for each point,
        ## but that would mean a lot of calls to the DB,
        ## and we should avoid that due to the latency.
        field_values = lTDB.getData(
                time, field_points,
                sinterp = 0, tinterp = 0,
                data_set = self.info['name'],
                getFunction = getFunction)
        print 'got values from DB, now interpolating'
        result = np.zeros((len(dorder), points.shape[0], field_values.shape[-1]), dtype = np.float32)
        bxi = 0
        if self.info['yperiodic']:
            ygrid[:] = 0
        bzi = 0
        for p in range(points.shape[0]):
            for o in range(len(dorder)):
                xb = np.array([self.bx[     bxi][dorder[o][0]][k](xfrac[p])
                               for k in range(len(self.bx[     bxi][dorder[o][0]]))]).astype(field_values.dtype)
                yb = np.array([self.by[ygrid[p]][dorder[o][1]][k](yfrac[p])
                               for k in range(len(self.by[ygrid[p]][dorder[o][1]]))]).astype(field_values.dtype)
                zb = np.array([self.bz[     bzi][dorder[o][2]][k](zfrac[p])
                               for k in range(len(self.bz[     bzi][dorder[o][2]]))]).astype(field_values.dtype)
                #print ['{0:6}'.format(yb[k]) for k in range(2*self.n + 2)]
                #print ['{0:6}'.format(field_points[p, 0, k, 0, 1]) for k in range(2*self.n + 2)]
                #print ['{0:6}'.format(field_values[p, 0, k, 0, 1]) for k in range(2*self.n + 2)]
                result[o, p] = np.einsum('kjil,i,j,k->l', field_values[p], xb, yb, zb)
        return result

