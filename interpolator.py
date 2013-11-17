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
        self.dx = self.info['xnodes'][1:] - self.info['xnodes'][0:self.info['nx']-1]
        self.dy = self.info['ynodes'][1:] - self.info['ynodes'][0:self.info['ny']-1]
        self.dz = self.info['znodes'][1:] - self.info['znodes'][0:self.info['nz']-1]
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
        xgrid = np.searchsorted(self.info['xnodes'], points[:, 0]) - 1
        ygrid = np.searchsorted(self.info['ynodes'], points[:, 1]) - 1
        zgrid = np.searchsorted(self.info['znodes'], points[:, 2]) - 1
        xfrac = (points[:, 0] - self.info['xnodes'][xgrid])/self.dx[xgrid]
        yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.dy[ygrid]
        zfrac = (points[:, 2] - self.info['znodes'][zgrid])/self.dz[zgrid]
        field_points = np.zeros((points.shape[0], 2*self.n+2, 2*self.n+2, 2*self.n+2, 3), dtype = np.float32)
        xb = np.zeros((points.shape[0], 2*self.n+2), dtype = np.float32)
        yb = np.zeros((points.shape[0], 2*self.n+2), dtype = np.float32)
        zb = np.zeros((points.shape[0], 2*self.n+2), dtype = np.float32)
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 0] = self.info['xnodes'][np.newaxis, np.newaxis, xgrid[p]-self.n:xgrid[p]+self.n+2]
            field_points[p, :, :, :, 1] = self.info['ynodes'][np.newaxis, ygrid[p]-self.n:ygrid[p]+self.n+2, np.newaxis]
            field_points[p, :, :, :, 2] = self.info['znodes'][zgrid[p]-self.n:zgrid[p]+self.n+2, np.newaxis, np.newaxis]
        print 'computed points where field is needed, now getting values from DB'
        field_values = lTDB.getData(
                time, field_points,
                sinterp = 0, tinterp = 0,
                data_set = self.info['name'],
                getFunction = getFunction)
        print 'got values from DB, now interpolating'
        result = np.zeros((len(dorder), points.shape[0], field_values.shape[-1]), dtype = np.float32)
        bxi = 0
        if self.info['yuniform']:
            ygrid[p] = 0
        bzi = 0
        for p in range(points.shape[0]):
            for o in range(len(dorder)):
                xb[p] = np.array([self.bx[     bxi][dorder[o][0]][k](xfrac[p])
                                  for k in range(self.spline['x'].N)])
                yb[p] = np.array([self.by[ygrid[p]][dorder[o][1]][k](yfrac[p])
                                  for k in range(self.spline['y'].N)])
                zb[p] = np.array([self.bz[     bzi][dorder[o][2]][k](zfrac[p])
                                  for k in range(self.spline['z'].N)])
                result[o, p] = sum(sum(sum(field_values[p, k, j, i]*xb[p, i]
                                           for i in range(self.spline['x'].N)
                                                                  )*yb[p, j]
                                       for j in range(self.spline['y'].N)
                                                                  )*zb[p, k]
                                   for k in range(self.spline['z'].N))
        return result

