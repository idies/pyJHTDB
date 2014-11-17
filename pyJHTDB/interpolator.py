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
import ctypes as ct
import distutils
import distutils.command
import distutils.command.build_ext
import distutils.core
import distutils.dist
import tempfile

import pyJHTDB
import pyJHTDB.generic_splines as gs

if pyJHTDB.found_scipy:
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
        pickle_file = os.path.join(
                pyJHTDB.lib_folder,
                info['name'] + '_spline_interpolator_n{0}_m{1}.pickle.gz'.format(n, m))
        if os.path.exists(pickle_file):
            self.spline = pickle.load(gzip.open(pickle_file))
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
                        gzip.open(pickle_file,
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
            getFunction = 'getVelocityAndPressure'):
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
                yb = np.zeros((len(dorder), len(self.by[ygrid[p]][dorder[0][1]])), field_values.dtype)
                zb = np.zeros((len(dorder), len(self.bz[     bxi][dorder[0][2]])), field_values.dtype)
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
    def write_coefficients(self):
        for coord in ['x', 'y', 'z']:
            for order in range(self.m+1):
                text_file = open(
                        (self.info['name']
                        + '_' + coord
                        + 'spline_m{0}q{1:0>2}_d{2}_coeff.csv'.format(self.m, self.n*2 + 2, order)),
                        'w')
                if self.info[coord + 'periodic']:
                    for point in range(len(self.spline[coord].beta[0][order])):
                        text_file.write('0, {0}'.format(self.spline[coord].neighbour_list[0][point]))
                        for c in self.spline[coord].beta[0][order][point].coef:
                            text_file.write(', {0}'.format(c))
                        text_file.write('\r\n')
                else:
                    for node in range(len(self.spline[coord].beta)):
                        for point in range(len(self.spline[coord].beta[node][order])):
                            if (self.spline[coord].beta[node][order][point].coef.shape[0] > 1
                                 or (not (self.spline[coord].beta[node][order][point].coef[0] == 0.0))):
                                text_file.write('{0}, {1}'.format(node, self.spline[coord].neighbour_list[node][point]))
                                for c in self.spline[coord].beta[node][order][point].coef:
                                    text_file.write(', {0}'.format(c))
                                if self.spline[coord].beta[node][order][point].coef.shape[0] < self.m*2 + 2 - order:
                                    for tcounter in range(self.m*2 + 2
                                            - order - self.spline[coord].beta[node][order][point].coef.shape[0]):
                                        text_file.write(', 0')
                                text_file.write('\r\n')
                text_file.close()
        return None
    def generate_clib(
            self,
            cfile_name = None):
        try:
            self.clib = np.ctypeslib.load_library(
                    'lib' + os.path.basename(self.cfile_name),
                    pyJHTDB.lib_folder)
        except:
            self.write_cfile(cfile_name = cfile_name)
            builder = distutils.command.build_ext.build_ext(
                    distutils.dist.Distribution({'name' : os.path.basename(self.cfile_name)}))
            builder.extensions = [
                    distutils.core.Extension(
                        'lib' + os.path.basename(self.cfile_name),
                        sources = [self.cfile_name + '.c'])]
            builder.build_lib = os.path.abspath(pyJHTDB.lib_folder)
            builder.build_temp = tempfile.gettempdir()
            builder.swig_opts = []
            builder.verbose = True
            builder.run()
            self.clib = np.ctypeslib.load_library(
                    'lib' + os.path.basename(self.cfile_name),
                    pyJHTDB.lib_folder)
        return None
    def cinterpolate(
            self,
            x = None,
            f = None,
            diff = [0, 0, 0],
            field_offset = [0, 0, 0]):
        diff = np.array(diff).astype(np.int32)
        field_offset = np.array(field_offset).astype(np.int32)
        field_size   = np.array(f.shape[:-1]).astype(np.int32)
        assert(diff.shape[0] == 3 and
               len(diff.shape) == 1)
        assert(field_offset.shape[0] == 3 and
               len(field_offset.shape) == 1)
        assert(f.flags['C_CONTIGUOUS'] and
               f.dtype == np.float32)
        y = np.ascontiguousarray(x.reshape(-1, 3), np.float32)
        node_array = np.zeros(y.shape, np.int32)
        node_array[:, 0] = np.floor(y[:, 0] / self.info['dx'])
        if self.info['yperiodic']:
            node_array[:, 1] = np.floor(y[:, 1] / self.info['dy'])
        else:
            node_array[:, 1] = np.searchsorted(self.info['ynodes'], y[:, 1], side = 'right') - 1
        node_array[:, 2] = np.floor(y[:, 2] / self.info['dz'])
        frac_array = y.copy()
        frac_array[:, 0] = y[:, 0] / self.info['dx'] - node_array[:, 0]
        if self.info['yperiodic']:
            frac_array[:, 1] = y[:, 1] / self.info['dy'] - node_array[:, 1]
        else:
            frac_array[:, 1] = (y[:, 1] - self.info['ynodes'][node_array[:, 1]]) / self.info['dy'][node_array[:, 1]]
        frac_array[:, 2] = y[:, 2] / self.info['dz'] - node_array[:, 2]
        s = np.ascontiguousarray(np.zeros((y.shape[0], f.shape[-1]), np.float32))
        getattr(self.clib, 'interpolate_' + self.base_cname)(
                frac_array.ctypes.data_as(ct.POINTER(ct.c_float)),
                node_array.ctypes.data_as(ct.POINTER(ct.c_int)),
                ct.c_int(y.shape[0]),
                diff.ctypes.data_as(ct.POINTER(ct.c_int)),
                f.ctypes.data_as(ct.POINTER(ct.c_float)),
                field_offset.ctypes.data_as(ct.POINTER(ct.c_int)),
                field_size.ctypes.data_as(ct.POINTER(ct.c_int)),
                ct.c_int(f.shape[-1]),
                s.ctypes.data_as(ct.POINTER(ct.c_float)))
        return s.reshape(tuple(list(x.shape[:-1]) + [f.shape[-1]]))
    def write_cfile(
            self,
            cfile_name = None,  #'spline_m{0}q{1:0>2}'.format(self.m, self.n*2 + 2),
            base_cname = None): #'m{0}q{1:0>2}'.format(self.m, self.n*2 + 2)):
        if type(cfile_name) == type(None):
            self.cfile_name = (pyJHTDB.lib_folder +
                          self.info['name'] + '_' +
                          'spline_m{0}q{1:0>2}'.format(self.m, self.n*2 + 2))
        else:
            self.cfile_name = cfile_name
        if type(base_cname) == type(None):
            self.base_cname = 'm{0}q{1:0>2}'.format(self.m, self.n*2 + 2)
        else:
            self.base_name = base_name
        if os.path.exists(self.cfile_name + '.c'):
            return None
        def write_interp1D(bname, fname):
            tmp_txt  = '('
            for i in range(self.n*2+1):
                tmp_txt += '\n' + bname[i] + '*' + fname[i] + ' + '
            tmp_txt += '\n' + bname[self.n*2+1] + '*' + fname[self.n*2+1] + ')'
            return tmp_txt
        cfile = open(self.cfile_name + '.c', 'w')
        ### headers
        cfile.write(
                '#include <assert.h>\n' +
                '#include <stdlib.h>\n' +
                '#include <stdio.h>\n' +
                '\n')
        ### functions to compute beta polynomials
        for coord in ['x', 'y', 'z']:
            ## beta polynomial implementation
            cfile.write(
                    self.spline[coord].write_cfunction(
                        cprefix = coord,
                        csuffix = '_' + self.base_cname)
                    + '\n')
        ### write 3D interpolation
        src_txt = (
                'int interpolate_' + self.base_cname + '('
              + 'float *fractions, '
              + 'int *nodes, '
              + 'int npoints, '
              + 'int *diff, '
              + 'float *field, '
              + 'int *field_offset, '
              + 'int *field_size, '
              + 'int field_components, '
              + 'float *result)\n')
        src_txt += '{\n'
        # various variables
        src_txt += (
#                'fprintf(stderr, "entering interpolate %d %d %d\\n", field_offset[0], field_offset[1], field_offset[2]);\n' +
                'int point;\n' +
                'int component;\n' +
                'int i0, i1, i2;\n' +
                'float bx[{0}], by[{0}], bz[{0}];\n'.format(self.n*2+2) +
                'int ix[{0}], iy[{0}], iz[{0}];\n'.format(self.n*2+2))
        # loop over points
        src_txt += 'for (point = 0; point < npoints; point++)\n{\n'
#        src_txt += 'fprintf(stderr, "inside point loop, point is %d\\n", point);\n'
        # get polynomials
        src_txt += 'xbeta_' + self.base_cname + '(diff[0], fractions[point*3+0], bx);\n'
        if self.info['yperiodic']:
            src_txt += 'ybeta_' + self.base_cname + '(diff[1], fractions[point*3+1], by);\n'
        else:
            src_txt += 'ybeta_' + self.base_cname + '(nodes[3*point+1], diff[1], fractions[point*3+1], by);\n'
        src_txt += 'zbeta_' + self.base_cname + '(diff[2], fractions[point*3+2], bz);\n'
        src_txt += 'xindices_' + self.base_cname + '(nodes[3*point+0], ix);\n'
        src_txt += 'yindices_' + self.base_cname + '(nodes[3*point+1], iy);\n'
        src_txt += 'zindices_' + self.base_cname + '(nodes[3*point+2], iz);\n'
        # loop over components
        src_txt += 'for (component = 0; component < field_components; component++)\n{\n'
#        src_txt += 'fprintf(stderr, "inside component loop, component is %d\\n", component);\n'
        bx = ['bx[{0}]'.format(i) for i in range(self.n*2 + 2)]
        by = ['by[{0}]'.format(i) for i in range(self.n*2 + 2)]
        bz = ['bz[{0}]'.format(i) for i in range(self.n*2 + 2)]
        src_txt += (
                'i0 = nodes[3*point + 0] - field_offset[0];\n' +
                'i1 = nodes[3*point + 1] - field_offset[1];\n' +
                'i2 = nodes[3*point + 2] - field_offset[2];\n' +
                'if (i0 < 0 || i1 < 0 || i2 < 0)' +
                '{\n' +
                    'fprintf(stderr, "negative indices in interpolate %d %d %d\\n", i0, i1, i2);\n' +
                    'fprintf(stderr, "exiting interpolate now, results are most likely nonsensical\\n");\n' +
                    'return EXIT_FAILURE;\n' +
                '}\n')
        fzname = []
        for i in range(self.n*2 + 2):
            fyname = []
            for j in range(self.n*2 + 2):
                fxname = []
                for k in range(self.n*2 + 2):
                    fxname.append(
                           ('field[(((i2+iz[{0}])*field_size[1]' +
                                 ' + (i1+iy[{1}]))*field_size[2]' +
                                 ' + (i0+ix[{2}]))*field_components + component]').format(i, j, k))
                fyname.append(write_interp1D(bx, fxname) + '\n')
            fzname.append(write_interp1D(by, fyname) + '\n')
        src_txt += 'result[field_components*point + component] = ' + write_interp1D(bz, fzname) + ';\n'
        src_txt += '}\n'                                            # close component loop
        src_txt += '}\n'                                            # close point loop
        src_txt += 'return EXIT_SUCCESS;\n}\n'                      # close function
        cfile.write(src_txt)
        cfile.close()
        return None
    if pyJHTDB.found_scipy:
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

