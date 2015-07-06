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

import os
import sys
if sys.version_info < (3,):
    import cPickle as pickle
else:
    import pickle
import numpy as np
import gzip
import ctypes as ct

import distutils
import distutils.command
import distutils.command.build_ext
import distutils.core
import distutils.dist
import distutils.log
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
            compute_fast_beta = False,
            nx = None,
            mx = None,
            ny = None,
            my = None,
            nz = None,
            mz = None,
            initialize = True,
            cformula_unroll = False):
        self.nx = int(np.floor(n)) if (type(nx) == type(None)) else int(np.floor(nx))
        self.ny = int(np.floor(n)) if (type(ny) == type(None)) else int(np.floor(ny))
        self.nz = int(np.floor(n)) if (type(nz) == type(None)) else int(np.floor(nz))
        # for backwards compatibility, put in an n member as well
        # I'm using the max, since the point is to use the value for buffers of 3D
        # fields.
        self.n = max(self.nx, self.ny, self.nz)
        self.mx = int(np.floor(m)) if (type(mx) == type(None)) else int(np.floor(mx))
        self.my = int(np.floor(m)) if (type(my) == type(None)) else int(np.floor(my))
        self.mz = int(np.floor(m)) if (type(mz) == type(None)) else int(np.floor(mz))
        self.m = min(self.mx, self.my, self.mz)
        self.info = info
        self.clib_loaded = False
        self.cformula_unroll = cformula_unroll
        self.initialized = False
        if initialize:
            self.initialize(compute_fast_beta = compute_fast_beta)
        return None
    def initialize(self, compute_fast_beta = False):
        pickle_file = {
                'x' : os.path.join(
                        pyJHTDB.lib_folder,
                        (self.info['name'] + '_spline_' +
                         'xn{0}m{1}'.format(self.nx, self.mx) +
                         '.py{0}{1}'.format(sys.version_info[0],
                                            sys.version_info[1]) +
                         '.pickle.gz')),
                'y' : os.path.join(
                        pyJHTDB.lib_folder,
                        (self.info['name'] + '_spline_' +
                         'yn{0}m{1}'.format(self.ny, self.my) +
                         '.py{0}{1}'.format(sys.version_info[0],
                                            sys.version_info[1]) +
                         '.pickle.gz')),
                'z' : os.path.join(
                        pyJHTDB.lib_folder,
                        (self.info['name'] + '_spline_' +
                         'zn{0}m{1}'.format(self.nz, self.mz) +
                         '.py{0}{1}'.format(sys.version_info[0],
                                            sys.version_info[1]) +
                         '.pickle.gz'))}
        self.spline = {}
        for coord in ['x', 'y', 'z']:
            if os.path.exists(pickle_file[coord]):
                self.spline[coord] = pickle.load(gzip.open(pickle_file[coord]))
            else:
                if self.info[coord + 'uniform'] and self.info[coord + 'periodic']:
                    self.spline[coord] = gs.generic_spline_1D(
                            self.info[coord + 'nodes'][:2],
                            max_deriv = getattr(self, 'm' + coord),
                            neighbours = getattr(self, 'n' + coord),
                            period = self.info['l' + coord])
                else:
                    self.spline[coord] = gs.generic_spline_1D(
                            self.info[coord + 'nodes'],
                            max_deriv = getattr(self, 'm' + coord),
                            neighbours = getattr(self, 'n' + coord))
                self.spline[coord].compute_derivs()
                self.spline[coord].compute_beta()
                pickle.dump(self.spline[coord],
                            gzip.open(pickle_file[coord],
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
        self.initialized = True
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
        if (not self.initialized):
            self.initialize()
        field_points = np.zeros((points.shape[0], 2*self.nz+2, 2*self.ny+2, 2*self.nx+2, 3), dtype = np.float32)
        xgrid = np.floor(points[:, 0] / self.info['dx']).astype(np.int) % self.info['nx']
        xfrac = (points[:, 0] - self.info['xnodes'][xgrid])/self.info['dx']
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 0] = (self.info['xnodes'][xgrid[p]]
                    + np.array(list(range(-self.nx, self.nx+2)))*self.info['dx'])[np.newaxis, np.newaxis, :]
        if self.info['yperiodic']:
            ygrid = np.floor(points[:, 1] / self.info['dy']).astype(np.int) % self.info['ny']
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy']
            for p in range(points.shape[0]):
                field_points[p, :, :, :, 1] = (self.info['ynodes'][ygrid[p]]
                        + np.array(list(range(-self.ny, self.ny+2)))*self.info['dy'])[np.newaxis, :, np.newaxis]
        else:
            ygrid = np.searchsorted(self.info['ynodes'], points[:, 1]).astype(np.int) - 1
            yfrac = (points[:, 1] - self.info['ynodes'][ygrid])/self.info['dy'][ygrid]
            for p in range(points.shape[0]):
                if ygrid[p] < 0 or ygrid[p] > self.info['ny'] - 1:
                    return None
                elif ygrid[p] < self.ny:
                    field_points[p, :, :2*self.ny+1, :, 1] = self.info['ynodes'][np.newaxis, :2*self.ny+1, np.newaxis]
                elif ygrid[p] >= self.info['ny'] - self.ny - 1:
                    field_points[p, :, :2*self.ny+1, :, 1] = self.info['ynodes'][np.newaxis, self.info['ny'] - (2*self.ny+1):, np.newaxis]
                else:
                    field_points[p, :, :, :, 1] = self.info['ynodes'][np.newaxis, ygrid[p]-self.ny:ygrid[p]+self.ny+2, np.newaxis]
        zgrid = np.floor(points[:, 2] / self.info['dz']).astype(np.int) % self.info['nz']
        zfrac = (points[:, 2] - self.info['znodes'][zgrid])/self.info['dz']
        for p in range(points.shape[0]):
            field_points[p, :, :, :, 2] = (self.info['znodes'][zgrid[p]]
                    + np.array(list(range(-self.nz, self.nz+2)))*self.info['dz'])[:, np.newaxis, np.newaxis]
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
                result[:, p] = np.einsum('okjil,oi,oj,ok->ol', field_values[None, p], xb, yb, zb)
        return result
    def write_coefficients(self):
        if (not self.initialized):
            self.initialize()
        for coord in ['x', 'y', 'z']:
            for order in range(self.m+1):
                text_file = open(
                        (self.info['name']
                        + '_' + coord
                        + 'spline_m{0}q{1:0>2}_d{2}_coeff.csv'.format(
                                                                getattr(self, 'm' + coord),
                                                                getattr(self, 'n' + coord)*2 + 2,
                                                                order)),
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
        self.write_cfile(cfile_name = cfile_name)
        try:
            self.clib = np.ctypeslib.load_library(
                    'lib' +
                    'py{0}{1}'.format(sys.version_info[0],
                                      sys.version_info[1]) +
                    os.path.basename(self.cfile_name),
                    pyJHTDB.lib_folder)
        except:
            builder = distutils.command.build_ext.build_ext(
                    distutils.dist.Distribution({'name' : os.path.basename(self.cfile_name)}))
            builder.extensions = [
                    distutils.core.Extension(
                        'lib' +
                        'py{0}{1}'.format(sys.version_info[0],
                                          sys.version_info[1]) +
                        os.path.basename(self.cfile_name),
                        sources = [self.cfile_name + '.c'])]
            builder.build_lib = os.path.abspath(pyJHTDB.lib_folder)
            builder.build_temp = tempfile.gettempdir()
            builder.swig_opts = []
            distutils.log.set_verbosity(1)
            builder.run()
            self.clib = np.ctypeslib.load_library(
                    'lib' +
                    'py{0}{1}'.format(sys.version_info[0],
                                      sys.version_info[1]) +
                    os.path.basename(self.cfile_name),
                    pyJHTDB.lib_folder)
        self.clib_loaded = True
        return None
    def cinterpolate(
            self,
            x = None,
            f = None,
            diff = [0, 0, 0],
            field_offset = [0, 0, 0],
            debug = False):
        if not self.clib_loaded:
            self.generate_clib()
        diff = np.array(diff).astype(np.int32)
        field_offset = np.array(field_offset).astype(np.int32)
        field_size   = np.array(f.shape[:-1]).astype(np.int32)
        assert(diff.shape[0] == 3 and
               len(diff.shape) == 1)
        assert(field_offset.shape[0] == 3 and
               len(field_offset.shape) == 1)
        assert(f.flags['C_CONTIGUOUS'] and
               f.dtype == np.float32 and
               len(f.shape) == 4)
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
        if debug:
            print(node_array)
            print(field_offset)
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
            cfile_name = None,
            base_cname = None):
        if type(cfile_name) == type(None):
            self.cfile_name = (pyJHTDB.lib_folder +
                          self.info['name'] + '_spline' +
                          '_xm{0}q{1}'.format(self.mx, self.nx*2 + 2) +
                          '_ym{0}q{1}'.format(self.my, self.ny*2 + 2) +
                          '_zm{0}q{1}'.format(self.mz, self.nz*2 + 2))
        else:
            self.cfile_name = cfile_name
        base_xname = 'xm{0}q{1}'.format(self.mx, self.nx*2 + 2)
        base_yname = 'ym{0}q{1}'.format(self.my, self.ny*2 + 2)
        base_zname = 'zm{0}q{1}'.format(self.mz, self.nz*2 + 2)
        if type(base_cname) == type(None):
            self.base_cname = base_xname + '_' + base_yname + '_' + base_zname
        else:
            self.base_cname = base_cname
        if os.path.exists(self.cfile_name + '.c'):
            return None
        if (not self.initialized):
            self.initialize()
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
                        csuffix = (
                            '_' + coord +
                            'm{0}q{1}'.format(
                                getattr(self, 'm' + coord),
                                getattr(self, 'n' + coord)*2 + 2)))
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
                'int point;\n' +
                'int component;\n' +
                'int i0, i1, i2;\n' +
                'float bx[{0}], by[{1}], bz[{2}];\n'.format(
                    self.nx*2+2,
                    self.ny*2+2,
                    self.nz*2+2) +
                'int ix[{0}], iy[{1}], iz[{2}];\n'.format(
                    self.nx*2+2,
                    self.ny*2+2,
                    self.nz*2+2))
        if not self.cformula_unroll:
            src_txt += 'int xcounter, ycounter, zcounter;\n'
        # loop over points
        src_txt += 'for (point = 0; point < npoints; point++)\n{\n'
        # get polynomials
        src_txt += 'xbeta_' + base_xname + '(diff[0], fractions[point*3+0], bx);\n'
        if self.info['yperiodic']:
            src_txt += 'ybeta_' + base_yname + '(diff[1], fractions[point*3+1], by);\n'
        else:
            src_txt += 'ybeta_' + base_yname + '(nodes[3*point+1], diff[1], fractions[point*3+1], by);\n'
        src_txt += 'zbeta_' + base_zname + '(diff[2], fractions[point*3+2], bz);\n'
        src_txt += 'xindices_' + base_xname + '(nodes[3*point+0], ix);\n'
        src_txt += 'yindices_' + base_yname + '(nodes[3*point+1], iy);\n'
        src_txt += 'zindices_' + base_zname + '(nodes[3*point+2], iz);\n'
        # loop over components
        src_txt += 'for (component = 0; component < field_components; component++)\n{\n'
        bx = ['bx[{0}]'.format(i) for i in range(self.nx*2 + 2)]
        by = ['by[{0}]'.format(i) for i in range(self.ny*2 + 2)]
        bz = ['bz[{0}]'.format(i) for i in range(self.nz*2 + 2)]
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
        if self.cformula_unroll:
            def write_interp1D(bname, fname, q):
                tmp_txt  = '('
                for i in range(q-1):
                    tmp_txt += '\n' + bname[i] + '*' + fname[i] + ' + '
                tmp_txt += '\n' + bname[q-1] + '*' + fname[q-1] + ')'
                return tmp_txt
            fzname = []
            for i in range(self.nz*2 + 2):
                fyname = []
                for j in range(self.ny*2 + 2):
                    fxname = []
                    for k in range(self.nx*2 + 2):
                        fxname.append(
                               ('field[(((i2+iz[{0}])*field_size[1]' +
                                     ' + (i1+iy[{1}]))*field_size[2]' +
                                     ' + (i0+ix[{2}]))*field_components + component]').format(i, j, k))
                    fyname.append(write_interp1D(bx, fxname, self.nx*2 + 2) + '\n')
                fzname.append(write_interp1D(by, fyname, self.ny*2 + 2) + '\n')
            src_txt += ('result[field_components*point + component] = ' +
                        write_interp1D(bz, fzname, self.nz*2+2) +
                        ';\n')
        else:
            src_txt += (
                    'result[field_components*point + component] = 0;\n' +
                    'for (zcounter = 0; zcounter < {0}; zcounter++)\n'.format(self.nz*2+2) +
                    'for (ycounter = 0; ycounter < {0}; ycounter++)\n'.format(self.ny*2+2) +
                    'for (xcounter = 0; xcounter < {0}; xcounter++)\n'.format(self.nx*2+2) +
                    'result[field_components*point + component] += ' +
                    'bz[zcounter]*by[ycounter]*bx[xcounter]*' +
                    'field[(((i2 + iz[zcounter]) *field_size[1] + ' +
                            '(i1 + iy[ycounter]))*field_size[2] + ' +
                            '(i0 + ix[xcounter]))*field_components + component];\n')
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
            if (not self.initialized):
                self.initialize()
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

