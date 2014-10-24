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
import numpy as np
import ctypes
import inspect

import pyJHTDB

class ThresholdInfo(ctypes.Structure):
    _fields_ = [('x', ctypes.c_int),
                ('y', ctypes.c_int),
                ('z', ctypes.c_int),
                ('value', ctypes.c_float)]

class libJHTDB(object):
    def __init__(self,
            auth_token = pyJHTDB.auth_token):
        self.libname = 'libJHTDB'
        lib_location = os.path.dirname(inspect.getfile(pyJHTDB))
        self.lib = np.ctypeslib.load_library(self.libname, os.path.abspath(os.path.join(lib_location, os.path.pardir)))
        self.authToken = ctypes.c_char_p(auth_token.encode('ascii'))
        self.connection_on = False
        self.hdf5_file_list = []
        self.hdf5_file_desc = {}
        return None
    def initialize(self, exit_on_error = True):
        #initialize gSOAP
        self.lib.soapinit()
        if exit_on_error:
            #enable exit on error
            self.lib.turblibSetExitOnError(ctypes.c_int(1))
        self.connection_on = True
        return None
    def finalize(self):
        #free gSOAP resources
        self.lib.soapdestroy()
        self.connection_on = False
        return None
    def add_hdf5_file(self, filename):
        if pyJHTDB.found_h5py and (not filename in self.hdf5_file_list):
            self.hdf5_file_list.append(filename)
            data = pyJHTDB.h5py.File(filename + '.h5', mode = 'r')
            self.hdf5_file_desc[filename] = {}
            for key in ['_contents', '_dataset', '_size', '_start']:
                self.hdf5_file_desc[filename][key] = data[key][:]
            data.close()
            return self.lib.turblibAddLocalSource(ctypes.c_char_p((filename + '.h5').encode('ascii')))
        else:
            return 0
    def getData(self,
            time, point_coords,
            sinterp = 0, tinterp = 0,
            data_set = 'isotropic1024coarse',
            getFunction = 'getVelocity',
            make_modulo = False):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (point_coords.shape[-1] == 3):
            print ('wrong number of values for coordinates in getData')
            sys.exit()
            return None
        if not (point_coords.dtype == np.float32):
            print('point coordinates in getData must be floats. stopping.')
            sys.exit()
            return None
        npoints = point_coords.shape[0]
        for i in range(1, len(point_coords.shape)-1):
            npoints *= point_coords.shape[i]
        if make_modulo:
            pcoords = np.zeros(point_coords.shape, np.float64)
            pcoords[:] = point_coords
            np.mod(pcoords, 2*np.pi, point_coords)
        get_data = getattr(self.lib, getFunction)
        if getFunction in ['getVelocity',
                           'getForce',
                           'getMagneticField',
                           'getMagneticFieldDebug',
                           'getBunit',
                           'getVectorPotential',
                           'getPressureGradient',
                           'getVelocityLaplacian',
                           'getMagneticFieldLaplacian',
                           'getVectorPotentialLaplacian',
                           'getPressureGradient']:
            result_dim = 3
        elif getFunction in ['getVelocityAndPressure']:
            result_dim = 4
        elif getFunction in ['getPressureHessian']:
            result_dim = 6
        elif getFunction in ['getVelocityGradient',
                             'getMagneticFieldGradient',
                             'getVectorPotentialGradient']:
            result_dim = 9
        elif getFunction in ['getVelocityHessian',
                             'getMagneticFieldHessian',
                             'getVectorPotentialHessian']:
            result_dim = 18
        else:
            print(('wrong result type requested in getData\n'
                 + 'maybe it\'s just missing from the list?'))
            sys.exit()
            return None
        newshape = list(point_coords.shape[0:len(point_coords.shape)-1])
        newshape.append(result_dim)
        result_array = np.empty(newshape, dtype=np.float32)
        get_data(self.authToken,
                 ctypes.c_char_p(data_set.encode('ascii')),
                 ctypes.c_float(time),
                 ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                 point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                 result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array
    def getBoxFilter(self,
            time, point_coords,
            data_set = 'isotropic1024coarse',
            make_modulo = False,
            field = 'velocity',
            filter_width = 7*2*np.pi / 1024):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (point_coords.shape[-1] == 3):
            print ('wrong number of values for coordinates in getBoxFilter')
            sys.exit()
            return None
        if not (point_coords.dtype == np.float32):
            print('point coordinates in getBoxFilter must be floats. stopping.')
            sys.exit()
            return None
        npoints = point_coords.shape[0]
        for i in range(1, len(point_coords.shape)-1):
            npoints *= point_coords.shape[i]
        if make_modulo:
            pcoords = np.zeros(point_coords.shape, np.float64)
            pcoords[:] = point_coords
            np.mod(pcoords, 2*np.pi, point_coords)
        result_array = point_coords.copy()
        self.lib.getBoxFilter(self.authToken,
                 ctypes.c_char_p(data_set.encode('ascii')),
                 ctypes.c_char_p(field.encode('ascii')),
                 ctypes.c_float(time),
                 ctypes.c_float(filter_width),
                 ctypes.c_int(npoints),
                 point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                 result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array
    def getThreshold(
            self,
            data_set = 'isotropic1024coarse',
            field = 'vorticity',
            time = 0.1,
            threshold = 0.0,
            cx = 0, cy = 0, cz = 0,
            nx = 4, ny = 4, nz = 4,
            sinterp = 40,
            tinterp = 0):
        result = ctypes.POINTER(ThresholdInfo)()
        result_size = ctypes.c_int()
        call_result = self.lib.getThreshold(
                self.authToken,
                ctypes.c_char_p(data_set.encode('ascii')),
                ctypes.c_char_p(field.encode('ascii')),
                ctypes.c_float(time),
                ctypes.c_float(threshold),
                ctypes.c_int32(sinterp),
                ctypes.c_int32(cx),
                ctypes.c_int32(cy),
                ctypes.c_int32(cz),
                ctypes.c_int32(nx),
                ctypes.c_int32(ny),
                ctypes.c_int32(nz),
                ctypes.byref(result),
                ctypes.byref(result_size))
        data_type = np.dtype([('x', np.int32),
                              ('y', np.int32),
                              ('z', np.int32),
                              ('value', np.float32)])
        data = np.empty((result_size.value,), data_type)
        for i in range(result_size.value):
            data[i]['x'] = result[i].x
            data[i]['y'] = result[i].y
            data[i]['z'] = result[i].z
            data[i]['value'] = result[i].value
        return data
    def getPosition(self,
            starttime = 0.0,
            endtime = 0.1,
            dt = 0.0004,
            point_coords = None,
            sinterp = 4,
            steps_to_keep = 1,
            data_set = 'isotropic1024coarse'):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (point_coords.shape[-1] == 3):
            print ('wrong number of values for coordinates in getPosition')
            sys.exit()
            return None
        npoints = point_coords.shape[0]
        result_array = np.empty((npoints, 3), dtype=np.float32)
        traj_array = np.empty((steps_to_keep+1, npoints, 3), dtype = np.float32)
        traj_array[0] = point_coords
        time_array = np.linspace(starttime, endtime, steps_to_keep+1).astype(np.float32)
        time_array[0] = starttime
        evolver = self.lib.getPosition
        if data_set == 'custom':
            evolver = self.lib.getCustomPosition
        print('starting integration loop, dataset is ', data_set)
        for tstep in range(1, steps_to_keep + 1):
            print('at time step {0} out of {1}'.format(tstep, steps_to_keep))
            pcoords = traj_array[tstep - 1].copy()
            evolver(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time_array[tstep-1]),
                    ctypes.c_float(time_array[tstep  ]),
                    ctypes.c_float(dt),
                    ctypes.c_int(sinterp), ctypes.c_int(npoints),
                    pcoords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
            print('got next position for time step {0}'.format(tstep))
            traj_array[tstep] = result_array
        return traj_array, time_array
    def getFilteredPosition(self,
            starttime = 0.0,
            endtime = 0.1,
            dt = 0.0004,
            filterwidth = (2*np.pi / 1024) * 5,
            point_coords = None,
            sinterp = 4,
            steps_to_keep = 1,
            data_set = 'isotropic1024coarse'):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (point_coords.shape[-1] == 3):
            print ('wrong number of values for coordinates in getPosition')
            sys.exit()
            return None
        npoints = point_coords.shape[0]
        result_array = np.empty((npoints, 3), dtype=np.float32)
        traj_array = np.empty((steps_to_keep+1, npoints, 3), dtype = np.float32)
        traj_array[0] = point_coords
        time_array = np.empty((steps_to_keep+1), dtype = np.float32)
        integration_time = (endtime - starttime) / steps_to_keep
        time_array[0] = starttime
        evolver = self.lib.getFilteredPosition
        for tstep in range(1, steps_to_keep + 1):
            print('at time step {0} out of {1}'.format(tstep, steps_to_keep))
            pcoords = traj_array[tstep - 1].copy()
            evolver(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(starttime + (tstep - 1)*integration_time),
                    ctypes.c_float(starttime +  tstep     *integration_time),
                    ctypes.c_float(dt),
                    ctypes.c_float(filterwidth),
                    ctypes.c_int(npoints),
                    pcoords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
            traj_array[tstep] = result_array
            time_array[tstep] = starttime + tstep * integration_time
        return traj_array, time_array
    def getBlines(self,
            time = 0.0,
            ds = 0.0004,
            S = 0.1,
            points = np.array([[0, 0, 0]]).astype(np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024'):
        """ use Heun solver to integrate magnetic field lines.
            better solvers can be written in the future.
        """
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (points.shape[1] == 3 and len(points.shape) == 2):
            print ('wrong shape of initial condition in getBlines')
            sys.exit()
            return None
        iterations = int(S / ds)
        npoints = points.shape[0]
        history = np.zeros(
                (iterations+1, npoints, 3),
                dtype = np.float32)
        history[0] = points
        result_array = points.copy()
        get_data = getattr(self.lib, 'getMagneticField')
        def getBunit(point_coords):
            get_data(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time),
                    ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                    point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
            bsize = np.maximum(np.sqrt(np.sum(result_array**2, axis = 1)), 1e-18)
            return result_array / bsize[:, None]
        for s in range(1, iterations + 1):
            bhat0 = getBunit(history[s-1])
            y = history[s-1] + ds * bhat0
            history[s] = history[s-1] + .5*ds * (bhat0 + getBunit(y))
        return history
    def getBline(self,
            time = 0.0,
            ds = 0.0004,
            S = 0.1,
            point = np.array([0, 0, 0]).astype(np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024',
            out_of_domain = None):
        """ use Heun solver to integrate magnetic 1 magnetic field line, within a given domain.
            better solvers can be written in the future.
        """
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (point.shape[0] == 3 and len(point.shape) == 1):
            print ('wrong shape of initial condition in getBline')
            sys.exit()
            return None
        iterations = int(S / abs(ds))
        history = np.zeros(
                (iterations+1, 3),
                dtype = np.float32)
        history[0] = point
        result_array = point.copy()
        get_data = getattr(self.lib, 'getMagneticField')
        def getBunit(point_coords):
            #print point_coords
            get_data(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time),
                    ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(1),
                    point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
            #print 'get_data finished'
            bsize = max(np.sqrt(np.sum(result_array**2)), 1e-18)
            return result_array / bsize
        for s in range(1, iterations + 1):
            if not out_of_domain == None:
                if out_of_domain(history[s-1]):
                    return history[:s]
            y = history[s-1].copy()
            bhat0 = getBunit(y)
            y = history[s-1] + ds * bhat0
            history[s] = history[s-1] + .5*ds * (bhat0 + getBunit(y))
        return history
    def get2wayBline(self,
            time = 0.0,
            ds = 0.0004,
            S = 0.1,
            point = np.array([0, 0, 0]).astype(np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024',
            out_of_domain = None):
        l0 = self.getBline(
                time = time,
                ds = ds,
                S = S,
                point = point,
                sinterp = sinterp,
                tinterp = tinterp,
                data_set = data_set,
                out_of_domain = out_of_domain)
        l1 = self.getBline(
                time = time,
                ds = -ds,
                S = S,
                point = point,
                sinterp = sinterp,
                tinterp = tinterp,
                data_set = data_set,
                out_of_domain = out_of_domain)
        return np.concatenate((l1[::-1], l0[1:]), axis = 0)
    def getBlineAlt(self,
            time, nsteps, ds,
            x0,
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024'):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (x0.shape[1] == 3 and len(x0.shape) == 2):
            print(('wrong shape of initial condition in getBlineAlt, ', x0.shape))
            sys.exit()
            return None
        npoints = x0.shape[0]
        result_array = np.empty((nsteps+1, npoints, 3), dtype=np.float32)
        result_array[0] = x0
        self.lib.getBline(
                self.authToken,
                ctypes.c_char_p(data_set.encode('ascii')),
                ctypes.c_float(time),
                ctypes.c_int(nsteps),
                ctypes.c_float(ds),
                ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array
    def getBlineSphereBounded(self,
            time = 0,
            S = 0.05,
            ds = 0.005,
            x0 = np.array([[0, 0, 0]], dtype = np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024',
            origin = [0, 0, 0],
            radius = 1,
            both_ways = False):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (x0.shape[1] == 3 and len(x0.shape) == 2):
            print(('wrong shape of initial condition in getBlineAlt, ', x0.shape))
            sys.exit()
            return None
        nsteps = int(S / abs(ds))
        npoints = x0.shape[0]
        result_array = np.empty((npoints, nsteps+1, 3), dtype=np.float32)
        length_array = np.empty((npoints), dtype = np.int32)
        result_array[:, 0] = x0
        self.lib.getSphericalBoundedBline(
                self.authToken,
                ctypes.c_char_p(data_set.encode('ascii')),
                ctypes.c_float(time),
                ctypes.c_int(nsteps),
                ctypes.c_float(ds),
                ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ctypes.c_float(origin[0]),
                ctypes.c_float(origin[1]),
                ctypes.c_float(origin[2]),
                ctypes.c_float(radius))
        bline_list = [result_array[p, :length_array[p]].copy()
                      for p in range(x0.shape[0])]
        if both_ways:
            result_array[:, 0] = x0
            self.lib.getSphericalBoundedBline(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time),
                    ctypes.c_int(nsteps),
                    ctypes.c_float(-ds),
                    ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_float(origin[0]),
                    ctypes.c_float(origin[1]),
                    ctypes.c_float(origin[2]),
                    ctypes.c_float(radius))
            for l in range(len(bline_list)):
                bline_list[l] = np.concatenate((result_array[l, 1:length_array[l]][::-1].copy(), bline_list[l]), axis = 0)
        return bline_list
    def getBlineSphereBoundedDebug(self,
            time = 0,
            S = 0.05,
            ds = 0.005,
            x0 = np.array([[0, 0, 0]], dtype = np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024',
            origin = [0, 0, 0],
            radius = 1,
            both_ways = False):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (x0.shape[1] == 3 and len(x0.shape) == 2):
            print(('wrong shape of initial condition in getBlineAlt, ', x0.shape))
            sys.exit()
            return None
        nsteps = int(S / abs(ds))
        npoints = x0.shape[0]
        result_array = np.empty((npoints, nsteps+1, 3), dtype=np.float32)
        length_array = np.empty((npoints), dtype = np.int32)
        result_array[:, 0] = x0
        self.lib.getSphericalBoundedBlineDebug(
                self.authToken,
                ctypes.c_char_p(data_set.encode('ascii')),
                ctypes.c_float(time),
                ctypes.c_int(nsteps),
                ctypes.c_float(ds),
                ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ctypes.c_float(origin[0]),
                ctypes.c_float(origin[1]),
                ctypes.c_float(origin[2]),
                ctypes.c_float(radius))
        bline_list = [result_array[p, :length_array[p]].copy()
                      for p in range(x0.shape[0])]
        if both_ways:
            result_array[:, 0] = x0
            self.lib.getSphericalBoundedBlineDebug(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time),
                    ctypes.c_int(nsteps),
                    ctypes.c_float(-ds),
                    ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_float(origin[0]),
                    ctypes.c_float(origin[1]),
                    ctypes.c_float(origin[2]),
                    ctypes.c_float(radius))
            for l in range(len(bline_list)):
                bline_list[l] = np.concatenate((result_array[l, 1:length_array[l]][::-1].copy(), bline_list[l]), axis = 0)
        return bline_list
    def getBlineRectBounded(self,
            time = 0,
            S = 0.05,
            ds = 0.005,
            x0 = np.array([[0, 0, 0]], dtype = np.float32),
            sinterp = 4,
            tinterp = 0,
            data_set = 'mhd1024',
            xmin = 0, xmax = 1,
            ymin = 0, ymax = 1,
            zmin = 0, zmax = 1,
            both_ways = False):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()
        if not (x0.shape[1] == 3 and len(x0.shape) == 2):
            print(('wrong shape of initial condition in getBlineAlt, ', x0.shape))
            sys.exit()
            return None
        nsteps = int(S / abs(ds))
        npoints = x0.shape[0]
        result_array = np.empty((npoints, nsteps+1, 3), dtype=np.float32)
        length_array = np.empty((npoints), dtype = np.int32)
        result_array[:, 0] = x0
        self.lib.getRectangularBoundedBline(
                self.authToken,
                ctypes.c_char_p(data_set.encode('ascii')),
                ctypes.c_float(time),
                ctypes.c_int(nsteps),
                ctypes.c_float(ds),
                ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ctypes.c_float(xmin), ctypes.c_float(xmax),
                ctypes.c_float(ymin), ctypes.c_float(ymax),
                ctypes.c_float(zmin), ctypes.c_float(zmax))
        bline_list = [result_array[p, :length_array[p]].copy()
                      for p in range(x0.shape[0])]
        if both_ways:
            result_array[:, 0] = x0
            self.lib.getRectangularBoundedBline(
                    self.authToken,
                    ctypes.c_char_p(data_set.encode('ascii')),
                    ctypes.c_float(time),
                    ctypes.c_int(nsteps),
                    ctypes.c_float(-ds),
                    ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                    result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                    length_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_float(xmin), ctypes.c_float(xmax),
                    ctypes.c_float(ymin), ctypes.c_float(ymax),
                    ctypes.c_float(zmin), ctypes.c_float(zmax))
            for l in range(len(bline_list)):
                bline_list[l] = np.concatenate((result_array[l, 1:length_array[l]][::-1].copy(), bline_list[l]), axis = 0)
        return bline_list

