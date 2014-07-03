"""
    Python tools and wrappers for the Johns Hopkins Turbulence Databases C library.
"""

import os
import sys
import numpy as np
import ctypes
import inspect

class libTDB:
    def __init__(self,
            libname = 'libTDB',
            libdir = '.',
            srcdir = None,
            auth_token = 'edu.jhu.pha.turbulence.testing-201302'):
        self.libname = libname
        self.libdir = libdir
        if srcdir == None:
            self.srcdir = libdir
        else:
            self.srcdir = srcdir
        self.lib = np.ctypeslib.load_library(libname, libdir)
        self.authToken = ctypes.c_char_p(auth_token)
        self.connection_on = False
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
        return self.lib.turblibAddLocalSource(ctypes.c_char_p(filename + '.h5'))
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
            print 'point coordinates in getData must be floats. stopping.'
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
                           'getBunit',
                           'getVectorPotential',
                           'getPressureGradient',
                           'getVelocityLaplacian',
                           'getMagneticFieldLaplacian',
                           'getVectorPotentialLaplacian',
                           'getPressureGradient',
                           'getVelocitySoap',
                           'getForceSoap',
                           'getMagneticFieldSoap',
                           'getVectorPotentialSoap',
                           'getPressureGradientSoap',
                           'getVelocityLaplacianSoap',
                           'getMagneticFieldLaplacianSoap',
                           'getVectorPotentialLaplacianSoap',
                           'getPressureGradientSoap']:
            result_dim = 3
        elif getFunction in ['getVelocityAndPressure',
                             'getVelocityAndPressureSoap']:
            result_dim = 4
        elif getFunction in ['getPressureHessian',
                             'getPressureHessianSoap']:
            result_dim = 6
        elif getFunction in ['getVelocityGradient',
                             'getMagneticFieldGradient',
                             'getVectorPotentialGradient',
                             'getVelocityGradientSoap',
                             'getMagneticFieldGradientSoap',
                             'getVectorPotentialGradientSoap']:
            result_dim = 9
        elif getFunction in ['getVelocityHessian',
                             'getMagneticFieldHessian',
                             'getVectorPotentialHessian',
                             'getVelocityHessianSoap',
                             'getMagneticFieldHessianSoap',
                             'getVectorPotentialHessianSoap']:
            result_dim = 18
        else:
            print ('wrong result type requested in getData\n'
                 + 'maybe it\'s just missing from the list?')
            sys.exit()
            return None
        newshape = list(point_coords.shape[0:len(point_coords.shape)-1])
        newshape.append(result_dim)
        result_array = np.empty(newshape, dtype=np.float32)
        get_data(self.authToken,
                 ctypes.c_char_p(data_set),
                 ctypes.c_float(time),
                 ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                 point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                 result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array
    def getPosition(self,
            starttime, endtime, lag_dt,
            point_coords,
            sinterp = 4,
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
        self.lib.getPosition(
                self.authToken,
                ctypes.c_char_p(data_set),
                ctypes.c_float(starttime),
                ctypes.c_float(endtime),
                ctypes.c_float(lag_dt),
                ctypes.c_int(sinterp), ctypes.c_int(npoints),
                point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array
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
                    ctypes.c_char_p(data_set),
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
                    ctypes.c_char_p(data_set),
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
    def expand(self):
        repo_dir = os.path.dirname(inspect.getfile(libTDB))
        os.system('gcc -O3 -fPIC -Wall -c '
                + '-DCUTOUT_SUPPORT '
                + '-I' + self.srcdir + ' '
                + repo_dir + '/local_tools.c '
                + '-o ' + repo_dir + '/local_tools.o')
        linkcommand = ('gcc -dynamiclib '
                + self.libdir + '/stdsoap2.o '
                + self.libdir + '/soapC.o '
                + self.libdir + '/soapClient.o '
                + self.libdir + '/turblib.o '
                + repo_dir + '/local_tools.o '
                + '-o ' + repo_dir + '/libTDBe.so '
                + '-lhdf5 ')
        os.system(linkcommand)
        if self.connection_on:
            self.finalize()
            self.connection_on = True
        self.lib = np.ctypeslib.load_library('libTDBe', repo_dir)
        if self.connection_on:
            self.initialize()
        return None
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
            print ('wrong shape of initial condition in getBlineAlt, ', x0.shape)
            sys.exit()
            return None
        npoints = x0.shape[0]
        result_array = np.empty((nsteps+1, npoints, 3), dtype=np.float32)
        result_array[0] = x0
        self.lib.getBline(
                self.authToken,
                ctypes.c_char_p(data_set),
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
            print ('wrong shape of initial condition in getBlineAlt, ', x0.shape)
            sys.exit()
            return None
        nsteps = int(S / abs(ds))
        npoints = x0.shape[0]
        result_array = np.empty((npoints, nsteps+1, 3), dtype=np.float32)
        length_array = np.empty((npoints), dtype = np.int32)
        result_array[:, 0] = x0
        self.lib.getSphericalBoundedBline(
                self.authToken,
                ctypes.c_char_p(data_set),
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
                    ctypes.c_char_p(data_set),
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
            print ('wrong shape of initial condition in getBlineAlt, ', x0.shape)
            sys.exit()
            return None
        nsteps = int(S / abs(ds))
        npoints = x0.shape[0]
        result_array = np.empty((npoints, nsteps+1, 3), dtype=np.float32)
        length_array = np.empty((npoints), dtype = np.int32)
        result_array[:, 0] = x0
        self.lib.getRectangularBoundedBline(
                self.authToken,
                ctypes.c_char_p(data_set),
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
                    ctypes.c_char_p(data_set),
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

