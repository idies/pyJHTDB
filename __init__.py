"""
    Python tools and wrappers for the Johns Hopkins Turbulence Databases C library.
"""

import sys
import numpy
import ctypes

class libTDB:
    def __init__(self,
            libname = 'libTDB',
            libdir = '.',
            auth_token = 'edu.jhu.pha.turbulence.testing-201302'):
        self.lib = numpy.ctypeslib.load_library(libname, libdir)
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
        if not (point_coords.dtype == numpy.float32):
            print 'point coordinates in getData must be floats. stopping.'
            sys.exit()
            return None
        npoints = point_coords.shape[0]
        for i in range(1, len(point_coords.shape)-1):
            npoints *= point_coords.shape[i]
        if make_modulo:
            pcoords = numpy.zeros(point_coords.shape, numpy.float64)
            pcoords[:] = point_coords
            numpy.mod(pcoords, 2*numpy.pi, point_coords)
        get_data = getattr(self.lib, getFunction)
        if getFunction in ['getVelocity',
                           'getForce',
                           'getMagneticField',
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
        result_array = numpy.empty(newshape, dtype=numpy.float32)
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
        result_array = numpy.empty((npoints, 3), dtype=numpy.float32)
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
            points = numpy.array([[0, 0, 0]]).astype(numpy.float32),
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
        history = numpy.zeros(
                (iterations+1, npoints, 3),
                dtype = numpy.float32)
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
            bsize = numpy.maximum(numpy.sqrt(numpy.sum(result_array**2, axis = 1)), 1e-18)
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
            point = numpy.array([0, 0, 0]).astype(numpy.float32),
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
        history = numpy.zeros(
                (iterations+1, 3),
                dtype = numpy.float32)
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
            bsize = max(numpy.sqrt(numpy.sum(result_array**2)), 1e-18)
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

