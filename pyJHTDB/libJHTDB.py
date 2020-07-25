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
import h5py

import pyJHTDB
from pyJHTDB.dbinfo import interpolation_code

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
        self.lib = np.ctypeslib.load_library(
            self.libname,
            os.path.abspath(os.path.join(lib_location, os.path.pardir)))
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

    def add_token(self,token):
        self.authToken = ctypes.c_char_p(token.encode('ascii'))

#    def add_hdf5_file(self, filename):
#        if pyJHTDB.found_h5py and (not filename in self.hdf5_file_list):
#            self.hdf5_file_list.append(filename)
#            data = pyJHTDB.h5py.File(filename + '.h5', mode = 'r')
#            self.hdf5_file_desc[filename] = {}
#            for key in ['_contents', '_dataset', '_size', '_start']:
#                self.hdf5_file_desc[filename][key] = data[key][:]
#            data.close()
#            return self.lib.turblibAddLocalSource(ctypes.c_char_p((filename + '.h5').encode('ascii')))
#        else:
#            return 0

    def getData(self,
                time, point_coords,
                sinterp=0, tinterp=0,
                data_set='isotropic1024coarse',
                getFunction='getVelocity',
                make_modulo=False):
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
        if (type(sinterp) == str):
            sinterp = interpolation_code[sinterp]
        if (type(tinterp) == str):
            tinterp = interpolation_code[tinterp]
        npoints = point_coords.shape[0]
        for i in range(1, len(point_coords.shape) - 1):
            npoints *= point_coords.shape[i]
        if make_modulo:
            pcoords = np.zeros(point_coords.shape, np.float64)
            pcoords[:] = point_coords
            np.mod(pcoords, 2 * np.pi, point_coords)
        if not getFunction[0:3] == 'get':
            getFunction = 'get' + getFunction
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
                           'getPressureGradient',
                           'getTemperatureGradient']:
            result_dim = 3
        elif getFunction in ['getVelocityAndPressure', 'getVelocityAndTemperature']:
            result_dim = 4
        elif getFunction in ['getPressureHessian', 'getTemperatureHessian']:
            result_dim = 6
        elif getFunction in ['getVelocityGradient',
                             'getMagneticFieldGradient',
                             'getVectorPotentialGradient']:
            result_dim = 9
        elif getFunction in ['getVelocityHessian',
                             'getMagneticFieldHessian',
                             'getVectorPotentialHessian']:
            result_dim = 18
        elif getFunction in ['getPressure', 'getTemperature']:
            result_dim = 1
        elif getFunction in ['getInvariant']:
            result_dim = 2
        else:
            print(('wrong result type requested in getData\n'
                   + 'maybe it\'s just missing from the list?'))
            sys.exit()
            return None
        newshape = list(point_coords.shape[0:len(point_coords.shape) - 1])
        newshape.append(result_dim)
        result_array = np.empty(newshape, dtype=np.float32)
        get_data(self.authToken,
                 ctypes.c_char_p(data_set.encode('ascii')),
                 ctypes.c_float(time),
                 ctypes.c_int(sinterp), ctypes.c_int(tinterp), ctypes.c_int(npoints),
                 point_coords.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                 result_array.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))))
        return result_array

    def getRawData(
            self,
            time=0,
            start=np.array([0, 0, 0], dtype=np.int),
            size=np.array([8, 8, 8], dtype=np.int),
            data_set='channel',
            getFunction='Velocity'):

        print(('This function is no longer supported. Please use getbigCutout instead.'))
        sys.exit()
        return None

        # if not self.connection_on:
        #     print('you didn\'t connect to the database')
        #     sys.exit()
        # if getFunction in ['Velocity',
        #                    'MagneticField',
        #                    'VectorPotential']:
        #     result_dim = 3
        # elif getFunction in ['Pressure', 'Temperature']:
        #     result_dim = 1
        # else:
        #     print(('wrong result type requested in getRawData\n'
        #            + 'maybe it\'s just missing from the list?'))
        #     sys.exit()
        #     return None
        # getFunction = 'getRaw' + getFunction
        # get_data = getattr(self.lib, getFunction)
        # result_array = np.empty(tuple(list(size[::-1]) + [result_dim]), dtype=np.float32)
        # get_data(self.authToken,
        #          ctypes.c_char_p(data_set.encode('ascii')),
        #          ctypes.c_int(time),
        #          ctypes.c_int(start[0]),
        #          ctypes.c_int(start[1]),
        #          ctypes.c_int(start[2]),
        #          ctypes.c_int(size[0]),
        #          ctypes.c_int(size[1]),
        #          ctypes.c_int(size[2]),
        #          result_array.ctypes.data_as(ctypes.POINTER(ctypes.c_char)))
        # return result_array
    
    def getCutout(
            self,
            data_set='isotropic1024coarse',
            field='u',
            time_step=int(0),
            start=np.array([1, 1, 1], dtype=np.int),
            end=np.array([8, 8, 8], dtype=np.int),
            step=np.array([1, 1, 1], dtype=np.int),
            filter_width=1):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()

        time_step=int(time_step)
        if field in ['u', 'a', 'b']:
            result_dim = 3
        elif field in ['p', 'd', 't']:
            result_dim = 1
        else:
            print(('wrong result type requested in getCutout\n'
                   + 'maybe it\'s just missing from the list?'))
            sys.exit()
            return None
        
        tempa=np.arange(start[0], end[0]+1, step[0])
        tempb=np.arange(start[1], end[1]+1, step[1])
        tempc=np.arange(start[2], end[2]+1, step[2])
        real_size=np.array([np.size(tempa), np.size(tempb), np.size(tempc)], dtype=np.int)
        
        getFunction = 'getCutout'
        get_data = getattr(self.lib, getFunction)
        result_array = np.empty(tuple(list(real_size[::-1]) + [result_dim]), dtype=np.float32)
        try:
            get_data(self.authToken,
                 ctypes.c_char_p(data_set.encode('ascii')),
                 ctypes.c_char_p(field.encode('ascii')),
                 ctypes.c_int(time_step),
                 ctypes.c_int(start[0]),
                 ctypes.c_int(start[1]),
                 ctypes.c_int(start[2]),
                 ctypes.c_int(end[0]),
                 ctypes.c_int(end[1]),
                 ctypes.c_int(end[2]),
                 ctypes.c_int(step[0]),
                 ctypes.c_int(step[1]),
                 ctypes.c_int(step[2]),
                 ctypes.c_int(filter_width),
                 result_array.ctypes.data_as(ctypes.POINTER(ctypes.c_char)))
            
        except Exception as es:
            print(es)
            raise

        return result_array

    def getbigCutout(
            self,
            data_set='isotropic1024coarse',
            fields='u',
            t_start=int(1),
            t_end=int(1),
            t_step=int(1),
            start=np.array([1, 1, 1], dtype=np.int),
            end=np.array([8, 8, 8], dtype=np.int),
            step=np.array([1, 1, 1], dtype=np.int),
            filter_width=1,
            filename='N/A'):
            #hdf5_output=True):
        if not self.connection_on:
            print('you didn\'t connect to the database')
            sys.exit()

        if (filename.lower()=='n/a' or filename.lower()=='na'):
            hdf5_output=False
        else:
            hdf5_output=True

        idx_t=np.arange(t_start, t_end+1, t_step)
        idx_x=np.arange(start[0], end[0]+1, step[0])
        idx_y=np.arange(start[1], end[1]+1, step[1])
        idx_z=np.arange(start[2], end[2]+1, step[2])
        nnt=np.size(idx_t)
        nnx=np.size(idx_x)
        nny=np.size(idx_y)
        nnz=np.size(idx_z)

        npoints=nnx*nny*nnz
        tem=0

        for field in fields:
            if field == 'u':
                tem = tem + 3
            elif field == 'a':
                tem = tem + 3
            elif field == 'b':
                tem = tem + 3
            elif field == 'p':
                tem = tem + 1
            elif field == 'd':
                tem = tem + 1
            elif field == 't':
                tem = tem + 1 
            else:
                print(('wrong field type requested in getCutout\n'
                    + 'maybe it\'s just missing from the list?'))
                sys.exit()
                return None

        if (npoints*nnt*tem>(1024**3)*4): #a full snapshot of 1024^3 with u and p
            print(('The file size would exceed our limit 16GB. Please reduce the file size.'))
            sys.exit()
            return None

        if (hdf5_output):
            nl = '\r\n'
            hdf5_file, xdmf_file, shape=self.hdf5_init(filename, data_set,t_start,t_end,t_step,start,end,step,filter_width,idx_x,idx_y,idx_z)

        for field in fields:
            if field == 'u':
                VarName="Velocity"
                dim = 3
            elif field == 'a':
                VarName="VectorPotential"
                dim = 3
            elif field == 'b':
                VarName="MagneticField"
                dim = 3
            elif field == 'p':
                VarName="Pressure"
                dim = 1
            elif field == 'd':
                VarName="Density"
                dim = 1
            elif field == 't':
                VarName="Temperature"
                dim = 1 
            else:
                print(('wrong field type requested in getCutout\n'
                    + 'maybe it\'s just missing from the list?'))
                sys.exit()
                return None

            split_no=int(np.ceil(npoints/(192000000/dim)))
            tmp=np.array_split(np.arange(npoints).reshape(nnx,nny,nnz), split_no)
            
            if (hdf5_output):
                    print(f"    <Grid Name=\"{VarName}\" GridType=\"Collection\" CollectionType=\"Temporal\">{nl}", file=xdmf_file)

            for time_step in np.arange(t_start, t_end+1, t_step):

                result=np.zeros((nnz,nny,nnx,dim),dtype='float32')

                for t in range(split_no):
                    xyzs0 = np.unravel_index(tmp[t][0,0,0], (nnx,nny,nnz))
                    xyze0 = np.unravel_index(tmp[t][-1,-1,-1], (nnx,nny,nnz))
                    xyzs1 = (idx_x[xyzs0[0]], idx_y[xyzs0[1]], idx_z[xyzs0[2]])
                    xyze1 = (idx_x[xyze0[0]], idx_y[xyze0[1]], idx_z[xyze0[2]])
            
                    temp = self.getCutout(
                        data_set=data_set, field=field, time_step=time_step,
                        start=np.array(xyzs1, dtype = np.int),
                        end=np.array(xyze1, dtype = np.int),
                        step=np.array(step, dtype = np.int),
                        filter_width=filter_width)
            
                    result[xyzs0[2]:xyze0[2]+1, xyzs0[1]:xyze0[1]+1, xyzs0[0]:xyze0[0]+1,:] = temp

                if (hdf5_output):
                    self.hdf5_writing(filename,result,data_set,VarName,dim,time_step,hdf5_file,xdmf_file,shape)
                    
            if (hdf5_output):
                    print(f"    </Grid>{nl}", file=xdmf_file)

        if (hdf5_output):
            self.hdf5_end(hdf5_file,xdmf_file)

        return result

    def hdf5_init(
            self,
            filename,
            data_set,
            t_start,
            t_end,
            t_step,
            start,
            end,
            step,
            filter_width,
            idx_x,idx_y,idx_z):

        idx_x=idx_x-1
        idx_y=idx_y-1
        idx_z=idx_z-1

        if data_set in ["channel","channel5200", "transition_bl"]:
                
            if data_set == "channel":
                ygrid = pyJHTDB.dbinfo.channel['ynodes']
                dx=8.0*np.pi/2048
                dz=3.0*np.pi/1536
                x_offset=0

            elif data_set == "channel5200":
                ygrid = pyJHTDB.dbinfo.channel5200['ynodes']
                dx=8.0*np.pi/10240.0
                dz=3.0*np.pi/7680.0
                x_offset=0

            elif data_set == "transition_bl":
                ygrid = pyJHTDB.dbinfo.transition_bl['ynodes']
                dx=0.292210466240511
                dz=0.117244748412311
                x_offset=30.218496172581567
                        
            xcoor=idx_x*dx+x_offset
            ycoor=ygrid[idx_y]
            zcoor=idx_z*dz
        else:
            if data_set in ["isotropic1024coarse", "isotropic1024fine", "mhd1024", "mixing"]:
                dx=2.0*np.pi/1024.0
            elif data_set in ["isotropic4096", "rotstrat4096"]:
                dx=2.0*np.pi/4096.0
            elif data_set in ["isotropic8192"]:
                dx=2.0*np.pi/8192.0
                
            xcoor=idx_x*dx
            ycoor=idx_y*dx
            zcoor=idx_z*dx

        #filename=data_set
        fh = h5py.File(filename+'.h5', 'x', driver='core', block_size=16, backing_store=True)
        fh.attrs["dataset"] = np.string_(data_set)
        #fh.attrs["timeStep"] = time_step
        fh.attrs["t_start"] = t_start
        fh.attrs["t_end"] = t_end
        fh.attrs["t_step"] = t_step
        fh.attrs["x_start"] = start[0]
        fh.attrs["y_start"] = start[1]
        fh.attrs["z_start"] = start[2]
        fh.attrs["x_end"] = end[0]
        fh.attrs["y_end"] = end[1]
        fh.attrs["z_end"] = end[2]
        fh.attrs["x_step"] = step[0]
        fh.attrs["y_step"] = step[1]
        fh.attrs["z_step"] = step[2]
        fh.attrs["filterWidth"] = filter_width

        shape = [0]*3
        shape[0] = np.size(idx_z)
        shape[1] = np.size(idx_y)
        shape[2] = np.size(idx_x)

        dset = fh.create_dataset("xcoor", (shape[2],), maxshape=(shape[2],))
        dset[...]=xcoor
        dset = fh.create_dataset("ycoor", (shape[1],), maxshape=(shape[1],))
        dset[...]=ycoor
        dset = fh.create_dataset("zcoor", (shape[0],), maxshape=(shape[0],))
        dset[...]=zcoor

        nl = '\r\n'

        tf=open(filename+".xmf", "w")
        print(f"<?xml version=\"1.0\" ?>{nl}", file=tf)
        print(f"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>{nl}", file=tf)
        print(f"<Xdmf Version=\"2.0\">{nl}", file=tf)
        print(f"  <Domain>{nl}", file=tf)

        return fh, tf, shape

    def hdf5_writing(
            self,
            filename,
            result,
            data_set,
            VarName,
            dim,
            time_step,
            fh,
            tf,
            shape):

        H5_ds_name='{0}_{1:04d}'.format(VarName,time_step)
        dset = fh.create_dataset(H5_ds_name,
            (shape[0], shape[1], shape[2], dim), maxshape=(shape[0], shape[1], shape[2], dim))
        dset[...]=result

        if (dim==3):
            Attribute_Type="Vector"
        elif (dim==1):
            Attribute_Type="Scalar"
                
        #filename=data_set
        nl = '\r\n'
        
        print(f"      <Grid Name=\"Structured Grid\" GridType=\"Uniform\">{nl}", file=tf)
        print(f"        <Time Value=\"{time_step}\" />{nl}", file=tf)
        print(f"        <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"{shape[0]} {shape[1]} {shape[2]}\"/>{nl}", file=tf)
        print(f"        <Geometry GeometryType=\"VXVYVZ\">{nl}", file=tf)
        print(f"          <DataItem Name=\"Xcoor\" Dimensions=\"{shape[2]}\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">{nl}", file=tf)
        print(f"            {filename}.h5:/xcoor{nl}", file=tf)
        print(f"          </DataItem>{nl}", file=tf)
        print(f"          <DataItem Name=\"Ycoor\" Dimensions=\"{shape[1]}\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">{nl}", file=tf)
        print(f"            {filename}.h5:/ycoor{nl}", file=tf)
        print(f"          </DataItem>{nl}", file=tf)
        print(f"          <DataItem Name=\"Zcoor\" Dimensions=\"{shape[0]}\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">{nl}", file=tf)
        print(f"            {filename}.h5:/zcoor{nl}", file=tf)
        print(f"          </DataItem>{nl}", file=tf)
        print(f"        </Geometry>{nl}", file=tf)

        print(f"{nl}", file=tf)

        print(f"        <Attribute Name=\"{VarName}\" AttributeType=\"{Attribute_Type}\" Center=\"Node\">{nl}", file=tf)
        print(f"          <DataItem Dimensions=\"{shape[0]} {shape[1]} {shape[2]} {dim}\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">{nl}", file=tf)
        print(f"            {filename}.h5:/{H5_ds_name}{nl}", file=tf)
        print(f"          </DataItem>{nl}", file=tf)
        print(f"        </Attribute>{nl}", file=tf)

        print(f"      </Grid>{nl}", file=tf)

        print(f"{nl}", file=tf)

        return

    def hdf5_end(
            self,
            fh,
            tf):

        fh.close()

        nl = '\r\n'

        print(f"  </Domain>{nl}", file=tf)
        print(f"</Xdmf>{nl}", file=tf)
        tf.close()

        return

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

    def getBoxFilterGradient(self,
            time, point_coords,
            data_set = 'isotropic1024coarse',
            make_modulo = False,
            field = 'velocity',
            filter_width = 7*2*np.pi / 1024,
            FD_spacing = 4*2*np.pi / 1024):
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
        field=field.lower()
        if field in ['velocity', 'magnetic', 'vectorPotential']:
            result_dim = 9
        elif field in ['pressure', 'temperature']:
            result_dim = 3
        result_array = np.empty([npoints,result_dim], dtype=np.float32)
        self.lib.getBoxFilterGradient(self.authToken,
                 ctypes.c_char_p(data_set.encode('ascii')),
                 ctypes.c_char_p(field.encode('ascii')),
                 ctypes.c_float(time),
                 ctypes.c_float(filter_width),
                 ctypes.c_float(FD_spacing),
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
            x_start = 1, y_start = 1, z_start = 1,
            x_end = 4, y_end = 4, z_end = 4,
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
                ctypes.c_int32(x_start),
                ctypes.c_int32(y_start),
                ctypes.c_int32(z_start),
                ctypes.c_int32(x_end),
                ctypes.c_int32(y_end),
                ctypes.c_int32(z_end),
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
            print(('wrong shape of initial condition in getBlineSphereBounded, ', x0.shape))
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
            print(('wrong shape of initial condition in getBlineSphereBoundedDebug, ', x0.shape))
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
            print(('wrong shape of initial condition in getBlineRectBounded, ', x0.shape))
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

