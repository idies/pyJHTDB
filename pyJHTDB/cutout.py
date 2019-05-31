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
if sys.version_info[0] == 2:
    import urllib
elif sys.version_info[0] == 3:
    import urllib.request, urllib.parse, urllib.error

import numpy as np
import h5py

def get_cutout(
        filename = 'tst',
        t0 = 0, tl = 1,
        x0 = 0, xl = 16,
        y0 = 0, yl = 16,
        z0 = 0, zl = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.jhu.pha.turbulence.testing-201302',
        base_website = 'dsp033.pha.jhu.edu',
        file_format = 'hdf5'):
    url = ('http://{0}/jhtdb/getcutout/{1}/{2}/{3}/'.format(
                base_website, auth_token, data_set, data_type)
         + '{0},{1}/'.format(t0, tl)
         + '{0},{1}/'.format(x0, xl)
         + '{0},{1}/'.format(y0, yl)
         + '{0},{1}/'.format(z0, zl)
         + '{0}/'.format(file_format))
    print(url)
    if data_type in ['u', 'b', 'a']:
        ncomponents = 3
    elif data_type in ['p']:
        ncomponents = 1
    elif data_type in ['ub']:
        ncomponents = 6
    print('Retrieving h5 file, size {0} MB = {1} MiB.'.format(
            xl*yl*zl*ncomponents * 4. / 10**6,
            xl*yl*zl*ncomponents * 4. / 2**20))
    if sys.version_info[0] == 2:
        urllib.urlretrieve(url, filename + '.h5')
    elif sys.version_info[0] == 3:
        urllib.request.urlretrieve(url, filename + '.h5')
    # check if file downloaded ok
    data = h5py.File(filename + '.h5', mode = 'r')
    data.close()
    print('Data downloaded and ' + filename + '.h5 written successfuly.')
    return None

def get_big_cutout(
        filename = 'tst',
        t0 = 0, tl = 1,
        x0 = 0, xl = 32,
        y0 = 0, yl = 32,
        z0 = 0, zl = 32,
        chunk_xdim = 16,
        chunk_ydim = 16,
        chunk_zdim = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.jhu.pha.turbulence.testing-201302',
        base_website = 'dsp033.pha.jhu.edu',
        file_format = 'hdf5'):
    big_data_file = h5py.File(filename + '.h5', mode='w')
    xchunk_list = [chunk_xdim for n in range(int(xl / chunk_xdim))]
    if not (xl % chunk_xdim == 0):
        xchunk_list.append(xl % chunk_xdim)
    ychunk_list = [chunk_ydim for n in range(int(yl / chunk_ydim))]
    if not (yl % chunk_ydim == 0):
        ychunk_list.append(yl % chunk_ydim)
    zchunk_list = [chunk_zdim for n in range(int(zl / chunk_zdim))]
    if not (zl % chunk_zdim == 0):
        zchunk_list.append(zl % chunk_zdim)
    for current_data_type in data_type:
        if current_data_type in ['u', 'b', 'a']:
            ncomponents = 3
        elif current_data_type in ['p']:
            ncomponents = 1
        big_data = []
        for time in range(t0, t0 + tl):
            big_data.append(big_data_file.create_dataset(
                    current_data_type + '{0:0>5}'.format(time*10),
                    (zl, yl, xl, ncomponents),
                    np.float32,
                    compression = 'lzf')) ### is compression a good idea?
        for cz in range(len(zchunk_list)):
            for cy in range(len(ychunk_list)):
                for cx in range(len(xchunk_list)):
                    for time in range(t0, t0 + tl):
                        tmp_filename = (filename
                                      + '_{0:0>2x}{1:0>2x}{2:0>2x}_{3}'.format(cz, cy, cx, current_data_type))
                        if not os.path.exists(tmp_filename + '.h5'):
                            get_cutout(
                                    tmp_filename,
                                    t0 = time, tl = tl,
                                    x0 = x0+cx*chunk_xdim, y0 = y0+cy*chunk_ydim, z0 = z0+cz*chunk_zdim,
                                    xl = xchunk_list[cx], yl = ychunk_list[cy], zl = zchunk_list[cz],
                                    data_set = data_set,
                                    data_type = current_data_type,
                                    auth_token = auth_token,
                                    base_website = base_website,
                                    file_format = file_format)
                        new_file = h5py.File(tmp_filename + '.h5', mode='r')
                        new_data = new_file[current_data_type + '{0:0>5}'.format(time*10)]
                        big_data[time - t0][cz*chunk_zdim:cz*chunk_zdim+zchunk_list[cz],
                                            cy*chunk_ydim:cy*chunk_ydim+ychunk_list[cy],
                                            cx*chunk_xdim:cx*chunk_xdim+xchunk_list[cx], :] = new_data
    big_data_file.create_dataset(
            '_contents',
            new_file['_contents'].shape,
            new_file['_contents'].dtype)
    big_data_file['_contents'][:] = new_file['_contents'][:]
    if data_type == 'ub':
        big_data_file['_contents'][0] = 5
    big_data_file.create_dataset(
            '_dataset',
            new_file['_dataset'].shape,
            new_file['_dataset'].dtype)
    big_data_file['_dataset'][:] = new_file['_dataset'][:]
    big_data_file.create_dataset(
            '_size',
            new_file['_size'].shape,
            new_file['_size'].dtype)
    big_data_file['_size'][0] = new_file['_size'][0]
    big_data_file['_size'][1] = xl
    big_data_file['_size'][2] = yl
    big_data_file['_size'][3] = zl
    big_data_file.create_dataset(
            '_start',
            new_file['_start'].shape,
            new_file['_start'].dtype)
    big_data_file['_start'][0] = new_file['_start'][0]
    big_data_file['_start'][1] = x0
    big_data_file['_start'][2] = y0
    big_data_file['_start'][3] = z0
    new_file.close()
    big_data_file.close()
    return None

def main():
    # dumb test
    for data_type in ['u', 'b', 'ub']:
        get_cutout(
            data_set = 'mhd1024',
            data_type = data_type)
        f0 = h5py.File('tst.h5', mode='r')
        print((data_type, f0['_contents'][:]))
        f0.close()
    return None

if __name__ == '__main__':
    main()

