import os
import urllib
import numpy as np
import h5py

def get_cutout(
        filename = 'tst',
        t0 = 0, tl = 1,
        x0 = 0, xl = 32,
        y0 = 0, yl = 32,
        z0 = 0, zl = 32,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.jhu.pha.turbulence.testing-201302',
        base_website = 'turbulence.pha.jhu.edu'):
    url = ('http://' + base_website + '/download.aspx/'
         + auth_token + '/'
         + data_set + '/' + data_type + '/'
         + '{0},{1}/'.format(t0, tl)
         + '{0},{1}/'.format(x0, xl)
         + '{0},{1}/'.format(y0, yl)
         + '{0},{1}/'.format(z0, zl))
    print url
    if data_type in ['u', 'b', 'a']:
        ncomponents = 3
    elif data_type in ['p']:
        ncomponents = 1
    elif data_type in ['ub']:
        ncomponents = 6
    print 'Retrieving h5 file, size {0} MB = {1} MiB.'.format(
            xl*yl*zl*ncomponents * 4. / 10**6,
            xl*yl*zl*ncomponents * 4. / 2**20)
    urllib.urlretrieve(url, filename + '.h5')
    # check if file downloaded ok
    data = h5py.File(filename + '.h5', mode = 'r')
    data.close()
    print 'Data downloaded and ' + filename + '.h5 written successfuly.'
    return None

def get_big_cutout(
        filename = 'tst',
        t0 = 0, tl = 1,
        x0 = 0, xl = 32,
        y0 = 0, yl = 32,
        z0 = 0, zl = 32,
        chunk_dim = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.jhu.pha.turbulence.testing-201302',
        base_website = 'turbulence.pha.jhu.edu'):
    if ((xl % chunk_dim != 0) or
        (yl % chunk_dim != 0) or
        (zl % chunk_dim != 0)):
        print 'in get_big_cutout, each dimension except time must be a multiple of chunk_dim'
        return None
    big_data_file = h5py.File(filename + '.h5', mode='w')
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
        for cz in range(zl/chunk_dim):
            for cy in range(yl/chunk_dim):
                for cx in range(xl/chunk_dim):
                    for time in range(t0, t0 + tl):
                        tmp_filename = filename + '_{0:0>2x}{1:0>2x}{2:0>2x}_{3}'.format(cz, cy, cx, current_data_type)
                        if not os.path.exists(tmp_filename + '.h5'):
                            get_cutout(
                                    tmp_filename,
                                    t0 = time, tl = tl,
                                    x0 = cx*chunk_dim, y0 = cy*chunk_dim, z0 = cz*chunk_dim,
                                    xl = chunk_dim, yl = chunk_dim, zl = chunk_dim,
                                    data_set = data_set,
                                    data_type = current_data_type,
                                    auth_token = auth_token,
                                    base_website = base_website)
                        new_file = h5py.File(tmp_filename + '.h5', mode='r')
                        new_data = new_file[current_data_type + '{0:0>5}'.format(time*10)]
                        big_data[time - t0][cz*chunk_dim:(cz+1)*chunk_dim,
                                            cy*chunk_dim:(cy+1)*chunk_dim,
                                            cx*chunk_dim:(cx+1)*chunk_dim, :] = new_data
    big_data_file.create_dataset(
            '_contents',
            new_file['_contents'].shape,
            new_file['_contents'].dtype)
    big_data_file['_contents'][:] = new_file['_contents'][:]
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
    get_cutout(filename = 'tst0')
    get_big_cutout(filename = 'tst1')
    f0 = h5py.File('tst0.h5', mode='r')
    f1 = h5py.File('tst1.h5', mode='r')
    print np.abs(f0['u00000'][:] - f1['u00000'][:]).max()
    return None

if __name__ == '__main__':
    main()

