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
    if data_type in ['u', 'b', 'a']:
        ncomponents = 3
    elif data_type in ['p']:
        ncomponents = 1
    print 'Retrieving h5 file, size {0} MB = {1} MiB.'.format(
            xl*yl*zl*ncomponents * 4. / 10**6,
            xl*yl*zl*ncomponents * 4. / 2**20)
    urllib.urlretrieve(url, filename + '.h5')
    print 'Data downloaded and ' + filename + '.h5 written successfuly.'
    return None

def get_big_cutout(
        filename = 'tst',
        time = 0,
        cube_dim = 32,
        chunk_dim = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.jhu.pha.turbulence.testing-201302',
        base_website = 'turbulence.pha.jhu.edu'):
    if cube_dim % chunk_dim != 0:
        print 'in get_big_cutout, cube_dim must be a multiple of chunk_dim'
        return None
    if data_type in ['u', 'b', 'a']:
        ncomponents = 3
    elif data_type in ['p']:
        ncomponents = 1
    big_data_file = h5py.File(filename + '.h5', mode='w')
    big_cube = big_data_file.create_dataset(
            data_type + '{0:0>5}'.format(time*10),
            (cube_dim, cube_dim, cube_dim, ncomponents),
            np.float32,
            compression = 'lzf') ### is compression a good idea?
    for cz in range(cube_dim/chunk_dim):
        for cy in range(cube_dim/chunk_dim):
            for cx in range(cube_dim/chunk_dim):
                if not os.path.exists(filename + '_{0:0>2x}{1:0>2x}{2:0>2x}.h5'.format(cz, cy, cx)):
                    get_cutout(
                            filename + '_{0:0>2x}{1:0>2x}{2:0>2x}'.format(cz, cy, cx),
                            t0 = time, tl = 1,
                            x0 = cx*chunk_dim, y0 = cy*chunk_dim, z0 = cz*chunk_dim,
                            xl = chunk_dim, yl = chunk_dim, zl = chunk_dim,
                            data_set = data_set,
                            data_type = data_type,
                            auth_token = auth_token,
                            base_website = base_website)
                new_file = h5py.File(filename + '_{0:0>2x}{1:0>2x}{2:0>2x}.h5'.format(cz, cy, cx), mode='r')
                new_data = new_file[data_type + '{0:0>5}'.format(time*10)]
                big_cube[cz*chunk_dim:(cz+1)*chunk_dim,
                         cy*chunk_dim:(cy+1)*chunk_dim,
                         cx*chunk_dim:(cx+1)*chunk_dim, :] = new_data
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

