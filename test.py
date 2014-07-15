import pyJHTDB
import pyJHTDB.cutout

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import h5py

def clean_2D_field(
        field_2D,
        dpi = 100,
        figname = 'tst',
        cmap = cm.jet,
        img_type = 'pdf'):
    fig = plt.figure(
                figsize=(field_2D.shape[1]*1./dpi,
                         field_2D.shape[0]*1./dpi))
    ax = fig.add_axes([.0, .0, 1., 1.], frameon=False)
    ax.set_axis_off()
    im = ax.imshow(field_2D,
            interpolation='none',
            cmap = cmap)
    fig.savefig(
            figname + '.' + img_type,
            dpi = dpi,
            format = img_type)
    return None

def main():
    pyJHTDB.cutout.get_big_cutout(
            cube_dim = 32,
            chunk_dim = 16,
            data_type = 'u',
            filename = 'tmp',
            base_website = 'test.turbulence.pha.jhu.edu')
    data = h5py.File('tmp.h5', mode = 'r')
    print data['u00000'].shape
    energy = (data['u00000'][0, :, :, 0]**2
            + data['u00000'][0, :, :, 1]**2
            + data['u00000'][0, :, :, 2]**2)
    clean_2D_field(energy, figname = 'tst_0yx')
    energy = (data['u00000'][:, 0, :, 0]**2
            + data['u00000'][:, 0, :, 1]**2
            + data['u00000'][:, 0, :, 2]**2)
    clean_2D_field(energy, figname = 'tst_z0x')
    energy = (data['u00000'][:, :, 0, 0]**2
            + data['u00000'][:, :, 0, 1]**2
            + data['u00000'][:, :, 0, 2]**2)
    clean_2D_field(energy, figname = 'tst_zy0')
    return None

if __name__ == '__main__':
    main()

