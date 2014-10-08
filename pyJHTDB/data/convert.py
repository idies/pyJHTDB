import h5py
import numpy as np

def main():
    cgrid = h5py.File('channel_grid.h5', mode = 'r')

    #using single precision since all getF functions need single precision
    np.save('channel_xgrid', cgrid['x'][:].astype(np.float32))
    np.save('channel_ygrid', cgrid['y'][:].astype(np.float32))
    np.save('channel_zgrid', cgrid['z'][:].astype(np.float32))
    return None

if __name__ == '__main__':
    main()

