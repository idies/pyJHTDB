import numpy as np

isotropic1024coarse = {'name'   : 'isotropic1024coarse',
                       'nx'     : 1024,
                       'ny'     : 1024,
                       'nz'     : 1024,
                       'lx'     : 2*np.pi,
                       'ly'     : 2*np.pi,
                       'lz'     : 2*np.pi,
                       'dx'     : np.pi/512,
                       'dy'     : np.pi/512,
                       'dz'     : np.pi/512,
                       'xnodes' : (np.pi/512)*np.array(range(1024)),
                       'ynodes' : (np.pi/512)*np.array(range(1024)),
                       'znodes' : (np.pi/512)*np.array(range(1024))}

channel = {'name'   : 'channel',
           'nx'     : 2048,
           'ny'     : 512,
           'nz'     : 1536,
           'lx'     : 8*np.pi,
           'ly'     : 2.0,
           'lz'     : 6*np.pi,
           'dx'     : np.pi/256,
           'dy'     :    1./256,                         ## WRONG, FIXME
           'dz'     : np.pi/256,
           'xnodes' : (np.pi/256)*np.array(range(2048)),
           'ynodes' : (   1./256)*np.array(range( 512)), ## WRONG, FIXME
           'znodes' : (np.pi/256)*np.array(range(1536))}

