#! /usr/bin/env python

import os
import sys
import ctypes
import math

homefolder = os.path.expanduser('~')

try:
    import numpy
except ImportError:
    print 'Numpy is needed for this to work.'
    print 'You should be able to find installation instructions at http://www.scipy.org'
    sys.exit()
try:
    import matplotlib.pyplot as plt
except ImportError:
    print 'matplotlib is needed for contour plots.'
    print 'You should be able to find installation instructions at http://matplotlib.sourceforge.net'
    plt = None

import pyJHTDB
import pyJHTDB.dbinfo

def turbc_clone(N=10):
    #time = 0.364
    #turbc.c has a fixed time, but it makes more sense to have a random one
    time = 0.002 * numpy.random.randint(1024)
    # points must be created with the single precision data type
    points = numpy.empty((N, 3), dtype = 'float32')
    # [:,:] is there to force the conversion from the return type of random_sample to single precision
    points[:,:] = 2*math.pi*numpy.random.random_sample(size = (N, 3))[:,:]
    spatialInterp  = 6  # 6 point Lagrange
    temporalInterp = 0  # no time interpolation
    FD4Lag4        = 40 # 4 point Lagrange interp for derivatives

    # mhdc has starttime .364 and endtime .376
    startTime = 0.002 * numpy.random.randint(1024)
    endTime = startTime + 0.012
    lag_dt = 0.0004

    print 'Coordinates of {0} points where variables are requested:'.format(N)
    for p in range(N):
        print p, points[p]
    print 'Data is requested at time {0}'.format(time)

    # load shared library
    lTDB = pyJHTDB.libTDB(libname = 'libTDB',
                          srcdir = homefolder + '/ext_installs/turblib')
    #initialize webservices
    lTDB.initialize()

    print 'Requesting velocity at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getVelocitySoap')
    for p in range(N):
        print p, result[p]
    print 'Requesting forcing at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getForce')
    for p in range(N):
        print p, result[p]
    print 'Requesting velocity and pressure at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getVelocityAndPressureSoap')
    for p in range(N):
        print p, result[p][0:2], 'p = {0}'.format(result[p][3])
    print 'Requesting velocity gradient at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityGradientSoap')
    for p in range(N):
        print ('{0}: duxdx = {1}, duxdy = {2}, duxdz = {3}, '.format(p, result[p][0], result[p][1], result[p][2])
                  + 'duydx = {0}, duydy = {1}, duydz = {2}, '.format(result[p][3], result[p][4], result[p][5])
                  + 'duzdx = {0}, duzdy = {1}, duzdz = {2}'.format(result[p][6], result[p][7], result[p][8]))
    print 'Requesting velocity hessian at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityHessianSoap')
    for p in range(N):
        print ('{0}: '.format(p)
            + 'd2uxdxdx = {0}, d2uxdxdy = {1}, d2uxdxdz = {2}, '.format(result[p][ 0], result[p][ 1], result[p][ 2])
            + 'd2uxdydy = {0}, d2uxdydz = {1}, d2uxdzdz = {2}, '.format(result[p][ 3], result[p][ 4], result[p][ 5])
            + 'd2uydxdx = {0}, d2uydxdy = {1}, d2uydxdz = {2}, '.format(result[p][ 6], result[p][ 7], result[p][ 8])
            + 'd2uydydy = {0}, d2uydydz = {1}, d2uydzdz = {2}, '.format(result[p][ 9], result[p][10], result[p][11])
            + 'd2uzdxdx = {0}, d2uzdxdy = {1}, d2uzdxdz = {2}, '.format(result[p][12], result[p][13], result[p][14])
            + 'd2uzdydy = {0}, d2uzdydz = {1}, d2uzdzdz = {2}, '.format(result[p][15], result[p][16], result[p][17]))
    print 'Requesting velocity laplacian at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityLaplacianSoap')
    for p in range(N):
        print ('{0}: '.format(p)
            + 'grad2ux = {0}, grad2uy = {1}, grad2uz = {2}, '.format(result[p][0], result[p][1], result[p][2]))
    print 'Requesting pressure gradient at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getPressureGradientSoap')
    for p in range(N):
        print ('{0}: '.format(p)
            + 'dpdx = {0}, dpdy = {1}, dpdz = {2}, '.format(result[p][0], result[p][1], result[p][2]))
    print 'Requesting pressure hessian at {0} points...'.format(N)
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityHessianSoap')
    for p in range(N):
        print ('{0}: '.format(p)
            + 'd2pdxdx = {0}, d2pdxdy = {1}, d2pdxdz = {2}, '.format(result[p][0], result[p][1], result[p][2])
            + 'd2pdydy = {0}, d2pdydz = {1}, d2pdzdz = {2}, '.format(result[p][3], result[p][4], result[p][5]))

#    print 'Requesting position at {0} points, starting at time {1} and ending at time {2}...'.format(N, startTime, endTime)
#    result = pyJHTDB.getPosition(startTime, endTime, lag_dt, points, sinterp = spatialInterp)
#    print 'Coordinates of {0} points at startTime:'.format(N)
#    for p in range(N):
#        print p, points[p]
#    print 'Coordinates of {0} points at endTime:'.format(N)
#    for p in range(N):
#        print p, result[p]

    ##  only if matplotlib is present (i.e. plt is not None)
    if plt:
        ken_contours(
                'kin_en_contours',
                lTDB,
                spatialInterp = spatialInterp,
                temporalInterp = temporalInterp,
                time = 0.002 * numpy.random.randint(1024),
                spacing = 2 * math.pi * 2.**(-10),
                nx = 64, ny = 64,
                xoff = 2*math.pi * numpy.random.rand(),
                yoff = 2*math.pi * numpy.random.rand(),
                zoff = 2*math.pi * numpy.random.rand())

    #finalize webservices
    lTDB.finalize()
    return None

def ken_contours(
        figname,
        lTDB,
        levels = 30,
        spatialInterp = 4, temporalInterp = 0,
        time = 0.002 * 512,
        spacing = math.pi * 2.**(-9),
        nx = 64, ny = 64,
        xoff = .0, yoff = .0, zoff = .0):
    """
        Generate a simple contour plot
        see http://matplotlib.sourceforge.net/examples/pylab_examples/contour_demo.html
        for information on how to make prettier plots.

        This function assumes the webservices have already been initialized,
        so call pyJHTDB.init() before calling it, and pyJHTDB.finalize() afterwards
    """
    x = spacing * numpy.arange(0, nx, 1, dtype = 'float32') + xoff
    y = spacing * numpy.arange(0, ny, 1, dtype = 'float32') + yoff
    points = numpy.empty((nx, ny, 3), dtype = 'float32')
    points[:, :, 0] = x[:, numpy.newaxis]
    points[:, :, 1] = y[numpy.newaxis, :]
    points[:, :, 2] = zoff

    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getVelocitySoap')

    energy = .5*(numpy.sqrt(result[:,:,0]**2 + result[:,:,1]**2 + result[:,:,2]**2)).transpose()
    fig = plt.figure(figsize=(6.,6.))
    ax = fig.add_axes([.0, .0, 1., 1.])
    contour = ax.contour(x, y, energy, levels)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.clabel(contour, inline=1, fontsize=10)
    plt.title('Energy contours, t = {0:.5}, z = {1:.3}'.format(time, zoff))
    fig.savefig(figname + '.eps', format = 'eps', bbox_inches = 'tight')
    return None

def spectra_check(
        info = pyJHTDB.dbinfo.isotropic1024coarse,
        libTDB = None):
    if not libTDB:
        print ('no library given')
        return None
    nlines = 32
    print (nlines, info['nx'], 3)
    print (nlines, info['ny'], 3)
    print (nlines, info['nz'], 3)
    lines  = [numpy.zeros((nlines, info['nx'], 3), dtype = numpy.float32),
              numpy.zeros((nlines, info['ny'], 3), dtype = numpy.float32),
              numpy.zeros((nlines, info['nz'], 3), dtype = numpy.float32)]
    lines[0][:, :, 0] = info['dx']*numpy.arange(.0, info['nx'], 1)
    lines[0][:, :, 1] = info['dy']*numpy.random.randint(0, info['ny'], size = (nlines))[:, numpy.newaxis]
    lines[0][:, :, 2] = info['dz']*numpy.random.randint(0, info['nz'], size = (nlines))[:, numpy.newaxis]
    lines[1][:, :, 0] = info['dx']*numpy.random.randint(0, info['nx'], size = (nlines))[:, numpy.newaxis]
    lines[1][:, :, 1] = info['dy']*numpy.arange(.0, info['ny'], 1)
    lines[1][:, :, 2] = info['dz']*numpy.random.randint(0, info['nz'], size = (nlines))[:, numpy.newaxis]
    lines[2][:, :, 0] = info['dx']*numpy.random.randint(0, info['nx'], size = (nlines))[:, numpy.newaxis]
    lines[2][:, :, 1] = info['dy']*numpy.random.randint(0, info['ny'], size = (nlines))[:, numpy.newaxis]
    lines[2][:, :, 2] = info['dz']*numpy.arange(.0, info['nz'], 1)
    fig = plt.figure(figsize=(12.,6.))
    axu = fig.add_axes([.05, .1, .4, .8])
    axp = fig.add_axes([.55, .1, .4, .8])
    coordname = ['x', 'y', 'z']
    for i in range(3):
        result = libTDB.getData(.0, lines[i],
            sinterp = 0, tinterp = 0,
            data_set = info['name'], getFunction = 'getVelocityAndPressureSoap')
        spec = numpy.fft.rfft(
                numpy.sum(result[:, :, :3]**2, axis = 2), axis = 1)
        axu.plot(numpy.average(numpy.abs(spec), axis = 0), label = '$i = ' + coordname[i] + '$')
        spec = numpy.fft.rfft(
                result[:, :, 3]**2, axis = 1)
        axp.plot(numpy.average(numpy.abs(spec), axis = 0), label = '$i = ' + coordname[i] + '$')
    axu.set_ylabel('$\\langle u_x^2 + u_y^2 + u_z^2 \\rangle$')
    axp.set_ylabel('$\\langle p^2 \\rangle$')
    for ax in [axu, axp]:
        ax.set_xlabel('$k_i$')
        ax.legend(loc = 'best')
        ax.set_xscale('log')
        ax.set_yscale('log')
    fig.savefig('spec.pdf', format = 'pdf')
    return None

def contour_check(
        info = pyJHTDB.dbinfo.isotropic1024coarse,
        libTDB = None):
    if not libTDB:
        print ('no library given')
        return None
    nplanes = 1
    planes  = [numpy.zeros((nplanes, info['ny'], info['nz'], 3), dtype = numpy.float32),
               numpy.zeros((nplanes, info['nz'], info['nx'], 3), dtype = numpy.float32),
               numpy.zeros((nplanes, info['nx'], info['ny'], 3), dtype = numpy.float32)]
    planes[0][:, :, :, 0] = info['xnodes'][::(info['nx']/nplanes), numpy.newaxis, numpy.newaxis]
    planes[0][:, :, :, 1] = info['ynodes'][         numpy.newaxis,             :, numpy.newaxis]
    planes[0][:, :, :, 2] = info['znodes'][         numpy.newaxis, numpy.newaxis,             :]
    #planes[1][:, :, :, 0] = info['dx']*numpy.arange(.0, info['nx'], 1)
    #planes[1][:, :, :, 1] = info['dy']*numpy.random.randint(0, info['ny'], size = (nplanes))[:, numpy.newaxis, numpy.newaxis]
    #planes[1][:, :, :, 2] = info['dz']*numpy.arange(.0, info['nz'], 1)
    #planes[2][:, :, :, 0] = info['dx']*numpy.arange(.0, info['nx'], 1)
    #planes[2][:, :, :, 1] = info['dy']*numpy.arange(.0, info['ny'], 1)
    #planes[2][:, :, :, 2] = info['dz']*numpy.random.randint(0, info['nz'], size = (nplanes))[:, numpy.newaxis, numpy.newaxis]
    coordname = ['yz', 'zx', 'xy']
    for i in range(1):
        result = libTDB.getData(.1, planes[i],
            sinterp = 0, tinterp = 0,
            data_set = info['name'], getFunction = 'getVelocityAndPressureSoap')
        fig = plt.figure(figsize = (10.24,10.24))
        ax = fig.add_axes([0, 0, 1, 1], frameon = False)
        ax.set_axis_off()
        ax.imshow(result[0, :, :, 0])
        fig.savefig('plane_' + coordname[i] + '_0.png', format = 'png', dpi = 100)
    return None

def main():
    # load shared library
    lTDB = pyJHTDB.libTDB(libname = 'libTDB',
                        libdir = 'turblib',
                        auth_token = 'edu.jhu.pha.turbulence.testing-201302')
    #initialize webservices
    lTDB.initialize()

    spectra_check(libTDB = lTDB)
    contour_check(libTDB = lTDB,
                  info = pyJHTDB.dbinfo.channel)

    #finalize webservices
    lTDB.finalize()
    return None

#import cProfile
#import pstats

if __name__ == '__main__':
    turbc_clone()
    #main()
    #cProfile.run('main()', 'profile')
    #p = pstats.Stats('profile')
    #p.sort_stats('time').print_stats(10)
    #p.sort_stats('cumulative').print_stats(10)

