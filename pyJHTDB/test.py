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
import ctypes
import math
import numpy

import pyJHTDB
import pyJHTDB.dbinfo
import pyJHTDB.interpolator

if pyJHTDB.found_matplotlib:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
else:
    print('matplotlib is needed for contour plots.'
        + 'You should be able to find installation instructions at http://matplotlib.sourceforge.net')

if pyJHTDB.found_h5py:
    import h5py
else:
    print('h5py is needed for working with cutouts.')

def test_plain(N=10):
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

    print('Coordinates of {0} points where variables are requested:'.format(N))
    for p in range(N):
        print('{0}: {1}'.format(p, points[p]))
    print('Data is requested at time {0}'.format(time))

    # load shared library
    lTDB = pyJHTDB.libJHTDB()
    #initialize webservices
    lTDB.initialize()

    print('Requesting velocity at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getVelocity')
    for p in range(N):
        print('{0}: {1}'.format(p, result[p]))
    print('Requesting forcing at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getForce')
    for p in range(N):
        print('{0}: {1}'.format(p, result[p]))
    print('Requesting velocity and pressure at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = spatialInterp, tinterp = temporalInterp,
            getFunction = 'getVelocityAndPressure')
    for p in range(N):
        print('{0}: v = {1}, p = {2:+}'.format(p, result[p][0:3], result[p][3]))
    print('Requesting velocity gradient at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityGradient')
    for p in range(N):
        print('{0}: '.format(p) +
              'duxdx = {0:+e}, duxdy = {1:+e}, duxdz = {2:+e}\n   '.format(result[p][0], result[p][1], result[p][2]) +
              'duydx = {0:+e}, duydy = {1:+e}, duydz = {2:+e}\n   '.format(result[p][3], result[p][4], result[p][5]) +
              'duzdx = {0:+e}, duzdy = {1:+e}, duzdz = {2:+e}'.format(result[p][6], result[p][7], result[p][8]))
    print('Requesting velocity hessian at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityHessian')
    for p in range(N):
        print('{0}: '.format(p) +
              'd2uxdxdx = {0:+e}, d2uxdxdy = {1:+e}, d2uxdxdz = {2:+e}\n   '.format(result[p][ 0], result[p][ 1], result[p][ 2])
            + 'd2uxdydy = {0:+e}, d2uxdydz = {1:+e}, d2uxdzdz = {2:+e}\n   '.format(result[p][ 3], result[p][ 4], result[p][ 5])
            + 'd2uydxdx = {0:+e}, d2uydxdy = {1:+e}, d2uydxdz = {2:+e}\n   '.format(result[p][ 6], result[p][ 7], result[p][ 8])
            + 'd2uydydy = {0:+e}, d2uydydz = {1:+e}, d2uydzdz = {2:+e}\n   '.format(result[p][ 9], result[p][10], result[p][11])
            + 'd2uzdxdx = {0:+e}, d2uzdxdy = {1:+e}, d2uzdxdz = {2:+e}\n   '.format(result[p][12], result[p][13], result[p][14])
            + 'd2uzdydy = {0:+e}, d2uzdydz = {1:+e}, d2uzdzdz = {2:+e}'.format(result[p][15], result[p][16], result[p][17]))
    print('Requesting velocity laplacian at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityLaplacian')
    for p in range(N):
        print('{0}: '.format(p) +
              'grad2ux = {0:+e}, grad2uy = {1:+e}, grad2uz = {2:+e}, '.format(result[p][0], result[p][1], result[p][2]))
    print('Requesting pressure gradient at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getPressureGradient')
    for p in range(N):
        print('{0}: '.format(p)
            + 'dpdx = {0:+e}, dpdy = {1:+e}, dpdz = {2:+e}, '.format(result[p][0], result[p][1], result[p][2]))
    print('Requesting pressure hessian at {0} points...'.format(N))
    result = lTDB.getData(time, points,
            sinterp = FD4Lag4, tinterp = temporalInterp,
            getFunction = 'getVelocityHessian')
    for p in range(N):
        print('{0}: '.format(p) +
              'd2pdxdx = {0:+e}, d2pdxdy = {1:+e}, d2pdxdz = {2:+e}\n   '.format(result[p][0], result[p][1], result[p][2])
            + 'd2pdydy = {0:+e}, d2pdydz = {1:+e}, d2pdzdz = {2:+e}'.format(result[p][3], result[p][4], result[p][5]))

#    print 'Requesting position at {0} points, starting at time {1} and ending at time {2}...'.format(N, startTime, endTime)
#    result = pyJHTDB.getPosition(startTime, endTime, lag_dt, points, sinterp = spatialInterp)
#    print 'Coordinates of {0} points at startTime:'.format(N)
#    for p in range(N):
#        print p, points[p]
#    print 'Coordinates of {0} points at endTime:'.format(N)
#    for p in range(N):
#        print p, result[p]

    ##  only if matplotlib is present
    if pyJHTDB.found_matplotlib:
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
            getFunction = 'getVelocity')

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

if pyJHTDB.found_matplotlib:

    def spectra_check(
            info = pyJHTDB.dbinfo.isotropic1024coarse,
            lJHTDB = None):
        if not lJHTDB:
            print ('no library given')
            return None
        nlines = 32
        print(nlines, info['nx'], 3)
        print(nlines, info['ny'], 3)
        print(nlines, info['nz'], 3)
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
            result = lJHTDB.getData(.0, lines[i],
                sinterp = 0, tinterp = 0,
                data_set = info['name'], getFunction = 'getVelocityAndPressure')
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
            lJHTDB = None):
        if not lJHTDB:
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
            result = lJHTDB.getData(.1, planes[i],
                sinterp = 0, tinterp = 0,
                data_set = info['name'], getFunction = 'getVelocityAndPressure')
            fig = plt.figure(figsize = (10.24,10.24))
            ax = fig.add_axes([0, 0, 1, 1], frameon = False)
            ax.set_axis_off()
            ax.imshow(result[0, :, :, 0])
            fig.savefig('plane_' + coordname[i] + '_0.png', format = 'png', dpi = 100)
        return None

    def clean_2D_field(
            field_2D,
            dpi = 100,
            figname = 'tst',
            cmap = cm.get_cmap("jet"),
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

    def test_misc():
        # load shared library
        lJHTDB = pyJHTDB.libJHTDB()
        #initialize webservices
        lJHTDB.initialize()
        spectra_check(lJHTDB = lJHTDB)
        contour_check(lJHTDB = lJHTDB,
                      info = pyJHTDB.dbinfo.channel)
        #finalize webservices
        lJHTDB.finalize()
        return None

def test_rawData(
        info = pyJHTDB.dbinfo.channel,
        npoints = 256):

    start = numpy.array([0, 0, 0], dtype = int)
    width = numpy.array([npoints, npoints, 1], dtype = int)

    xg = info['xnodes'][0:width[0]]
    yg = info['ynodes'][0:width[1]]
    zg = info['znodes'][0:width[2]]
    x = numpy.zeros((npoints, npoints, 3), numpy.float32)
    x[:, :, 0] = xg[None, :]
    x[:, :, 1] = yg[:, None]
    x[:, :, 2] = zg[0]

    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    res4 = lJHTDB.getData(
            0,
            x,
            sinterp = 4,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'Velocity')
    res0 = lJHTDB.getRawData(
            0,
            start = start,
            size  = width,
            data_set = info['name'],
            getFunction = 'Velocity')
    lJHTDB.finalize()

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)
    c = ax.contour(res4[:, :, 0])
    ax.clabel(c)
    ax.set_title('Lag4 result from database')
    ax = fig.add_subplot(122)
    c = ax.contour(res0[:, :, 0, 0])
    ax.clabel(c)
    ax.set_title('rawData result from database')
    fig.savefig('tst.pdf', format = 'pdf')
    return None

def test_interp_1D(
        info = pyJHTDB.dbinfo.channel,
        m = 1,
        q = 4,
        npoints = 256):

    start = numpy.array([0, 0, 0], dtype = int)
    width = numpy.array([51+q, 37+q, 17+q], dtype = int)

    i = pyJHTDB.interpolator.spline_interpolator(
            info = info,
            n = (q - 2)/2,
            m = m)
    i.generate_clib()

    xg = numpy.linspace(info['xnodes'][i.nx+1],
                        info['xnodes'][width[0] - i.nx - 1],
                        npoints)
    if info['yperiodic']:
        yg = numpy.linspace(info['ynodes'][i.ny+1],
                            info['ynodes'][width[1] - i.ny - 1],
                            npoints)
    else:
        yg = numpy.linspace(info['ynodes'][0],
                            info['ynodes'][width[1] - i.ny - 1],
                            npoints)
    zg = numpy.linspace(info['znodes'][i.nz+1],
                        info['znodes'][width[2] - i.nz - 1],
                        npoints)
    x = numpy.zeros((npoints, 3), numpy.float32)
    x[:, 0] = xg[0]
    x[:, 1] = yg[:]
    x[:, 2] = zg[0]

    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    # get raw data to interpolate
    test_field = lJHTDB.getRawData(
            0,
            start = start,
            size  = width,
            data_set = info['name'],
            getFunction = 'Velocity')
    # get Lag8 velocity
    res0 = lJHTDB.getData(
            0,
            x,
            sinterp = 8,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocity')
    # get locally interpolated values
    res1 = i.cinterpolate(
            x, test_field)
    # get Lag4 gradient
    resd0 = lJHTDB.getData(
            0,
            x,
            sinterp = 44,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocityGradient')
    # get locally interpolated gradient
    resdx1 = i.cinterpolate(
            x,
            test_field,
            diff = [1, 0, 0])
    resdy1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 1, 0])
    resdz1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 0, 1])
    resd1 = resd0.copy()
    resd1[..., 0] = resdx1[..., 0]
    resd1[..., 1] = resdy1[..., 0]
    resd1[..., 2] = resdz1[..., 0]
    resd1[..., 3] = resdx1[..., 1]
    resd1[..., 4] = resdy1[..., 1]
    resd1[..., 5] = resdz1[..., 1]
    resd1[..., 6] = resdx1[..., 2]
    resd1[..., 7] = resdy1[..., 2]
    resd1[..., 8] = resdz1[..., 2]
    del resdx1, resdy1, resdz1
    lJHTDB.finalize()

    if pyJHTDB.found_matplotlib:
        def compare_results(
                fld0,
                fld1,
                figname = 'tst'):
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
            ax.plot(fld0, color = 'blue')
            ax.plot(fld1, color = 'red')
            fig.savefig(figname + '.pdf', format = 'pdf')

        compare_results(
                res0,
                res1,
                figname = 'tst')

        compare_results(
                resd0,
                resd1,
                figname = 'dtst')
    dist = (numpy.average(numpy.sqrt(numpy.sum((res0 - res1)**2, axis = 1))) /
            numpy.average(numpy.sqrt(numpy.sum((res0)**2, axis = 1))))
    print ('average distance for result {0}'.format(dist))
    ddist = (numpy.average(numpy.sqrt(numpy.sum((resd0 - resd1)**2, axis = 1))) /
             numpy.average(numpy.sqrt(numpy.sum((resd0)**2, axis = 1))))
    print ('average distance for dresult {0}'.format(ddist))
    return res0, res1, resd0, resd1

def test_interp_2D(
        info = pyJHTDB.dbinfo.channel,
        m = 1,
        q = 4,
        npoints = 256):

    start = numpy.array([0, 0, 0], dtype = int)
    width = numpy.array([91, 67, 11], dtype = int)

    i = pyJHTDB.interpolator.spline_interpolator(
            info = info,
            n = (q - 2)/2,
            m = m)
    i.generate_clib()

    xg = numpy.linspace(info['xnodes'][i.n+1], info['xnodes'][width[0] - i.n - 1], npoints)
    if info['yperiodic']:
        yg = numpy.linspace(info['ynodes'][i.n+1], info['ynodes'][width[1] - i.n - 1], npoints)
    else:
        yg = numpy.linspace(info['ynodes'][0], info['ynodes'][width[1] - i.n - 1], npoints)
    zg = numpy.linspace(info['znodes'][i.n+1], info['znodes'][width[2] - i.n - 1], npoints)
    x = numpy.zeros((npoints, npoints, 3), numpy.float32)
    x[:, :, 0] = xg[0] #None, :]
    x[:, :, 1] = yg[:, None]
    x[:, :, 2] = zg[0] #zg[:, None]

    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    # get raw data to interpolate
    test_field = lJHTDB.getRawData(
            0,
            start = start,
            size  = width,
            data_set = info['name'],
            getFunction = 'Velocity')
    # get Lag8 velocity
    res0 = lJHTDB.getData(
            0,
            x,
            sinterp = 8,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocity')
    # get locally interpolated values
    res1 = i.cinterpolate(
            x, test_field)
    # get Lag4 gradient
    resd0 = lJHTDB.getData(
            0,
            x,
            sinterp = 44,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocityGradient')
    # get locally interpolated gradient
    resdx1 = i.cinterpolate(
            x,
            test_field,
            diff = [1, 0, 0])
    resdy1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 1, 0])
    resdz1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 0, 1])
    resd1 = resd0.copy()
    resd1[..., 0] = resdx1[..., 0]
    resd1[..., 1] = resdy1[..., 0]
    resd1[..., 2] = resdz1[..., 0]
    resd1[..., 3] = resdx1[..., 1]
    resd1[..., 4] = resdy1[..., 1]
    resd1[..., 5] = resdz1[..., 1]
    resd1[..., 6] = resdx1[..., 2]
    resd1[..., 7] = resdy1[..., 2]
    resd1[..., 8] = resdz1[..., 2]
    del resdx1, resdy1, resdz1
    lJHTDB.finalize()

    def compare_results(
            fld0,
            fld1,
            figname = 'tst'):
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(121)
        c = ax.contour(
                #x[:, :, 0],
                #x[:, :, 1],
                fld0)
        ax.clabel(c)
        ax.set_title('result from database')
        ax = fig.add_subplot(122)
        c = ax.contour(
                #x[:, :, 0],
                #x[:, :, 1],
                fld1, c.levels)
        ax.clabel(c)
        ax.set_title('local M{0}Q{1:0>2} result'.format(i.m, i.n*2+2))
        fig.savefig(figname + '.pdf', format = 'pdf')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(fld0[:, 0], color = 'blue')
        ax.plot(fld1[:, 0], color = 'red')
        fig.savefig(figname + '_1D.pdf', format = 'pdf')

    compare_results(
            res0[:, :, 1],
            res1[:, :, 1],
            figname = 'tst')

    compare_results(
            resd0[:, :, 1],
            resd1[:, :, 1],
            figname = 'dtst')
    ddist = (numpy.average(numpy.sqrt(numpy.sum((resd0 - resd1)**2, axis = 1))) /
             numpy.average(numpy.sqrt(numpy.sum((resd0)**2, axis = 1))))
    print (ddist)

    return res0, res1, resd0, resd1

def test_divfree(
        info = pyJHTDB.dbinfo.channel,
        m = 1,
        q = 4,
        npoints = 256,
        dbinterp = 44):

    start = numpy.array([0, 0, 0], dtype = int)
    width = numpy.array([91, 67, 31], dtype = int)

    i = pyJHTDB.interpolator.spline_interpolator(
            info = info,
            n = (q - 2)/2,
            m = m)
    i.generate_clib()

    xg = [info['xnodes'][i.n+1], info['xnodes'][width[0] - i.n - 1]]
    if info['yperiodic']:
        yg = [info['ynodes'][i.n+1], info['ynodes'][width[1] - i.n - 1]]
    else:
        yg = [info['ynodes'][0], info['ynodes'][width[1] - i.n - 1]]
    zg = [info['znodes'][i.n+1], info['znodes'][width[2] - i.n - 1]]
    x = numpy.random.random(size = (npoints, 3)).astype(numpy.float32)
    x[:, 0] = xg[0] + x[:, 0]*(xg[1] - xg[0])
    x[:, 1] = yg[0] + x[:, 1]*(yg[1] - yg[0])
    x[:, 2] = zg[0] + x[:, 2]*(zg[1] - zg[0])

    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    # get raw data to interpolate
    test_field = lJHTDB.getRawData(
            0,
            start = start,
            size  = width,
            data_set = info['name'],
            getFunction = 'Velocity')
    # get Lag4 gradient
    resd0 = lJHTDB.getData(
            0,
            x,
            sinterp = dbinterp,
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocityGradient')
    # get locally interpolated gradient
    resdx1 = i.cinterpolate(
            x,
            test_field,
            diff = [1, 0, 0])
    resdy1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 1, 0])
    resdz1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 0, 1])
    resd1 = resd0.copy()
    resd1[..., 0] = resdx1[..., 0]
    resd1[..., 1] = resdy1[..., 0]
    resd1[..., 2] = resdz1[..., 0]
    resd1[..., 3] = resdx1[..., 1]
    resd1[..., 4] = resdy1[..., 1]
    resd1[..., 5] = resdz1[..., 1]
    resd1[..., 6] = resdx1[..., 2]
    resd1[..., 7] = resdy1[..., 2]
    resd1[..., 8] = resdz1[..., 2]
    del resdx1, resdy1, resdz1
    lJHTDB.finalize()

    dmagnitude = numpy.average(numpy.sqrt(numpy.sum((resd0)**2, axis = 1)))
    ddist = numpy.average(numpy.sqrt(numpy.sum((resd0 - resd1)**2, axis = 1))) / dmagnitude
    print ('average relative distance between dbinterp and M{0}Q{1} is {2}'.format(m, q, ddist))

    div0 = resd0[:, 0] + resd0[:, 4] + resd0[:, 8]
    div1 = resd1[:, 0] + resd1[:, 4] + resd1[:, 8]
    print('average divergence for dbinterp is {0}'.format(numpy.average(div0**2) / dmagnitude))
    print('average divergence for M{0}Q{1} is {2}'.format(m, q, numpy.average(div1**2) / dmagnitude))
    return None

def test_local_vs_db_interp(
        info = pyJHTDB.dbinfo.channel,
        time = 0.0,
        m = 1,
        q = 4,
        npoints = 256,
        dbinterp = [8, 44],
        start = numpy.array([0, 0, 0], dtype = int),
        width = numpy.array([91, 67, 31], dtype = int),
        messages_on = False,
        token = "edu.jhu.pha.turbulence.testing-201311"):

    i = pyJHTDB.interpolator.spline_interpolator(
            info = info,
            n = (q - 2)/2,
            m = m)
    i.generate_clib()

    # build point array
    xg = [info['xnodes'][start[0] + i.n+1], info['xnodes'][start[0] + width[0] - i.n - 1]]
    if info['yperiodic']:
        yg = [info['ynodes'][start[1] + i.n+1], info['ynodes'][start[1] + width[1] - i.n - 1]]
    else:
        yg = [info['ynodes'][start[1]        ], info['ynodes'][start[1] + width[1] - i.n - 1]]
    zg = [info['znodes'][start[2] + i.n+1], info['znodes'][start[2] + width[2] - i.n - 1]]
    x = numpy.random.random(size = (npoints, 3)).astype(numpy.float32)
    x[:, 0] = xg[0] + x[:, 0]*(xg[1] - xg[0])
    x[:, 1] = yg[0] + x[:, 1]*(yg[1] - yg[0])
    x[:, 2] = zg[0] + x[:, 2]*(zg[1] - zg[0])

    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    lJHTDB.add_token(token)

    # get raw data to interpolate
    test_field = lJHTDB.getRawData(
            int(time),
            start = start,
            size  = width,
            data_set = info['name'],
            getFunction = 'Velocity')
    # get DB field
    res0 = lJHTDB.getData(
            time,
            x,
            sinterp = dbinterp[0],
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocity')
    # get DB gradient
    resd0 = lJHTDB.getData(
            time,
            x,
            sinterp = dbinterp[1],
            tinterp = 0,
            data_set = info['name'],
            getFunction = 'getVelocityGradient')
    # get locally interpolated field
    res1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 0, 0],
            field_offset = start)
    # get locally interpolated gradient
    resdx1 = i.cinterpolate(
            x,
            test_field,
            diff = [1, 0, 0],
            field_offset = start)
    resdy1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 1, 0],
            field_offset = start)
    resdz1 = i.cinterpolate(
            x,
            test_field,
            diff = [0, 0, 1],
            field_offset = start)
    resd1 = resd0.copy()
    resd1[..., 0] = resdx1[..., 0]
    resd1[..., 1] = resdy1[..., 0]
    resd1[..., 2] = resdz1[..., 0]
    resd1[..., 3] = resdx1[..., 1]
    resd1[..., 4] = resdy1[..., 1]
    resd1[..., 5] = resdz1[..., 1]
    resd1[..., 6] = resdx1[..., 2]
    resd1[..., 7] = resdy1[..., 2]
    resd1[..., 8] = resdz1[..., 2]
    del resdx1, resdy1, resdz1
    lJHTDB.finalize()

    if messages_on:
        comp0 = ['ux', 'uy', 'uz']
        comp1 = ['dxux', 'dyux', 'dzux',
                 'dxuy', 'dyuy', 'dzuy',
                 'dxuz', 'dyuz', 'dzuz']

        print ('printing average relative distance between DB and local')
        print ('example point is {0}'.format(x[0]))
        print ('for direct interpolation using (DB) {0} and (local) M{1}Q{2}'.format(dbinterp[0], m, q))
        print ('printing average((DB) - (local)) / average(DB), (DB) at example point, abs((DB) - (local)) at example point ')
        for i in range(3):
            magnitude = numpy.average(numpy.abs(res0[:, i]))
            distance  = numpy.average(numpy.abs(res0[:, i] - res1[:, i])) / magnitude
            print (comp0[i] + ' ' +
                   '{0}, {1:+}, {2}'.format(distance, res0[0, i], numpy.abs(res0[0, i] - res1[0, i])))
        print ('for gradient interpolation using (DB) {0} and (local) M{1}Q{2}'.format(dbinterp[1], m, q))
        for i in range(9):
            magnitude = numpy.average(numpy.abs(resd0[:, i]))
            distance  = numpy.average(numpy.abs(resd0[:, i] - resd1[:, i])) / magnitude
            print (comp1[i] + ' ' +
                   '{0}, {1:+}, {2}'.format(distance, resd0[0, i], numpy.abs(resd0[0, i] - resd1[0, i])))
    return res0, res1, resd0, resd1

class LocalInterpTest:
    def __init__(
            self,
            info):
        self.info = info
        return None
    def set_up_field(
            self,
            xnodes = [0,  64],
            ynodes = [0, 128],
            znodes = [0,  64],
            buffer_size = 16):
        self.buffer_size = 16
        full_frame = h5py.File(
            '/stuff/data/{0}/{0}_t0000.h5'.format(self.info['name']),
            'r')
        self.xnodes = xnodes
        self.ynodes = ynodes
        self.znodes = znodes
        self.test_field = full_frame['u00000'][
            self.znodes[0]:self.znodes[1],
            self.ynodes[0]:self.ynodes[1],
            self.xnodes[0]:self.xnodes[1]].copy()
        self.y0buffer = min(self.buffer_size, self.ynodes[0])
        full_frame.close()
        return None
    def set_up_points(
            self,
            npoints = 2**5,
            yval = None):
        self.npoints = npoints
        if type(yval) == type(None):
            if self.info['yperiodic']:
                self.yindices = range(self.ynodes[0] + self.buffer_size,
                                      self.ynodes[1] - self.buffer_size)
            else:
                self.yindices = range(self.ynodes[0] + self.y0buffer,
                                      self.ynodes[1] - self.buffer_size)
            self.p = numpy.random.random(
                size = (npoints, self.yindices.count, 3)).astype(numpy.float32)
            if self.info['yperiodic']:
                self.p[..., 1] *= self.info['dy']
            else:
                self.p[..., 1] *= self.info['dy'][None, self.yindices]
            self.p[..., 1] += self.info['ynodes'][None, self.yindices]
        else:
            self.p = numpy.random.random(
                size = (npoints, 3)).astype(numpy.float32)
            if yval == 'random':
                if self.info['yperiodic']:
                    self.p[..., 1] = (
                        self.info['ynodes'][self.ynodes[0] + self.buffer_size] +
                        self.p[..., 1]*(
                            self.info['ynodes'][self.ynodes[1] - self.buffer_size] -
                            self.info['ynodes'][self.ynodes[0] + self.buffer_size]))
                else:
                    self.p[..., 1] = (
                        self.info['ynodes'][self.ynodes[0] + self.y0buffer] +
                        self.p[..., 1]*(
                            self.info['ynodes'][self.ynodes[1] - self.buffer_size] -
                            self.info['ynodes'][self.ynodes[0] + self.y0buffer]))
            else:
                self.p[..., 1] = yval
        self.p[..., 0] = (
            self.info['xnodes'][self.xnodes[0] + self.buffer_size] +
            self.p[..., 0]*(
                self.info['xnodes'][self.xnodes[1] - self.buffer_size] -
                self.info['xnodes'][self.xnodes[0] + self.buffer_size]))
        self.p[..., 2] = (
            self.info['znodes'][self.znodes[0] + self.buffer_size] +
            self.p[..., 2]*(
                self.info['znodes'][self.znodes[1] - self.buffer_size] -
                self.info['znodes'][self.znodes[0] + self.buffer_size]))
        return None
    def set_up_interpolators(
            self,
            pars = [[12, 3, 3, 2],
                    [12, 4, 4, 2],
                    [12, 5, 5, 2],
                    [12, 6, 6, 2],
                    [12, 7, 7, 2],
                    [12, 8, 8, 2],
                    [12, 9, 9, 2],
                    [12,10,10, 2]]):
        self.interp = []
        self.keys   = []
        for par in pars:
            self.keys.append(
                'nx{0:0>2}_ny{1:0>2}_nz{2:0>2}_m{3}'.format(
                    par[0], par[1], par[2], par[3]))
            self.interp.append(
                pyJHTDB.interpolator.spline_interpolator(
                    info = self.info,
                    nx = par[0],
                    ny = par[1],
                    nz = par[2],
                    m = par[3],
                    initialize = False,
                    cformula_unroll = False))
        return None
    def interpolate(
            self):
        self.uval = []
        self.gradu = []
        for k in range(len(self.keys)):
            self.uval.append(self.interp[k].cinterpolate(
                self.p,
                self.test_field,
                diff = [0, 0, 0],
                field_offset = [self.xnodes[0],
                                self.ynodes[0],
                                self.znodes[0]]))
            dxvel = self.interp[k].cinterpolate(
                self.p,
                self.test_field,
                diff = [1, 0, 0],
                field_offset = [self.xnodes[0],
                                self.ynodes[0],
                                self.znodes[0]])
            dyvel = self.interp[k].cinterpolate(
                self.p,
                self.test_field,
                diff = [0, 1, 0],
                field_offset = [self.xnodes[0],
                                self.ynodes[0],
                                self.znodes[0]])
            dzvel = self.interp[k].cinterpolate(
                self.p,
                self.test_field,
                diff = [0, 0, 1],
                field_offset = [self.xnodes[0],
                                self.ynodes[0],
                                self.znodes[0]])
            self.gradu.append(numpy.zeros(
                dxvel.shape[:-1] + (9,),
                dtype = dxvel.dtype))
            self.gradu[k][..., 0] = dxvel[..., 0]
            self.gradu[k][..., 1] = dyvel[..., 0]
            self.gradu[k][..., 2] = dzvel[..., 0]
            self.gradu[k][..., 3] = dxvel[..., 1]
            self.gradu[k][..., 4] = dyvel[..., 1]
            self.gradu[k][..., 5] = dzvel[..., 1]
            self.gradu[k][..., 6] = dxvel[..., 2]
            self.gradu[k][..., 7] = dyvel[..., 2]
            self.gradu[k][..., 8] = dzvel[..., 2]
        self.uval  = numpy.array(self.uval)
        self.gradu = numpy.array(self.gradu)
        return None
    def get_divergence(
            self):
        self.divu = []
        self.divu_upper = []
        self.divu_lower = []
        for k in range(len(self.keys)):
            factor = (numpy.sqrt(3) /
                numpy.average(numpy.sqrt(
                    self.gradu[k][..., 0]**2 +
                    self.gradu[k][..., 4]**2 +
                    self.gradu[k][..., 8]**2), axis = 0))
            self.divu.append(factor*numpy.average(numpy.abs(
                    self.gradu[k][..., 0] +
                    self.gradu[k][..., 4] +
                    self.gradu[k][..., 8]), axis = 0))
            self.divu_upper.append(factor*numpy.percentile(numpy.abs(
                    self.gradu[k][..., 0] +
                    self.gradu[k][..., 4] +
                    self.gradu[k][..., 8]), 90, axis = 0))
            self.divu_lower.append(factor*numpy.percentile(numpy.abs(
                    self.gradu[k][..., 0] +
                    self.gradu[k][..., 4] +
                    self.gradu[k][..., 8]), 10, axis = 0))
        self.divu = numpy.array(self.divu)
        self.divu_upper = numpy.array(self.divu_upper)
        self.divu_lower = numpy.array(self.divu_lower)
        return None

if __name__ == '__main__':
    test_plain()

