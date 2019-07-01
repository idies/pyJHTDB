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

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import pyJHTDB
import pyJHTDB.dbinfo
import pyJHTDB.interpolator
#import matplotlib_extra as me

def generate_interpolation_points(
        info,
        npoints = 2**5,
        randseeds = [1],
        dir_name = 'test',
        xfixed = None,
        yfixed = None,
        zfixed = None,
        parallelogram = None):
    dir_name = info['name'] + '_' + dir_name + '/'
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    if not parallelogram:
        points = np.zeros((npoints, 3), dtype = np.float32)
        for randseed in randseeds:
            np.random.seed(randseed)
            if xfixed:
                points[:, 0] = xfixed
            else:
                points[:, 0] = np.random.random(npoints)*info['lx'] + info['xnodes'][0]
            if yfixed:
                points[:, 1] = yfixed
            else:
                points[:, 1] = np.random.random(npoints)*info['ly'] + info['ynodes'][0]
            if zfixed:
                points[:, 2] = zfixed
            else:
                points[:, 2] = np.random.random(npoints)*info['lz'] + info['znodes'][0]
            np.save(dir_name
                    + 'points_{0:0>3x}_{1:0>3x}'.format(npoints,
                                                        randseed),
                       points)
    else:
        points = np.zeros((npoints, npoints, 3), dtype = np.float32)
        indices = np.arange(.0, npoints, 1).astype(np.float32)
        points[:, :] = (parallelogram['origin'][np.newaxis, np.newaxis, :]
                      + parallelogram['a'][np.newaxis, np.newaxis, :] * indices[np.newaxis, :, np.newaxis]
                      + parallelogram['b'][np.newaxis, np.newaxis, :] * indices[:, np.newaxis, np.newaxis])
        np.save(dir_name
                + 'points_' + parallelogram['name'],
                   points)
    return points

def interpolate_spline(
        lTDB, info,
        npoints = 2**5,
        spline_neighbours = 7,
        spline_continuity = 2,
        randseeds = [1],
        dir_name = 'test'):
    print(('just entered interpolate_spline, dataset is {0}, '
         + 'and (n, m) are ({1}, {2})').format(
             info['name'],
             spline_neighbours,
             spline_continuity))
    dir_name = info['name'] + '_' + dir_name + '/'
    s = None
    result_list = []
    for randseed in randseeds:
        points = np.load(dir_name
                          + 'points_{0:0>3x}_{1:0>3x}.npy'.format(npoints,
                                                                  randseed))
        file_name = (dir_name
                   + 'res_n{0}m{1}_{2:0>3x}_{3:0>3x}'.format(spline_neighbours,
                                                             spline_continuity,
                                                             npoints,
                                                             randseed))
        if not os.path.exists(file_name + '.npy'):
            if s == None:
                s = pyJHTDB.interpolator.spline_interpolator(
                        info,
                        n = spline_neighbours,
                        m = spline_continuity)
            ress = s(lTDB = lTDB,
                     time = .1,
                     points = points,
                     dorder = [(0, 0, 0),  # regular
                               (1, 0, 0),  # gradient 0
                               (0, 1, 0),  # gradient 1
                               (0, 0, 1),  # gradient 2
                               (2, 0, 0),  # hessian  0
                               (1, 1, 0),  # hessian  1
                               (1, 0, 1),  # hessian  2
                               (0, 2, 0),  # hessian  3
                               (0, 1, 1),  # hessian  4
                               (0, 0, 2)], # hessian  5
                     getFunction = 'getVelocityAndPressure')
            np.save(file_name, ress)
        else:
            ress = np.load(file_name + '.npy')
        result_list.append(ress)
    return np.concatenate(result_list, axis = 1)

def spline_error(
        lTDB, info,
        npoints = 2**4,
        nmax = 5,
        mmax = 3,
        randseeds = [1],
        dir_name = 'test',
        plots_on = False):
    results = []
    neighbours = []
    for m in range(2, mmax+1):
        neighbours.append(np.array(range((m+1)/2, nmax)))
        results.append(np.array([interpolate_spline(
                                         lTDB, info,
                                         npoints = npoints,
                                         spline_neighbours = n,
                                         spline_continuity = m,
                                         randseeds = randseeds,
                                         dir_name = dir_name) for n in range((m+1)/2, nmax+1)]))
    if plots_on:
        diff = [(0, 0, 0),
                (1, 0, 0),
                (0, 1, 0),
                (0, 0, 1),
                (2, 0, 0),
                (1, 1, 0),
                (1, 0, 1),
                (0, 2, 0),
                (0, 1, 1),
                (0, 0, 2)]
        def get_err(res):
            return (np.average(np.abs(
                        res[1:] - res[:res.shape[0]-1]), axis = 2)
                  / np.average(np.abs(
                        res[1]), axis = 1))
        fig1 = plt.figure(figsize=(12,6))
        ax1 = fig1.add_axes([.05, .1, .9, .8])
        fig2 = plt.figure(figsize=(12,6))
        ax2 = fig2.add_axes([.05, .1, .9, .8])
        ax1.set_title('relative interpolation error for $u$ averaged over {0} points'.format(npoints*len(randseeds)))
        ax2.set_title('relative interpolation error for $p$ averaged over {0} points'.format(npoints*len(randseeds)))
        marker = ['.', '*']
        line_style = [[(None, None), (0.0, 0.0, 0.0)],
                      [(None, None), (1.0, 0.0, 0.0)],
                      [(None, None), (0.0, 1.0, 0.0)],
                      [(None, None), (0.0, 0.0, 1.0)],
                      [       (1,1), (1.0, 0.0, 0.0)],
                      [       (2,2), (0.5, 0.5, 0.0)],
                      [       (2,2), (0.5, 0.0, 0.5)],
                      [       (1,1), (0.0, 1.0, 0.0)],
                      [       (2,2), (0.0, 0.5, 0.5)],
                      [       (1,1), (0.0, 0.0, 1.0)]]
        for i in range(len(results)):
            err = get_err(results[i])
            for deriv in range(len(diff)):
                ax1.plot(neighbours[i],
                    np.log2(np.average(err[:, deriv, :3], axis = 1)),
                    marker = marker[m%2],
                    dashes = line_style[deriv][0],
                    color  = line_style[deriv][1],
                    label = '$m = {0}$, deriv = ${1}$'.format(m, diff[deriv]))
                ax2.plot(neighbours[i],
                    np.log2(err[:, deriv, 3]),
                    marker = marker[m%2],
                    dashes = line_style[deriv][0],
                    color  = line_style[deriv][1],
                    label = '$m = {0}$, deriv = ${1}$'.format(m, diff[deriv]))
        ax1.set_xscale('log')
        ax1.set_xlim(.95, 25)
        ax1.set_ylim(-14, 0)
        ax1.set_xlabel('neighbours in formula')
        ax1.set_ylabel('$\\log_2 \\langle \\varepsilon \\rangle$')
        ax1.legend(loc = 'best', ncol = mmax-1)
        ax2.set_xscale('log')
        ax2.set_xlim(.95, 25)
        ax2.set_ylim(-14, 0)
        ax2.set_xlabel('neighbours in formula')
        ax2.legend(loc = 'best', ncol = mmax-1)
        dir_name = info['name'] + '_' + dir_name + '/'
        fig1.savefig(dir_name + 'interp_err_u.pdf', format='pdf')
        fig2.savefig(dir_name + 'interp_err_p.pdf', format='pdf')
        fig1.clear()
        fig2.clear()
    return None

def Lagrange_err_histograms(
        lTDB, info,
        npoints = 2**5,
        compute = False,
        spline_neighbours = 7,
        spline_continuity = 2,
        randseeds = [1],
        dir_name = 'test',
        plots_on = False):
    ress = interpolate_spline(
                lTDB, info,
                npoints = npoints,
                spline_neighbours = spline_neighbours,
                spline_continuity = spline_continuity,
                randseeds = randseeds,
                dir_name = dir_name)
    ress_bad = interpolate_spline(
                lTDB, info,
                npoints = npoints,
                spline_neighbours = spline_neighbours - 1,
                spline_continuity = spline_continuity,
                randseeds = randseeds,
                dir_name = dir_name)
    dir_name = info['name'] + '_' + dir_name + '/'
    res4 = []
    res6 = []
    res8 = []
    grad = []
    hess = []
    for randseed in randseeds:
        fnames = [dir_name + 'res_Lag4_{0:0>3x}_{1:0>3x}'.format(npoints, randseed),
                  dir_name + 'res_Lag6_{0:0>3x}_{1:0>3x}'.format(npoints, randseed),
                  dir_name + 'res_Lag8_{0:0>3x}_{1:0>3x}'.format(npoints, randseed),
                  dir_name + 'grad_Lag_{0:0>3x}_{1:0>3x}'.format(npoints, randseed),
                  dir_name + 'hess_Lag_{0:0>3x}_{1:0>3x}'.format(npoints, randseed)]
        points = np.load(dir_name
                          + 'points_{0:0>3x}_{1:0>3x}.npy'.format(npoints,
                                                                  randseed))
        if not os.path.exists(fnames[0] + '.npy'):
            res4.append(lTDB.getData(.1,
                    points,
                    sinterp = 4, tinterp = 0,
                    data_set = info['name'],
                    getFunction = 'getVelocityAndPressure'))
            np.save(fnames[0], res4[-1])
        else:
            res4.append(np.load(fnames[0] + '.npy'))
        if not os.path.exists(fnames[1] + '.npy'):
            res6.append(lTDB.getData(.1,
                    points,
                    sinterp = 6, tinterp = 0,
                    data_set = info['name'],
                    getFunction = 'getVelocityAndPressure'))
            np.save(fnames[1], res6[-1])
        else:
            res6.append(np.load(fnames[1] + '.npy'))
        if not os.path.exists(fnames[2] + '.npy'):
            res8.append(lTDB.getData(.1,
                    points,
                    sinterp = 8, tinterp = 0,
                    data_set = info['name'],
                    getFunction = 'getVelocityAndPressure'))
            np.save(fnames[2], res8[-1])
        else:
            res8.append(np.load(fnames[2] + '.npy'))
        if not os.path.exists(fnames[3] + '.npy'):
            grad.append(lTDB.getData(.1,
                    points,
                    sinterp = 44, tinterp = 0,
                    data_set = info['name'],
                    getFunction = 'getVelocityGradient'))
            np.save(fnames[3], grad[-1])
        else:
            grad.append(np.load(fnames[3] + '.npy'))
        if not os.path.exists(fnames[4] + '.npy'):
            hess.append(lTDB.getData(.1,
                    points,
                    sinterp = 44, tinterp = 0,
                    data_set = info['name'],
                    getFunction = 'getPressureHessian'))
            np.save(fnames[4], hess[-1])
        else:
            hess.append(np.load(fnames[4] + '.npy'))
    res4 = np.concatenate(res4, axis = 0)
    res6 = np.concatenate(res6, axis = 0)
    res8 = np.concatenate(res8, axis = 0)
    grad = np.concatenate(grad, axis = 0)
    hess = np.concatenate(hess, axis = 0)
    if plots_on:
        err4 = np.abs(ress[0] - res4)
        err6 = np.abs(ress[0] - res6)
        err8 = np.abs(ress[0] - res8)
        errs = np.abs(ress[0] - ress_bad[0])
        def get_denominator(array):
            tval = np.average(np.abs(array))
            return np.where(np.abs(array) > tval, np.abs(array), tval)
        for coord in range(4):
            tarr = get_denominator(ress[0, :, coord])
            err4[:, coord] /= tarr
            err6[:, coord] /= tarr
            err8[:, coord] /= tarr
            errs[:, coord] /= tarr
        #fig = plt.figure(figsize=(3,5))
        #ax = fig.add_axes([.1, .1, .8, .8])
        #ax.set_title('rel interp err averaged over {0} points'.format(res4.shape[0]))
        #neighbours = np.array([1, 2, 3])
        #err = np.array([err4, err6, err8])
        #ax.plot(neighbours,
        #    np.average(np.average(err[:, :, :3], axis = 1), axis = 1),
        #    label = 'Lagrange $u$',
        #    marker = '.')
        #ax.plot(neighbours,
        #    np.average(err[:, :, 3], axis = 1),
        #    label = 'Lagrange $p$',
        #    marker = '.')
        #ax.set_yscale('log')
        #ax.set_xscale('log')
        #ax.set_xlabel('neighbours in formula')
        #ax.legend(loc = 'best')
        #fig.savefig(dir_name + 'err_Lag.pdf', format='pdf')
        #fig.clear()
        errg = grad.copy()
        errg[:, 0] = np.abs(ress[1, :, 0] - grad[:, 0]) / get_denominator(ress[1, :, 0])
        errg[:, 1] = np.abs(ress[2, :, 0] - grad[:, 1]) / get_denominator(ress[2, :, 0])
        errg[:, 2] = np.abs(ress[3, :, 0] - grad[:, 2]) / get_denominator(ress[3, :, 0])
        errg[:, 3] = np.abs(ress[1, :, 1] - grad[:, 3]) / get_denominator(ress[1, :, 1])
        errg[:, 4] = np.abs(ress[2, :, 1] - grad[:, 4]) / get_denominator(ress[2, :, 1])
        errg[:, 5] = np.abs(ress[3, :, 1] - grad[:, 5]) / get_denominator(ress[3, :, 1])
        errg[:, 6] = np.abs(ress[1, :, 2] - grad[:, 6]) / get_denominator(ress[1, :, 2])
        errg[:, 7] = np.abs(ress[2, :, 2] - grad[:, 7]) / get_denominator(ress[2, :, 2])
        errg[:, 8] = np.abs(ress[3, :, 2] - grad[:, 8]) / get_denominator(ress[3, :, 2])
        errh = hess.copy()
        errh[:, 0] = np.abs(ress[4, :, 3] - hess[:, 0]) / get_denominator(ress[4, :, 3])
        errh[:, 1] = np.abs(ress[5, :, 3] - hess[:, 1]) / get_denominator(ress[5, :, 3])
        errh[:, 2] = np.abs(ress[6, :, 3] - hess[:, 2]) / get_denominator(ress[6, :, 3])
        errh[:, 3] = np.abs(ress[7, :, 3] - hess[:, 3]) / get_denominator(ress[7, :, 3])
        errh[:, 4] = np.abs(ress[8, :, 3] - hess[:, 4]) / get_denominator(ress[8, :, 3])
        errh[:, 5] = np.abs(ress[9, :, 3] - hess[:, 5]) / get_denominator(ress[9, :, 3])
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_axes([.05, .1, .9, .8])
        # epsilon added because otherwise there will be NaNs...
        epsilon = 2.**(-20)
        err_data = [[err4[:, :3].flatten() + epsilon, 'Lag4 $u$'                   , (0.2, 0.2, 0.2), 'dotted' , np.log2(np.average(errs[:, :3].flatten()))],
                    [err6[:, :3].flatten() + epsilon, 'Lag6 $u$'                   , (0.5, 0.5, 0.5), 'dotted' , np.log2(np.average(errs[:, :3].flatten()))],
                    [err8[:, :3].flatten() + epsilon, 'Lag8 $u$'                   , (0.8, 0.8, 0.8), 'dotted' , np.log2(np.average(errs[:, :3].flatten()))],
                    [err4[:,  3]           + epsilon, 'Lag4 $p$'                   , (0.2, 0.2, 0.2), 'dashdot', np.log2(np.average(errs[:, 3]))],
                    [err6[:,  3]           + epsilon, 'Lag6 $p$'                   , (0.5, 0.5, 0.5), 'dashdot', np.log2(np.average(errs[:, 3]))],
                    [err8[:,  3]           + epsilon, 'Lag8 $p$'                   , (0.8, 0.8, 0.8), 'dashdot', np.log2(np.average(errs[:, 3]))],
                    [errg[:,  0]           + epsilon, '$\\partial_x u_x$'          , (1.0, 0.0, 0.0), 'solid'  , np.log2(np.average(np.abs(ress[1, :, 0] - ress_bad[1, :, 0]) / get_denominator(ress[1, :, 0])))],
                    [errg[:,  1]           + epsilon, '$\\partial_y u_x$'          , (0.7, 0.5, 0.0), 'solid'  , np.log2(np.average(np.abs(ress[2, :, 0] - ress_bad[2, :, 0]) / get_denominator(ress[2, :, 0])))],
                    [errg[:,  2]           + epsilon, '$\\partial_z u_x$'          , (0.7, 0.0, 0.5), 'solid'  , np.log2(np.average(np.abs(ress[3, :, 0] - ress_bad[3, :, 0]) / get_denominator(ress[3, :, 0])))],
                    [errg[:,  3]           + epsilon, '$\\partial_x u_y$'          , (0.5, 0.7, 0.0), 'solid'  , np.log2(np.average(np.abs(ress[1, :, 1] - ress_bad[1, :, 1]) / get_denominator(ress[1, :, 1])))],
                    [errg[:,  4]           + epsilon, '$\\partial_y u_y$'          , (0.0, 1.0, 0.0), 'solid'  , np.log2(np.average(np.abs(ress[2, :, 1] - ress_bad[2, :, 1]) / get_denominator(ress[2, :, 1])))],
                    [errg[:,  5]           + epsilon, '$\\partial_z u_y$'          , (0.0, 0.7, 0.5), 'solid'  , np.log2(np.average(np.abs(ress[3, :, 1] - ress_bad[3, :, 1]) / get_denominator(ress[3, :, 1])))],
                    [errg[:,  6]           + epsilon, '$\\partial_x u_z$'          , (0.5, 0.0, 0.7), 'solid'  , np.log2(np.average(np.abs(ress[1, :, 2] - ress_bad[1, :, 2]) / get_denominator(ress[1, :, 2])))],
                    [errg[:,  7]           + epsilon, '$\\partial_y u_z$'          , (0.0, 0.5, 0.7), 'solid'  , np.log2(np.average(np.abs(ress[2, :, 2] - ress_bad[2, :, 2]) / get_denominator(ress[2, :, 2])))],
                    [errg[:,  8]           + epsilon, '$\\partial_z u_z$'          , (0.0, 0.0, 1.0), 'solid'  , np.log2(np.average(np.abs(ress[3, :, 2] - ress_bad[3, :, 2]) / get_denominator(ress[3, :, 2])))],
                    [errh[:,  0]           + epsilon, '$\\partial_x \\partial_x p$', (1.0, 0.0, 0.0), 'dashed' , np.log2(np.average(np.abs(ress[4, :, 3] - ress_bad[4, :, 3]) / get_denominator(ress[4, :, 3])))],
                    [errh[:,  1]           + epsilon, '$\\partial_x \\partial_y p$', (0.5, 0.5, 0.0), 'dashed' , np.log2(np.average(np.abs(ress[5, :, 3] - ress_bad[5, :, 3]) / get_denominator(ress[5, :, 3])))],
                    [errh[:,  2]           + epsilon, '$\\partial_x \\partial_z p$', (0.5, 0.0, 0.5), 'dashed' , np.log2(np.average(np.abs(ress[6, :, 3] - ress_bad[6, :, 3]) / get_denominator(ress[6, :, 3])))],
                    [errh[:,  3]           + epsilon, '$\\partial_y \\partial_y p$', (0.0, 1.0, 0.0), 'dashed' , np.log2(np.average(np.abs(ress[7, :, 3] - ress_bad[7, :, 3]) / get_denominator(ress[7, :, 3])))],
                    [errh[:,  4]           + epsilon, '$\\partial_y \\partial_z p$', (0.0, 0.5, 0.5), 'dashed' , np.log2(np.average(np.abs(ress[8, :, 3] - ress_bad[8, :, 3]) / get_denominator(ress[8, :, 3])))],
                    [errh[:,  5]           + epsilon, '$\\partial_z \\partial_z p$', (0.0, 0.0, 1.0), 'dashed' , np.log2(np.average(np.abs(ress[9, :, 3] - ress_bad[9, :, 3]) / get_denominator(ress[9, :, 3])))]]
        nbins = 15
        for ed in err_data[:3]:
            hist, bins = np.histogram(
                    np.log2(ed[0]),
                    bins = nbins)
            ax.plot(bins[1:],
                    np.cumsum(hist)/(3. * res4.shape[0]),
                    label = (ed[1]
                          + ', with $\\log_2 \\langle \\varepsilon \\rangle = '
                          + '{0:.2f}$'.format(np.log2(np.average(ed[0])))
                          + ' and $\\log_2 \\langle \\varepsilon_s \\rangle < {0:.2f}$'.format(ed[4])),
                    color = ed[2],
                    linestyle = ed[3])
        for ed in err_data[3:]:
            hist, bins = np.histogram(
                    np.log2(ed[0]),
                    bins = nbins)
            ax.plot(bins[1:],
                    np.cumsum(hist) / (1. * res4.shape[0]),
                    label = (ed[1]
                          + ', with $\\log_2 \\langle \\varepsilon \\rangle = '
                          + '{0:.2f}$'.format(np.log2(np.average(ed[0])))
                          + ' and $\\log_2 \\langle \\varepsilon_s \\rangle < {0:.2f}$'.format(ed[4])),
                    color = ed[2],
                    linestyle = ed[3])
        ax.set_xlim(-20, 16)
        ax.xaxis.set_ticks(range(-20, 8), minor = True)
        ax.set_xlabel('$\\log_2 \\varepsilon$')
        ax.yaxis.set_ticks([.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])
        ax.yaxis.set_ticks([.05, .15, .25, .35, .45, .55, .65, .75, .85, .95], minor = True)
        ax.legend(loc = 'best', fontsize=10)
        ax.set_title('cum histogram of $\\log_2 \\varepsilon$ over {0} points'.format(res4.shape[0]))
        fig.savefig(dir_name + 'err_Lag_histograms.pdf', format='pdf')
        fig.clear()
    return None

def estimate_interpolation_error(
        lTDB,
        info,
        npoints = 2**5,
        randseeds = range(2),
        dir_name = 'test',
        xfixed = None,
        yfixed = None,
        zfixed = None,
        nmax = 7,
        mmax = 2):
    generate_interpolation_points(
            info,
            npoints = npoints,
            randseeds = randseeds,
            dir_name = dir_name,
            xfixed = xfixed,
            yfixed = yfixed,
            zfixed = zfixed)
    spline_error(
            lTDB, info,
            nmax = nmax,
            mmax = mmax,
            npoints = npoints,
            randseeds = randseeds,
            dir_name = dir_name,
            plots_on = True)
    Lagrange_err_histograms(
            lTDB, info,
            npoints = npoints,
            spline_neighbours = 7,
            spline_continuity = 2,
            randseeds = randseeds,
            dir_name = dir_name,
            plots_on = True)
    return None

def test_interp(
        nbatches = 4,
        npoints = 2**5,
        isotropic1024coarse = False,
        mhd1024 = False,
        channel = False,
        nmax = 7,
        mmax = 2):
    # load shared library
    lJHTDB = pyJHTDB.libJHTDB()
    #initialize webservices
    lJHTDB.initialize()

    if isotropic1024coarse:
        estimate_interpolation_error(
                lJHTDB,
                pyJHTDB.dbinfo.isotropic1024coarse,
                npoints = npoints,
                randseeds = range(nbatches),
                nmax = nmax,
                mmax = mmax)
    if mhd1024:
        estimate_interpolation_error(
                lJHTDB,
                pyJHTDB.dbinfo.mhd1024,
                npoints = npoints,
                randseeds = range(nbatches),
                nmax = nmax,
                mmax = mmax)
    if channel:
        estimate_interpolation_error(
                lJHTDB,
                pyJHTDB.dbinfo.channel,
                npoints = npoints,
                randseeds = range(nbatches),
                dir_name = 'random',
                nmax = nmax,
                mmax = mmax)
        for ynode in [0, 64, 128, 192, 256]:
            yval = pyJHTDB.dbinfo.channel['ynodes'][ynode] + pyJHTDB.dbinfo.channel['dy'][ynode]/2
            estimate_interpolation_error(
                    lJHTDB,
                    pyJHTDB.dbinfo.channel,
                    npoints = npoints,
                    randseeds = range(nbatches),
                    dir_name = 'ynode_{0:0>3x}'.format(ynode),
                    yfixed = yval,
                    nmax = nmax,
                    mmax = mmax)

    #finalize webservices
    lJHTDB.finalize()
    return None

