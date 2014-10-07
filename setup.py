#! /usr/bin/env python

import shlib
from shlib import build_shlib

from distutils.core import setup
#from setuptools import setup # doesn't work because of build_shlib...
import datetime

TURBLIB_VERSION = '20140606'

########################################################################
#
# version
#
now = datetime.datetime.now()
date_name = '{0:0>4}{1:0>2}{2:0>2}'.format(now.year, now.month, now.day)
#
########################################################################



########################################################################
#
# check what's available on the system
#
from distutils.spawn import find_executable
h5cc_executable = find_executable('h5cc')
h5cc_present = not (h5cc_executable == None)
#
########################################################################



########################################################################
#
# get the turbulence library
#
import os
if not os.path.isdir('turblib-' + TURBLIB_VERSION):
    import urllib
    bla = urllib.urlretrieve('http://turbulence.pha.jhu.edu/download/turblib-' + TURBLIB_VERSION + '.tar.gz',
                             'turblib-' + TURBLIB_VERSION + '.tar.gz')
    import tarfile
    turblib = tarfile.open('turblib-' + TURBLIB_VERSION + '.tar.gz')
    turblib.extractall()
    turblib.close()
#
########################################################################



libraries = []
extra_compile_args = []
if h5cc_present:
    libraries.append('hdf5')
    extra_compile_args.append('-DCUTOUT_SUPPORT')
setup(
        name = 'pyJHTDB',
        version = date_name,
        packages = ['pyJHTDB',],
        package_data = {'pyJHTDB': ['data/channel_xgrid.npy',
                                    'data/channel_ygrid.npy',
                                    'data/channel_zgrid.npy']},
        shlibs = [build_shlib.SharedLibrary(
                     name         = 'JHTDB',
                     sources      = ['C/local_tools.c',
                                     'turblib-' + TURBLIB_VERSION + '/turblib.c',
                                     'turblib-' + TURBLIB_VERSION + '/soapC.c',
                                     'turblib-' + TURBLIB_VERSION + '/soapClient.c',
                                     'turblib-' + TURBLIB_VERSION + '/stdsoap2.c'],
                     libraries    = libraries,
                     include_dirs = ['turblib-' + TURBLIB_VERSION],
                     language     = 'c',
                     extra_compile_args = extra_compile_args)],
        install_requires = 'numpy>=1.6')

