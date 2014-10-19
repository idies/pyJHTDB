#! /usr/bin/env python
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


########################################################################
#
# some global settings
#
TURBLIB_VERSION = '20140606'
HDF5_ON = True
GROUP_URL = 'http://turbulence.pha.jhu.edu/'
GROUP_EMAIL = 'turbulence@pha.jhu.edu'
GROUP_NAME  = 'Johns Hopkins Turbulence Database Group'
AUTHOR = GROUP_NAME
AUTHOR_EMAIL = GROUP_EMAIL
#
########################################################################



########################################################################
#
# define version for pyJHTDB
# TODO:
#   1. VERSION should come from checkingwhen the sources were last
#      modified.
#   2. VERSION should contain information on whether or not it depends
#      on the hdf5 library.
#
import datetime
now = datetime.datetime.now()
date_name = '{0:0>4}{1:0>2}{2:0>2}'.format(now.year, now.month, now.day)
VERSION = date_name
#if HDF5_ON:
#    VERSION += '-hdf5'
#
########################################################################



########################################################################
#
# get the turbulence library
#
import os
if not os.path.isdir('turblib-' + TURBLIB_VERSION):
    import urllib
    urllib.urlretrieve('http://turbulence.pha.jhu.edu/download/turblib-'
                                + TURBLIB_VERSION
                                + '.tar.gz',
                             'turblib-'
                                + TURBLIB_VERSION
                                + '.tar.gz')
    import tarfile
    turblib = tarfile.open('turblib-' + TURBLIB_VERSION + '.tar.gz')
    turblib.extractall()
    turblib.close()
open('MANIFEST.in', 'w').write('graft turblib-' + TURBLIB_VERSION)
#
########################################################################



########################################################################
#
# check what's available on the system
#
import distutils.spawn
h5cc_executable = distutils.spawn.find_executable('h5cc')
h5cc_present = not (h5cc_executable == None)
#
########################################################################



libraries = []
macros = []
if h5cc_present and HDF5_ON:
    libraries.append('hdf5')
    macros.append(('CUTOUT_SUPPORT', '1'))

from setuptools import setup, Extension
libJHTDB = Extension(
        'libJHTDB',
        sources = ['C/local_tools.c',
                   'turblib-' + TURBLIB_VERSION + '/turblib.c',
                   'turblib-' + TURBLIB_VERSION + '/soapC.c',
                   'turblib-' + TURBLIB_VERSION + '/soapClient.c',
                   'turblib-' + TURBLIB_VERSION + '/stdsoap2.c'],
        include_dirs = ['turblib-' + TURBLIB_VERSION],
        define_macros = macros,
        libraries = libraries)

setup(
        name = 'pyJHTDB',
        version = VERSION,
        packages = ['pyJHTDB'],
        package_data = {'pyJHTDB': ['data/channel_xgrid.npy',
                                    'data/channel_ygrid.npy',
                                    'data/channel_zgrid.npy']},
        install_requires = ['numpy>=1.8', 'sympy>=0.7.4.1'],
        ext_modules = [libJHTDB],
        test_suite = 'tests',

        #### package description stuff goes here
        description = 'Python wrapper for the Johns Hopkins turbulence database library',
        long_description = open('README.rst', 'r').read(),
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        license = 'Apache Version 2.0',
        url = GROUP_URL,
        download_url = 'https://github.com/chichilalescu/pyJHTDB',
        classifiers = [
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Natural Language :: English',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics',
            ],
        )

