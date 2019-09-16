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
#      on the hdf5 library, since that can't be installed through pip.
#
import pyJHTDB
VERSION = pyJHTDB.version
#if HDF5_ON:
#    VERSION += '-hdf5'
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



########################################################################
#
# generate MANIFEST.in
#
open('MANIFEST.in',
     'w').write(
             'include turblib/turblib.c\n' +
             'include turblib/soapC.c\n' +
             'include turblib/soapClient.c\n' +
             'include turblib/stdsoap2.c\n' +
             'include turblib/*.h\n' +
             'include turblib/*.nsmap')
#
########################################################################



libraries = []
macros = []
#if h5cc_present and HDF5_ON:
#    libraries.append('hdf5')
#    macros.append(('CUTOUT_SUPPORT', '1'))

from setuptools import setup, Extension
libJHTDB = Extension(
        'libJHTDB',
        sources = ['C/local_tools.c',
                   'turblib/turblib.c',
                   'turblib/soapC.c',
                   'turblib/soapClient.c',
                   'turblib/stdsoap2.c'],
        include_dirs = ['turblib'],
        define_macros = macros,
        libraries = libraries)

setup(
        name = 'pyJHTDB',
        version = VERSION,
        packages = ['pyJHTDB'],
        ###scripts=['scripts/pyJHTDB-auth'],
        package_data = {'pyJHTDB': ['data/channel_xgrid.npy',
                                    'data/channel_ygrid.npy',
                                    'data/channel_zgrid.npy',
                                    'data/channel5200_ygrid.npy',
                                    'data/transition_bl_ygrid.npy']},
        python_requires='>3.6',
        install_requires = ['numpy>=1.15.0', 'scipy>=1.1.0', 'sympy>=1.2', 'h5py>=2.8.0', 'matplotlib>=3.0.0'],
        ext_modules = [libJHTDB],
        test_suite = 'tests',

        #### package description stuff goes here
        description = 'Python wrapper for the Johns Hopkins turbulence database library',
        long_description = open('README.md', 'r').read(),
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        license = 'Apache Version 2.0',
        url = GROUP_URL,
        download_url = 'https://github.com/idies/pyJHTDB',
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

