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

"""
Python tools and wrappers for the Johns Hopkins Turbulence Database
Cluster C library.
Contact: turbulence@pha.jhu.edu
Website: http://turbulence.pha.jhu.edu/

Although this particular Python wrapper is still a work in progress, it
is mature enough to be used in production work.

On first contact with this library, we recommend that you first run
"test_plain". To be more specific:

.. code:: python

    >>> from pyJHTDB import test_plain
    >>> test_plain()

The code that is executed can be found in "pyJHTDB/test.py", and it's
the simplest example of how to access the turbulence database.

"""

import os
import os.path
import sys
import numpy as np
import ctypes
import platform

from pkg_resources import get_distribution, DistributionNotFound

try:
    _dist = get_distribution('pyJHTDB')
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(_dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, 'pyJHTDB')):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install this project with setup.py'
else:
    __version__ = _dist.version

auth_token = 'edu.jhu.pha.turbulence.testing-201406'
homefolder = os.path.expanduser('~')
lib_folder = os.path.join(homefolder, '.config/', 'JHTDB/')
version_info = "20200827.0"
version = str(version_info)

# check if .config/JHTDB folder exists, create it if not
#if os.path.isdir(lib_folder):
#    if os.path.exists(os.path.join(lib_folder, 'auth_token.txt')):
#        auth_token = str(open(os.path.join(lib_folder, 'auth_token.txt'), 'r').readline().split()[0])
#    else:
#        open(os.path.join(lib_folder, 'auth_token.txt'), 'w').write(auth_token)
#else:
#    os.mkdir(lib_folder)
#    open(os.path.join(lib_folder, 'auth_token.txt'), 'w').write(auth_token)

try:
    import h5py
    found_h5py = True
except ImportError:
    found_h5py = False
    print('h5py not found. cutout functionality not available.')

try:
    import matplotlib
    found_matplotlib = True
except ImportError:
    found_matplotlib = False
    print('matplotlib not found. plotting functionality not available.')

try:
    import scipy
    found_scipy = True
except ImportError:
    found_scipy = False
    print('scipy not found. not all interpolation functionality available')

from .libJHTDB import *
from .test import test_plain
if pyJHTDB.found_matplotlib:
    from .test import test_misc
from .test import test_interp_1D as test_interpolator
from .generic_splines import main0 as test_gs
