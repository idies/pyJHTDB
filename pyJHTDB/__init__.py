"""
    Python tools and wrappers for the Johns Hopkins Turbulence Databases C library.
"""

import os
import sys
import numpy as np
import h5py
import ctypes
import inspect
import platform

from pkg_resources import get_distribution, DistributionNotFound
import os.path

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

auth_token = 'edu.jhu.pha.turbulence.testing-201302'

#__configuration_file__ = os.path.join(os.path.expanduser('~'), '/.config/pyJHTDB.cfg')
#__config__ = ConfigParser.ConfigParser()
#__config__.readfp(open(__configuration_file__))

from libJHTDB import *
