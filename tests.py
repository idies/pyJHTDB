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

import sys
sys.path[0] = ''

import argparse

parser = argparse.ArgumentParser(
    description = 'Test pyJHTDB installation.')
parser.add_argument(
    '-p',
    '--plain',
    dest = 'plain',
    action = 'store_true',
    help = 'run plain test, i.e. turbc clone.')
parser.add_argument(
    '--grid-splines',
    dest = 'grid_splines',
    action = 'store_true',
    help = 'run basic grid spline test.')
parser.add_argument(
    '--cutout',
    dest = 'cutout',
    action = 'store_true',
    help = 'run cutout test.')
parser.add_argument(
    '--misc',
    dest = 'misc',
    action = 'store_true',
    help = 'run misc test.')
parser.add_argument(
    '--interpolator',
    dest = 'interpolator',
    action = 'store_true',
    help = 'run interpolator test.')

opt = parser.parse_args()

import pyJHTDB

if opt.plain:
    pyJHTDB.test_plain()
if opt.grid_splines:
    pyJHTDB.test_gs()
if opt.interpolator:
    pyJHTDB.test_interpolator()

if opt.misc and pyJHTDB.found_matplotlib:
    pyJHTDB.test_misc()
if opt.cutout and pyJHTDB.found_h5py:
    pyJHTDB.test_cutout()

