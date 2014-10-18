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

import numpy
import scipy
import scipy.spatial

def points_on_sphere(
        N,
        origin = numpy.zeros(3),
        radius = 1.):
    """ Generate N evenly distributed points on the unit sphere centered at
        the origin. Uses the 'Golden Spiral'.
        Code by Chris Colbert from the numpy-discussion list.
    """
    phi = (1 + numpy.sqrt(5)) / 2   # the golden ratio
    long_incr = 2*numpy.pi / phi    # how much to increment the longitude

    dz = 2.0 / float(N)             # a unit sphere has diameter 2
    bands = numpy.arange(N)         # each band will have one point placed on it
    z = bands * dz - 1 + (dz/2)     # the height z of each band/point
    r = numpy.sqrt(1 - z*z)         # project onto xy-plane
    az = bands * long_incr          # azimuthal angle of point modulo 2 pi
    x = r * numpy.cos(az)
    y = r * numpy.sin(az)
    ## get triangles
    points = numpy.array([x, y, z])
    tri = scipy.spatial.ConvexHull(points.T)
    points = origin[None, :] + points.T*radius
    return points, tri.simplices

