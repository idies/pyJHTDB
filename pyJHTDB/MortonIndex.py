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

# taken from http://code.activestate.com/recipes/577558-interleave-bits-aka-morton-ize-aka-z-order-curve/

import numpy

def part1by2(x):
        n = x & 0x000003ff
        n = (n ^ (n << 16)) & 0xff0000ff
        n = (n ^ (n <<  8)) & 0x0300f00f
        n = (n ^ (n <<  4)) & 0x030c30c3
        n = (n ^ (n <<  2)) & 0x09249249
        return n

def unpart1by2(z):
        n = z & 0x09249249
        n = (n ^ (n >>  2)) & 0x030c30c3
        n = (n ^ (n >>  4)) & 0x0300f00f
        n = (n ^ (n >>  8)) & 0xff0000ff
        n = (n ^ (n >> 16)) & 0x000003ff
        return n

def grid3D_to_zindex(x):
    """return the corresponding Morton z-indices for an array of 3D indices

       :param x: input indices
       :type x: numpy array of integers, of shape (3, whatever),
          where n is the number of points
       :returns: z, array of integers, of shape (whatever)
    """
    return part1by2(x[0]) | (part1by2(x[1]) << 1) | (part1by2(x[2]) << 2)

def zindex_to_grid3D(z):
    """return the 3D indices corresponding to an array of Morton z-indices

       :param z: input indices
       :type z: numpy array of integers, of shape (whatever),
          where n is the number of points
       :returns: x, array of integers, of shape (3, whatever)
    """
    return numpy.array([unpart1by2(z     ),
                        unpart1by2(z >> 1),
                        unpart1by2(z >> 2)])
