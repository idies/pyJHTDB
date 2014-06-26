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

