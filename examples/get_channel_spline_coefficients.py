import pyJHTDB

# M1Q4
ii = pyJHTDB.interpolator.spline_interpolator(pyJHTDB.dbinfo.channel5200)
ii.write_coefficients()

# M2Q8
ii = pyJHTDB.interpolator.spline_interpolator(pyJHTDB.dbinfo.channel5200, m = 2, n = 3)
ii.write_coefficients()


# M2Q14
ii = pyJHTDB.interpolator.spline_interpolator(pyJHTDB.dbinfo.channel5200, m = 2, n = 6)
ii.write_coefficients()

