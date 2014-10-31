from pyJHTDB.interpolator import spline_interpolator as si
from pyJHTDB.dbinfo import isotropic1024coarse, mhd1024, channel

def write_coefficients(
        info = None,
        n = 1,
        m = 1):
    i = si(info, n = n, m = m)
    for coord in ['x', 'y', 'z']:
        for order in range(i.m+1):
            text_file = open(
                    (info['name']
                    + '_' + coord
                    + 'spline_m{0}q{1}_d{2}_coeff.txt'.format(i.m, i.n*2 + 2, order)),
                    'w')
            if info[coord + 'periodic']:
                for point in range(len(i.spline[coord].beta[0][order])):
                    text_file.write('0, {0}'.format(i.spline[coord].neighbour_list[0][point]))
                    for c in i.spline[coord].beta[0][order][point].coef:
                        text_file.write(', {0}'.format(c))
                    text_file.write('\n')
            else:
                for node in range(len(i.spline[coord].beta)):
                    for point in range(len(i.spline[coord].beta[node][order])):
                        if (i.spline[coord].beta[node][order][point].coef.shape[0] > 1
                             or (not( i.spline[coord].beta[node][order][point].coef[0] == 0.0))):
                            text_file.write('{0}, {1}'.format(node, i.spline[coord].neighbour_list[node][point]))
                            for c in i.spline[coord].beta[node][order][point].coef:
                                text_file.write(', {0}'.format(c))
                            text_file.write('\n')
            text_file.close()
    return None

def main():
    write_coefficients(info = channel)
    write_coefficients(info = channel, n = 2)
    write_coefficients(info = channel, n = 3, m = 2)
    return None

if __name__ == '__main__':
    main()

