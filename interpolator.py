import pyTDB
import pyTDB.generic_splines as gs

class data_interpolator:
    def __init__(
            self,
            info,
            n = 1,
            m = 1):
        if os.path.exists(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m)):
            return pickle.load(open(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m), 'r'))
        func = []
        for coord in ['x', 'y', 'z']:
            if info[coord + 'uniform']:
                func.append(generic_spline_1D(
                        info[coord + 'nodes'][:1],
                        max_deriv = m,
                        neighbours = n,
                        periodic = info[coord + 'periodic']))
            else:
                func.append(generic_spline_1D(
                        info[coord + 'nodes'],
                        max_deriv = m,
                        neighbours = n,
                        periodic = info[coord + 'periodic']))
            func[-1].compute_derivs()
            func[-1].compute_beta()
        interpolator = {'x': func[0],
                        'y': func[1],
                        'z': func[2]}
        pickle.dump(interpolator,
                    open(info['name'] + '_spline_interpolator_n{0}_m{1}.p'.format(n, m), 'w'))
        return interpolator
