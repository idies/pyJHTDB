import numpy as np
import sympy as sp
import copy

import matplotlib.pyplot as plt

def get_fornberg_coeffs(
        x = None,
        a = None):
    N = len(a) - 1
    d = []
    for m in range(N+1):
        d.append([])
        for n in range(N+1):
            d[m].append([])
            for j in range(N+1):
                d[m][n].append(sp.Rational(0))
    d[0][0][0] = sp.Rational(1)
    c1 = sp.Rational(1)
    for n in range(1, N+1):
        c2 = sp.Rational(1)
        for j in range(n):
            c3 = a[n] - a[j]
            c2 = c2*c3
            for m in range(n+1):
                d[m][n][j] = ((a[n] - x)*d[m][n-1][j] - m*d[m-1][n-1][j]) / c3
        for m in range(n+1):
            d[m][n][n] = (c1 / c2)*(m*d[m-1][n-1][n-1] - (a[n-1] - x)*d[m][n-1][n-1])
        c1 = c2
    coeffs = []
    for m in range(len(d)):
        coeffs.append([])
        for j in range(len(d)):
            coeffs[-1].append(d[m][N][j])
    return sp.Matrix(coeffs)

def get_alpha_polynomials(
        max_deriv = 2):
    alpha = []
    xi = sp.Symbol('xi')
    for l in range(max_deriv + 1):
        alpha.append(
                xi**l / sp.factorial(l)
              * (1 - xi)**(max_deriv + 1)
              * sum(sp.factorial(max_deriv + k) * xi**k / (sp.factorial(max_deriv)*sp.factorial(k))
                    for k in range(max_deriv - l + 1)))
    return (xi, sp.Matrix(alpha))

class uniform_spline_1D:
    def __init__(
            self,
            xvals,
            periodic = False,
            max_deriv = 1,
            neighbours = 1):
        self.x = xvals.copy()
        self.dx = self.x[1] - self.x[0]
        self.m = max_deriv
        self.n = neighbours
        self.N = 2*neighbours + 2
        self.periodic = periodic
        self.deriv_coeff = []
        self.beta = []
        self.xi, self.alpha0 = get_alpha_polynomials(max_deriv = self.m)
        if not self.periodic:
            prev_x = np.array([self.x[0] - (k + 1)*self.dx
                               for k in range(self.n)])[::-1]
            post_x = np.array([self.x[self.x.shape[0] - 1] + (k + 1)*self.dx
                               for k in range(self.n)])
        else:
            self.period = self.x[-1] - self.x[0] + self.dx
            prev_x = self.x[-self.n:] - self.period
            post_x = self.x[:self.n+1] + self.period
        self.tmpx = np.zeros((self.x.shape[0] + self.n + post_x.shape[0]), dtype = self.x.dtype)
        self.tmpx[:self.n] = prev_x[:]
        self.tmpx[self.n:self.n + self.x.shape[0]] = self.x[:]
        self.tmpx[self.n + self.x.shape[0]:] = post_x[:]
        return None
    def put_yvals(self, yvals):
        self.y = yvals.copy()
        if not self.periodic:
            tmpderiv = (self.y[1] - self.y[0])
            prev_y = np.array([self.y[0] - (k + 1)*tmpderiv
                               for k in range(self.n)])
            tmpderiv = (self.y[self.y.shape[0] - 1] - self.y[self.y.shape[0] - 2])
            post_y = np.array([self.y[self.y.shape[0] - 1] + (k + 1)*tmpderiv
                               for k in range(self.n)])
        else:
            prev_y = self.y[-self.n:] 
            post_y = self.y[: self.n+1]
        self.tmpy = np.zeros((self.y.shape[0] + self.n + post_y.shape[0]), dtype = self.y.dtype)
        self.tmpy[:self.n] = prev_y[:]
        self.tmpy[self.n:self.n + self.y.shape[0]] = self.y[:]
        self.tmpy[self.n + self.y.shape[0]:] = post_y[:]
        return None
    def __call__(self, x):
        if not self.periodic:
            ix = np.searchsorted(self.x, x) - 1
            if ix < 0:
                return self.y[0]
            elif ix >=  self.x.shape[0] - 1:
                return self.y[self.x.shape[0] - 1]
            xi = (x - self.x[ix]) / self.dx
            return (sum(self.beta[k]*self.tmpy[ix + k]
                    for k in range(self.N))).subs(sp.Symbol('xi'), xi)
        else:
            x = np.remainder(x, self.period)
            ix = np.searchsorted(self.tmpx, x) - 1
            xi = (x - self.tmpx[ix]) / self.dx
            return (sum(self.beta[k]*self.tmpy[(ix -self.n+ k)%self.tmpy.shape[0]]
                    for k in range(self.N))).subs(sp.Symbol('xi'), xi)
    def compute(self, x):
        return np.array([self.__call__(x[i]) for i in range(x.shape[0])])
    def compute_derivs(self):
        x0 = sp.Symbol('cd_temp_x0')
        a = [sp.Symbol('cd_temp_alpha_{0}'.format(k)) for k in range(self.N-1)]
        c0 = get_fornberg_coeffs(x0, a)
        ctmp = copy.deepcopy(c0)
        ctmp = ctmp.subs(x0, self.x[0])
        for k in range(self.N - 1):
            ctmp = ctmp.subs(a[k], self.tmpx[k])
        self.deriv_coeff = ctmp
        return None
    def compute_beta(self):
        btmp = [sum(self.deriv_coeff[l, 0]*self.alpha0[l]*self.dx**l
                    for l in range(self.m + 1))]
        for k in range(1, self.N-1):
            btmp.append(sum(self.deriv_coeff[l, k  ]*self.alpha0[l]*self.dx**l
                          + self.deriv_coeff[l, k-1]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx)**l
                            for l in range(self.m + 1)))
        btmp.append(sum(self.deriv_coeff[l, self.N-2]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx)**l
                        for l in range(self.m + 1)))
        self.beta = sp.Matrix(btmp)
        return None

class generic_spline_1D:
    def __init__(
            self,
            xvals,
            periodic = False,
            max_deriv = 1,
            neighbours = 1):
        self.x = xvals.copy()
        self.dx = self.x[1:] - self.x[:self.x.shape[0] - 1]
        self.m = max_deriv
        self.n = neighbours
        self.N = 2*neighbours + 2
        self.periodic = periodic
        self.deriv_coeff = []
        self.beta = []
        self.xi, self.alpha0 = get_alpha_polynomials(max_deriv = self.m)
        if not self.periodic:
            prev_x = np.array([self.x[0] - (k + 1)*self.dx[0]
                               for k in range(self.n)])[::-1]
            post_x = np.array([self.x[self.x.shape[0] - 1] + (k + 1)*self.dx[self.dx.shape[0] - 1]
                               for k in range(self.n)])
        else:
            self.dx = np.append(self.dx, [self.dx[0]])
            self.period = self.x[-1] - self.x[0] + self.dx[0]
            prev_x = self.x[-self.n:] - self.period
            post_x = self.x[:self.n+1] + self.period
        self.tmpx = np.zeros((self.x.shape[0] + self.n + post_x.shape[0]), dtype = self.x.dtype)
        self.tmpx[:self.n] = prev_x[:]
        self.tmpx[self.n:self.n + self.x.shape[0]] = self.x[:]
        self.tmpx[self.n + self.x.shape[0]:] = post_x[:]
        return None
    def put_yvals(self, yvals):
        self.y = yvals.copy()
        if not self.periodic:
            tmpderiv = (self.y[1] - self.y[0])
            prev_y = np.array([self.y[0] - (k + 1)*tmpderiv
                               for k in range(self.n)])
            tmpderiv = (self.y[self.y.shape[0] - 1] - self.y[self.y.shape[0] - 2])
            post_y = np.array([self.y[self.y.shape[0] - 1] + (k + 1)*tmpderiv
                               for k in range(self.n)])
        else:
            prev_y = self.y[-self.n:] 
            post_y = self.y[: self.n+1]
        self.tmpy = np.zeros((self.y.shape[0] + self.n + post_y.shape[0]), dtype = self.y.dtype)
        self.tmpy[:self.n] = prev_y[:]
        self.tmpy[self.n:self.n + self.y.shape[0]] = self.y[:]
        self.tmpy[self.n + self.y.shape[0]:] = post_y[:]
        return None
    def __call__(self, x):
        if not self.periodic:
            ix = np.searchsorted(self.x, x) - 1
            if ix < 0:
                return self.y[0]
            elif ix >=  self.x.shape[0] - 1:
                return self.y[self.x.shape[0] - 1]
            xi = (x - self.x[ix]) / self.dx[ix]
            return (sum(self.beta[ix][k]*self.tmpy[ix + k]
                    for k in range(self.N))).subs(sp.Symbol('xi'), xi)
        else:
            x = np.remainder(x, self.period)
            ix = np.searchsorted(self.tmpx, x) - 1
            xi = (x - self.tmpx[ix]) / self.dx[(ix-self.n)%self.dx.shape[0]]
            return (sum(self.beta[(ix-self.n)%len(self.beta)][k]*self.tmpy[(ix -self.n+ k)%self.tmpy.shape[0]]
                    for k in range(self.N))).subs(sp.Symbol('xi'), xi)
    def compute(self, x):
        return np.array([self.__call__(x[i]) for i in range(x.shape[0])])
    def compute_derivs(self):
        x0 = sp.Symbol('cd_temp_x0')
        a = [sp.Symbol('cd_temp_alpha_{0}'.format(k)) for k in range(self.N-1)]
        c0 = get_fornberg_coeffs(x0, a)
        for i in range(self.x.shape[0]):
            ctmp = copy.deepcopy(c0)
            ctmp = ctmp.subs(x0, self.x[i])
            for k in range(self.N - 1):
                ctmp = ctmp.subs(a[k], self.tmpx[i+k])
            self.deriv_coeff.append(ctmp)
        if self.periodic:
            i = self.x.shape[0]
            ctmp = copy.deepcopy(c0)
            ctmp = ctmp.subs(x0, self.tmpx[i+self.n])
            for k in range(self.N - 1):
                ctmp = ctmp.subs(a[k], self.tmpx[i+k])
            self.deriv_coeff.append(ctmp)
        return None
    def compute_beta(self):
        for i in range(self.x.shape[0]-1):
            btmp = [sum(self.deriv_coeff[i][l, 0]*self.alpha0[l]*self.dx[i]**l
                        for l in range(self.m + 1))]
            for k in range(1, self.N-1):
                btmp.append(sum(self.deriv_coeff[i  ][l, k  ]*self.alpha0[l]*self.dx[i]**l
                              + self.deriv_coeff[i+1][l, k-1]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx[i])**l
                                for l in range(self.m + 1)))
            btmp.append(sum(self.deriv_coeff[i+1][l, self.N-2]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx[i])**l
                            for l in range(self.m + 1)))
            self.beta.append(sp.Matrix(btmp))
        if self.periodic:
            i = self.x.shape[0] - 1
            btmp = [sum(self.deriv_coeff[i][l, 0]*self.alpha0[l]*self.dx[i]**l
                        for l in range(self.m + 1))]
            for k in range(1, self.N-1):
                btmp.append(sum(self.deriv_coeff[i  ][l, k  ]*self.alpha0[l]*self.dx[i]**l
                              + self.deriv_coeff[i+1][l, k-1]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx[i])**l
                                for l in range(self.m + 1)))
            btmp.append(sum(self.deriv_coeff[i+1][l, self.N-2]*self.alpha0[l].subs(self.xi, 1 - self.xi)*(-self.dx[i])**l
                            for l in range(self.m + 1)))
            self.beta.append(sp.Matrix(btmp))
        return None

def plot_uniform_weight_functions(
        n = 4,
        m = 2):
    x = np.arange(-n, n+1, 1)
    tst0 = uniform_spline_1D(
            x,
            max_deriv = m,
            neighbours = n)
    tst0.compute_derivs()
    tst0.compute_beta()
    xval = np.arange(x[0], x[-1], .05)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.set_title('Weight functions for {0} neighbours and {1} continuous derivatives'.format(n, m))
    for i in range(n+1):
        y = np.zeros(x.shape, x.dtype)
        y[i] = 1
        tst0.put_yvals(y)
        f = tst0.compute(xval)
        ax.plot(xval, f)
    fig.savefig('test.pdf', format = 'pdf')
    return None

def plot_generic_weight_functions(
        n = 4,
        m = 2):
    x = np.random.random(2*n + 1)
    x.sort()
    tst0 = generic_spline_1D(
            x,
            max_deriv = m,
            neighbours = n)
    tst0.compute_derivs()
    tst0.compute_beta()
    xval = []
    for i in range(x.shape[0]-1):
        xtmp = [x[i] + k*.1*(x[i+1] - x[i])
                         for k in range(10)]
        xval += xtmp
    xval = np.array(xval)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.set_title('Weight functions for {0} neighbours and {1} continuous derivatives'.format(n, m))
    for i in range(n+1):
        y = np.zeros(x.shape, x.dtype)
        y[i] = 1
        tst0.put_yvals(y)
        f = tst0.compute(xval)
        ax.plot(xval, f)
    fig.savefig('test.pdf', format = 'pdf')
    return None

def main0():
    plot_generic_weight_functions(n = 4, m = 2)
    return None

if __name__ == '__main__':
    main0()

