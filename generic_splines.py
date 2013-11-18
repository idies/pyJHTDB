import numpy as np
import sympy as sp
import copy
import pickle
import os
import sys

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

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
    return np.array(coeffs)

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

class generic_spline_1D:
    def __init__(
            self,
            xvals,
            period = None,
            max_deriv = 1,
            neighbours = 1):
        self.x = xvals.copy()
        self.dx = self.x[1:] - self.x[:self.x.shape[0] - 1]
        self.m = max_deriv
        self.n = neighbours
        self.N = 2*neighbours + 2
        self.periodic = (not period == None)
        self.deriv_coeff = []
        self.beta = []
        self.xi, self.alpha0 = get_alpha_polynomials(max_deriv = self.m)
        self.alpha0_coeff = []
        self.alpha1_coeff = []
        for l in range(self.m + 1):
            tcoeff0 = sp.Poly(self.alpha0[l], self.xi).all_coeffs()
            tcoeff1 = sp.Poly(self.alpha0[l].subs(self.xi, 1 - self.xi)*(-1)**l, self.xi).all_coeffs()
            tcoeff0.reverse()
            tcoeff1.reverse()
            self.alpha0_coeff.append(tcoeff0)
            self.alpha1_coeff.append(tcoeff1)
        self.alpha0_coeff = np.array(self.alpha0_coeff)
        self.alpha1_coeff = np.array(self.alpha1_coeff)
        if not self.periodic:
            prev_x = np.array([self.x[0] - (k + 1)*self.dx[0]
                               for k in range(self.n)])[::-1]
            post_x = np.array([self.x[self.x.shape[0] - 1] + (k + 1)*self.dx[self.dx.shape[0] - 1]
                               for k in range(self.n)])
            self.tmpx = np.zeros((self.x.shape[0] + self.n + post_x.shape[0]), dtype = self.x.dtype)
            self.tmpx[:self.n] = prev_x[:]
            self.tmpx[self.n:self.n + self.x.shape[0]] = self.x[:]
            self.tmpx[self.n + self.x.shape[0]:] = post_x[:]
        else:
            self.period = period
            self.tmpx = np.arange(-self.n, self.n+3, 1)*self.dx[0]
            self.dx = np.array([self.dx[0], self.dx[0]])
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
        shape_list = [self.y.shape[0] + self.n + post_y.shape[0]]
        for i in range(1, len(self.y.shape)):
            shape_list.append(self.y.shape[i])
        self.yshape = tuple(shape_list[1:])
        self.tmpy = np.zeros(tuple(shape_list), dtype = self.y.dtype)
        self.tmpy[:self.n] = prev_y[:]
        self.tmpy[self.n:self.n + self.y.shape[0]] = self.y[:]
        self.tmpy[self.n + self.y.shape[0]:] = post_y[:]
        return None
    def __call__(self, x, order = 0):
        if not self.periodic:
            ix = np.searchsorted(self.x, x) - 1
            if ix < 0:
                return self.y[0]
            elif ix >=  self.x.shape[0] - 1:
                return self.y[self.x.shape[0] - 1]
            xi = (x - self.x[ix]) / self.dx[ix]
            return (sum(self.fast_beta[ix][order][k](xi)*self.tmpy[ix + k]
                    for k in range(self.N)))
        else:
            x = np.remainder(x, self.period)
            ix = np.searchsorted(self.tmpx, x) - 1
            xi = (x - self.tmpx[ix]) / self.dx[(ix-self.n)%self.dx.shape[0]]
            return sum(self.fast_beta[(ix-self.n)%len(self.fast_beta)][order][k](xi)
                      *self.tmpy[(ix-self.n+k)%self.tmpy.shape[0]]
                       for k in range(self.N))
    def beta_values(self, xfrac = 0, xgrid = 0, order = 0):
        return np.array([self.fast_beta[xgrid][order][k](xfrac) for k in range(self.N)])
    def compute_derivs(self):
        for i in range(self.x.shape[0]):
            self.deriv_coeff.append(get_fornberg_coeffs(self.x[i], self.tmpx[i:i+self.N-1]))
        if self.periodic:
            i = self.x.shape[0]
            self.deriv_coeff.append(get_fornberg_coeffs(self.tmpx[i], self.tmpx[i:i+self.N-1]))
        return None
    def compute_beta(self):
        for i in range(len(self.deriv_coeff)-1):
            deltax = np.array([self.dx[i]**l for l in range(self.m + 1)])
            a0 = self.alpha0_coeff*deltax[:, np.newaxis]
            a1 = self.alpha1_coeff*deltax[:, np.newaxis]
            btmp = [np.polynomial.polynomial.Polynomial(
                        list(np.sum(self.deriv_coeff[i][:self.m+1, 0, np.newaxis]*a0, axis = 0)))]
            for k in range(1, self.N-1):
                btmp.append(np.polynomial.polynomial.Polynomial(np.sum(
                    self.deriv_coeff[i  ][:self.m+1, k  , np.newaxis]*a0
                  + self.deriv_coeff[i+1][:self.m+1, k-1, np.newaxis]*a1, axis = 0)))
            btmp.append(np.polynomial.polynomial.Polynomial(list(
                    np.sum(self.deriv_coeff[i+1][:self.m+1, self.N-2, np.newaxis]*a1, axis = 0))))
            diff_parameter = [tuple([0 for j in range(l)]) for l in range(self.m+1)]
            self.beta.append([[btmp[k].deriv(j)*self.dx[i]**(-j)
                               for k in range(self.N)]
                              for j in range(self.m+1)])
        return None
    def compute_fast_beta(self):
        self.fast_beta = []
        for i in range(len(self.beta)):
            self.fast_beta.append([[sp.utilities.lambdify((self.xi),
                sp.horner(sp.Poly((self.beta[i][j][k].coef[::-1]), self.xi)), np)
                                    for k in range(self.N)]
                                   for j in range(self.m + 1)])
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

