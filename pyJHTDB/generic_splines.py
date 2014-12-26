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

import numpy as np
import sympy as sp
import copy
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
        self.m = int(np.floor(max_deriv))
        self.n = int(np.floor(neighbours))
        self.N = 2*self.n + 2
        self.periodic = not (period == None)
        self.uniform = (self.x.shape[0] == 2)
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
        if self.periodic:
            self.period = period
            if self.uniform:
                self.tmpx = np.arange(-self.n, self.n+3, 1)*self.dx[0]
                self.dx = np.array([self.dx[0], self.dx[0]])
            else:
                prev_x = self.x[self.x.shape[0]-self.n:] - period
                post_x = self.x[:self.n+1] + period
                self.tmpx = np.zeros((self.x.shape[0] + self.n + post_x.shape[0]), dtype = self.x.dtype)
                self.tmpx[:self.n] = prev_x[:]
                self.tmpx[self.n:self.n + self.x.shape[0]] = self.x[:]
                self.tmpx[self.n + self.x.shape[0]:] = post_x[:]
                self.dx = np.append(self.dx, self.x[0] + period - self.x[-1])
        return None
    def put_yvals(self, yvals):
        self.y = yvals.copy()
        if self.periodic:
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
        else:
            self.yshape = self.y.shape[1:]
        return None
    def __call__(self, x, order = 0):
        if not self.periodic:
            ix = np.searchsorted(self.x, x) - 1
            if ix < 0:
                return self.y[0]
            elif ix >=  self.x.shape[0] - 1:
                return self.y[self.x.shape[0] - 1]
            xi = (x - self.x[ix]) / self.dx[ix]
            if ix < self.n:
                return sum(self.fast_beta[ix][order][k](xi)*self.y[k]
                           for k in range(self.N-1))
            elif ix >= self.x.shape[0] - self.n - 1:
                return sum(self.fast_beta[ix][order][k](xi)*self.y[self.x.shape[0] - self.N + k+1]
                           for k in range(self.N-1))
            return sum(self.fast_beta[ix][order][k](xi)*self.y[ix - self.n + k]
                       for k in range(self.N))
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
        if self.periodic:
            for i in range(self.x.shape[0]+1):
                self.deriv_coeff.append(get_fornberg_coeffs(self.tmpx[i+self.n], self.tmpx[i:i+self.N-1]))
        else:
            for i in range(self.n):
                self.deriv_coeff.append(get_fornberg_coeffs(self.x[i], self.x[:self.N-1]))
            for i in range(self.n, self.x.shape[0] - self.n):
                self.deriv_coeff.append(get_fornberg_coeffs(self.x[i], self.x[i-self.n:i+self.n+1]))
            for i in range(self.x.shape[0] - self.n, self.x.shape[0]):
                self.deriv_coeff.append(get_fornberg_coeffs(self.x[i], self.x[self.x.shape[0] - self.N + 1:]))
        return None
    def compute_beta(self):
        self.neighbour_list = []
        if self.periodic:
            for i in range(len(self.deriv_coeff)-1):
                self.neighbour_list.append(range(i-self.n, i+self.n+2))
                deltax = np.array([self.dx[i]**l for l in range(self.m + 1)])
                a0 = self.alpha0_coeff*deltax[:, np.newaxis]
                a1 = self.alpha1_coeff*deltax[:, np.newaxis]
                btmp = [np.polynomial.polynomial.Polynomial(
                            list(np.sum(self.deriv_coeff[i][:self.m+1, 0, np.newaxis]*a0, axis = 0)))]
                for k in range(1, self.N-1):
                    btmp.append(np.polynomial.polynomial.Polynomial(np.sum(
                        self.deriv_coeff[i  ][:self.m+1, k  , np.newaxis]*a0
                      + self.deriv_coeff[i+1][:self.m+1, k-1, np.newaxis]*a1 , axis = 0)))
                btmp.append(np.polynomial.polynomial.Polynomial(list(
                        np.sum(self.deriv_coeff[i+1][:self.m+1, self.N-2, np.newaxis]*a1, axis = 0))))
                self.beta.append([[btmp[k].deriv(j)*self.dx[i]**(-j)
                                   for k in range(self.N)]
                                  for j in range(self.m+1)])
        else:
            for i in range(self.n):
                self.neighbour_list.append(range(2*self.n+2))
                deltax = np.array([self.dx[i]**l for l in range(self.m + 1)])
                a0 = self.alpha0_coeff*deltax[:, np.newaxis]
                a1 = self.alpha1_coeff*deltax[:, np.newaxis]
                btmp = []
                for k in range(self.N-1):
                    btmp.append(np.polynomial.polynomial.Polynomial(np.sum(
                        self.deriv_coeff[i  ][:self.m+1, k, np.newaxis]*a0
                      + self.deriv_coeff[i+1][:self.m+1, k, np.newaxis]*a1 , axis = 0)))
                btmp.append(np.polynomial.polynomial.Polynomial([0]))
                self.beta.append([[btmp[k].deriv(j)*self.dx[i]**(-j)
                                   for k in range(self.N)]
                                  for j in range(self.m+1)])
            for i in range(self.n, len(self.deriv_coeff)-self.n-1):
                self.neighbour_list.append(range(i-self.n, i+self.n+2))
                deltax = np.array([self.dx[i]**l for l in range(self.m + 1)])
                a0 = self.alpha0_coeff*deltax[:, np.newaxis]
                a1 = self.alpha1_coeff*deltax[:, np.newaxis]
                btmp = [np.polynomial.polynomial.Polynomial(
                            list(np.sum(self.deriv_coeff[i][:self.m+1, 0, np.newaxis]*a0, axis = 0)))]
                for k in range(1, self.N-1):
                    btmp.append(np.polynomial.polynomial.Polynomial(np.sum(
                        self.deriv_coeff[i  ][:self.m+1, k  , np.newaxis]*a0
                      + self.deriv_coeff[i+1][:self.m+1, k-1, np.newaxis]*a1 , axis = 0)))
                btmp.append(np.polynomial.polynomial.Polynomial(list(
                        np.sum(self.deriv_coeff[i+1][:self.m+1, self.N-2, np.newaxis]*a1, axis = 0))))
                self.beta.append([[btmp[k].deriv(j)*self.dx[i]**(-j)
                                   for k in range(self.N)]
                                  for j in range(self.m+1)])
            for i in range(len(self.deriv_coeff)-self.n-1, len(self.deriv_coeff)-1):
                self.neighbour_list.append(range(len(self.deriv_coeff) - 2*self.n - 1, len(self.deriv_coeff)))
                deltax = np.array([self.dx[i]**l for l in range(self.m + 1)])
                a0 = self.alpha0_coeff*deltax[:, np.newaxis]
                a1 = self.alpha1_coeff*deltax[:, np.newaxis]
                btmp = []
                for k in range(self.N-1):
                    btmp.append(np.polynomial.polynomial.Polynomial(np.sum(
                        self.deriv_coeff[i  ][:self.m+1, k, np.newaxis]*a0
                      + self.deriv_coeff[i+1][:self.m+1, k, np.newaxis]*a1 , axis = 0)))
                btmp.append(np.polynomial.polynomial.Polynomial([0]))
                self.beta.append([[btmp[k].deriv(j)*self.dx[i]**(-j)
                                   for k in range(self.N)]
                                  for j in range(self.m+1)])
        return None
    def compute_fast_beta(self):
        self.fast_beta = []
        for i in range(len(self.beta)):
            self.fast_beta.append([[sp.utilities.lambdify((self.xi),
                sp.horner(sp.Poly((self.beta[i][j][k].coef[::-1]), self.xi)), np)
                                    for k in range(len(self.beta[i][j]))]
                                   for j in range(self.m + 1)])
        return None
    def write_cfunction(
            self,
            cprefix = None,
            csuffix = None,
            data_type = 'float'):
        src_txt = 'int ' + cprefix + 'beta' + csuffix + '('
        if not self.periodic:
            src_txt += 'int cell, '  # which cell are we in?
        src_txt += (
                'int diff, ' +       # which derivative should we use?
                data_type + ' t, ' +
                data_type + ' *bval)' +     # array where to place the values of the beta polynomials
                '\n{\n')
        # sanity check
        src_txt += 'assert(diff >= 0 && diff <= {0});\n'.format(self.m)
        def beta_cformulas(node):
            tmp_txt = (
                    'switch (diff)\n{\n')
            for diff in range(self.m+1):
                tmp_txt += 'case {0}:\n'.format(diff)
                for i in range(-self.n, self.n + 2):
                    tmp_txt += 'bval[{0}] = '.format(i+self.n)
                    end_paranthesis = ''
                    for k in range(self.beta[node][diff][i+self.n].coef.shape[0] - 1):
                        tmp_txt += '({0}) + t*('.format(self.beta[node][diff][i+self.n].coef[k])
                        end_paranthesis += ')'
                    tmp_txt += '{0}'.format(self.beta[node][diff][i+self.n].coef[-1])
                    tmp_txt += end_paranthesis + ';\n'
                tmp_txt += 'break;\n'
            tmp_txt += '\n}\n'      # end diff switch
            return tmp_txt
        if self.periodic:
            src_txt += beta_cformulas(0)
        else:
            src_txt += (
                    'switch (cell)\n{\n')
            for cell in range(len(self.beta)):
                src_txt += 'case {0}:\n'.format(cell)
                src_txt += beta_cformulas(cell)
                src_txt += 'break;\n'
            src_txt += ('\n}\n')    # end cell switch
        # end and return 0
        src_txt += 'return EXIT_SUCCESS;\n}\n'
        src_txt += 'int ' + cprefix + 'indices' + csuffix + '('
        src_txt += (
                'int cell, ' +      # which cell are we in?
                'int *index)' +    # array where to place the values of the beta polynomials
                '\n{\n')
        if self.periodic:
            for i in range(self.n*2 + 2):
                src_txt += 'index[{0}] = {1};\n'.format(i, i-self.n)
        else:
            src_txt += (
                    'switch (cell)\n{\n')
            for cell in range(self.n):
                src_txt += 'case {0}:\n'.format(cell)
                for i in range(len(self.neighbour_list[cell])):
                    src_txt += 'index[{0}] = {1};\n'.format(i, self.neighbour_list[cell][i] - cell)
                if len(self.neighbour_list[cell]) < self.n*2+2:
                    src_txt += 'index[{0}] = {1};\n'.format(self.n*2+1, self.neighbour_list[cell][-1] - cell)
                src_txt += 'break;\n'
            for cell in range(len(self.beta)-self.n, len(self.beta)):
                src_txt += 'case {0}:\n'.format(cell)
                for i in range(len(self.neighbour_list[cell])):
                    src_txt += 'index[{0}] = {1};\n'.format(i, self.neighbour_list[cell][i] - cell)
                if len(self.neighbour_list[cell]) < self.n*2+2:
                    src_txt += 'index[{0}] = {1};\n'.format(self.n*2+1, self.neighbour_list[cell][-1] - cell)
                src_txt += 'break;\n'
            cell = self.n+1
            src_txt += 'default:\n'.format(cell)
            for i in range(len(self.neighbour_list[cell])):
                src_txt += 'index[{0}] = {1};\n'.format(i, self.neighbour_list[cell][i] - cell)
            if len(self.neighbour_list[cell]) < self.n*2+2:
                src_txt += 'index[{0}] = {1};\n'.format(self.n*2+1, self.neighbour_list[cell][-1] - cell)
            src_txt += 'break;\n'
            src_txt += ('\n}\n')    # end cell switch
        # end and return 0
        src_txt += 'return EXIT_SUCCESS;\n}\n'
        return src_txt

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
    tst0.compute_fast_beta()
    xval = []
    for i in range(x.shape[0]-1):
        xtmp = [x[i] + k*.1*(x[i+1] - x[i])
                         for k in range(10)]
        xval += xtmp
    xval = np.array(xval)
    if plt:
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_axes([.1, .1, .8, .8])
        ax.set_title('Weight functions for {0} neighbours and {1} continuous derivatives'.format(n, m))
        for i in range(n+1):
            y = np.zeros(x.shape, x.dtype)
            y[i] = 1
            tst0.put_yvals(y)
            f = np.array([tst0(xvar) for xvar in xval])
            ax.plot(xval, f)
        fig.savefig('test.pdf', format = 'pdf')
    else:
        print('didn\'t find matplotlib, so I\'m just gonna print out the weight functions.')
        print('here are the points where I\'m computing them.')
        print(xval)
        print('and here are the weight functions.')
        for i in range(n+1):
            y = np.zeros(x.shape, x.dtype)
            y[i] = 1
            tst0.put_yvals(y)
            f = np.array([tst0(xvar) for xvar in xval])
            print(f)
    return None

def main0():
    plot_generic_weight_functions(n = 4, m = 2)
    return None

if __name__ == '__main__':
    main0()

