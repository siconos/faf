#!/usr/bin/env python

# jacobian of sqrt(x*x+y*y)

#
# Usage ./test.py [--ccode]
#
import sys
sys.path.append('..')

from sympy import Symbol, Matrix, Piecewise, sqrt, limit
from funcodegen import funcodegen
import numpy as np


x = Symbol('x', real=True)
y = Symbol('y', real=True)
t = Symbol('t', real=True)

v = Matrix([x, y])

f = Matrix([sqrt(x * x + y * y)])

J_ = f.jacobian(v)

EPSILON = np.finfo(float).eps
# EPSILON=0.


def lim0(expr):
    return limit(expr.subs(x, t).subs(y, t), t, 0)

J = Matrix(J_.shape[0], J_.shape[1],
           lambda i, j: Piecewise(
               (J_[i, j], v.norm() > EPSILON),
               (lim0(J_[i, j]), v.norm() <= EPSILON)))


print funcodegen('norm2d_jacobian', J, array_format='Fortran',
                 assertions=True, contracts=False, epsilon_inf=EPSILON,
                 epsilon_power=2, main_check=True, do_cse=True,
                 with_files_generation=False)
