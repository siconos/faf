#!/usr/bin/env sage

#
# print the Fischer Burmeister function and its derivatives in FAC A B Rnow pickle files
#

# Usage: ./fb1.py

from sympy import *
from sympy.core.numbers import  NaN
import pickle

import sage
from sage.all import maple, Maple
import numpy as np

init_printing()

maple = Maple(server="bizet.inria.fr")

#ZERO = Symbol('ZERO', real=True)
EPSILON = np.finfo(float).eps


mu = Symbol('mu', positive=True, real=True)
rn = Symbol('rn', real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)

xn =  un + mu * sqrt(ut1**2 + ut2**2)
xt1 = mu * ut1
xt2 = mu * ut2

yn = mu * rn
yt1 = rt1
yt2 = rt2

u = Matrix([un, ut1, ut2])
r = Matrix([rn, rt1, rt2])

x = Matrix([xn, xt1, xt2])
y = Matrix([yn, yt1, yt2])

ut = Matrix([ut1, ut2])
rt = Matrix([rt1, rt2])

assert x.shape == (3, 1)
assert y.shape == (3, 1)

# with spectral decomposition
import JordanAlgebra

Fnat = y - JordanAlgebra.projection(y - x)

FAC = Fnat

A_ = Fnat.jacobian(u)
B_ = Fnat.jacobian(r)

t = Symbol('t', positive=True, real=True)

class Memoize():
    def __init__(self, fun):
        self._fun = fun
        self._done = dict()

    def __call__(self, *args):
        if args in self._done:
            return self._done[args]
        else:
            try:
                r = self._fun(*args)
                self._done[args] = r
                return r
            except Exception as e:
                self._done[args] = e
                return e

def _maple_simplify(expr):

    return maple.simplify(expr)

maple_simplify = Memoize(_maple_simplify)

def _maple_piecewise(expr):

    p_expr = piecewise_fold(expr)
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        args = []
        for e, c in p_expr.args:
            args += [c, maple_piecewise(e)]
        return 'piecewise({0})'.format(', '.join([str(e) for e in args]))
    else:
        return p_expr

maple_piecewise = Memoize(_maple_piecewise)


def _mlimit(expr, var, lim, dir=None):

    if dir is None:
        dir = maple.right

    maple.set('Digits', 256)

    mexpr = maple_piecewise(expr)

    # failure if the limit is piecewise!
    return maple('limit({0}, {1}={2}, {3})'.format(mexpr,
                                                   var,
                                                   lim, dir))._sage_()._sympy_()


mlimit = Memoize(_mlimit)

maple.assume('mu', 'real')

maple.assume('un', 'real')
maple.assume('ut1', 'real')
maple.assume('ut2', 'real')

maple.assume('rn', 'real')
maple.assume('rt1', 'real')
maple.assume('rt2', 'real')

maple.assume('mu>0')

maple.assume('rn>=0')
maple.assume('t>0')



def msubs(e, syms):

    me = maple_piecewise(e)

    for sym1, sym2 in syms:

        me = maple('subs({0}={1}, {2})'.format(sym1, sym2, me))

    return me


# find a limit of e when ||ut|| -> 0
def utzero(e):

    if hasattr(e, 'subs'):
#        xe = msubs(e, [(ut1, t), (ut2, t)])
        xe = e.subs(ut1, t).subs(ut2, t)
        return mlimit(xe, t, 0)
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

# find a limit of e when (mu*ut1 - rt1)^2  + (mu*ut2 - rt2)^2 -> 0
def projzero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t*rt1/mu).subs(ut2, t*rt2/mu)
        return mlimit(xe, t, 0)


# ||u|| -> 0 ||r|| -> 0
def utrtzero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, 0).subs(rt1, 0).subs(rt2, t)
        return mlimit(xe, t, 0)


A = Matrix(3, 3, lambda i, j:
           Piecewise(
               (utrtzero(A_[i, j]), ut.norm() + rt.norm() <= EPSILON),
               (utzero(A_[i, j]), ut.norm() <= EPSILON),
               (A_[i, j], ut.norm() > EPSILON)))


B = Matrix(3, 3, lambda i, j:
           Piecewise(
               (utrtzero(B_[i, j]), ut.norm() + rt.norm() <= EPSILON), 
               (projzero(B_[i, j]), (mu*ut1 - rt1)**2 + (mu*ut2 - rt2)**2 <= EPSILON),
               (B_[i, j], (mu*ut1 - rt1)**2 + (mu*ut2 - rt2)**2 > EPSILON)))

def isnan(x):
    return isinstance(x, NaN)


print isnan(A_[2,2].subs(ut1, 0).subs(ut2, 0))
print isnan(B_[2,2].subs(ut1, rt1/mu).subs(ut2, rt2/mu))


# bad substitutions in sympy!!!
for i in range(3):
    for j in range(3):
        print i, j
        assert not isnan(A[i,j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(A[i,j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(A[i,j].subs('ut1', 0).subs('ut2', 0).subs('rt1', 0).subs('rt2', 0).subs('rn', 0).subs('un', 0))
        assert not isnan(A[i,j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))
        assert not isnan(B[i,j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(B[i,j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(B[i,j].subs('ut1', 0).subs('ut2', 0).subs('rt1', 0).subs('rt2', 0).subs('rn', 0).subs('un', 0))
        assert not isnan(B[i,j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))


with open('Fnat_A', 'w') as A_file:
    pickle.dump(A, A_file)

with open('Fnat_B', 'w') as B_file:
    pickle.dump(B, B_file)

with open('Fnat', 'w') as FAC_file:
    pickle.dump(Fnat, FAC_file)

Rnow = FAC.row_join(A).row_join(B)

with open('Fnat_Rnow', 'w') as Rnow_file:
    pickle.dump(Rnow, Rnow_file)

theta_Fnat = .5 * (Fnat[0]**2 + Fnat[1]**2 + Fnat[2]**2)

with open('theta_Fnat', 'w') as theta_Fnat_file:
    pickle.dump(theta_Fnat, theta_Fnat_file)

grad_theta_Fnat = Matrix([[theta_Fnat.diff(rn), theta_Fnat.diff(rt1), theta_Fnat.diff(rt2)]])

with open('grad_theta_Fnat', 'w') as grad_theta_Fnat_file:
    pickle.dump(grad_theta_Fnat, grad_theta_Fnat_file)

# ut_norm == 0
# xnxt_p_ynyt_norm == 0
# x_norm == 0
# y_norm == 0
# lambda_1 == 0
# lambda_2 == 0
