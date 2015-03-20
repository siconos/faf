#!/usr/bin/env sage

#
# print the Fischer Burmeister function and its derivatives in FAC A B Rnow pickle files
#

# Usage: ./fb1.py

from sympy import *
from sympy.core.numbers import  NaN
import pickle

import sage
from sage.all import maple


init_printing()

ZERO = Symbol('ZERO', real=True)

mu = Symbol('mu', positive=True, real=True)
rn = Symbol('rn', positive=True, real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', positive=True, real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)

fn = Function('fn')
ft1 = Function('ft1')
ft2 = Function('ft2')

#xn = fn(mu, un, ut1, ut2)
#xt1 = ft1(mu, un, ut1, ut2)
#xt2 = ft2(mu, un, ut1, ut2)

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

assert x.shape == (3, 1)
assert y.shape == (3, 1)

# with spectral decomposition

xt = Matrix(x[1:])
yt = Matrix(y[1:])

assert xt.shape == (2, 1)
assert yt.shape == (2, 1)


random1 = var('random1')
random2 = var('random2')

xnxt_p_ynyt = xn*xt + yn*yt

rand_v = Matrix([random1, random2])
rand_v_norm = sqrt(random1**2 + random2**2);
n_rand_v = rand_v / rand_v_norm

xnxt_p_ynyt_norm = sqrt(xnxt_p_ynyt[0]**2 + xnxt_p_ynyt[1]**2)
x_norm = sqrt(xn**2 + xt1**2 + xt2**2)
y_norm = sqrt(yn**2 + yt1**2 + yt2**2)

_xu_1 = Piecewise(((xn * xt + yn * yt)[0, 0] / xnxt_p_ynyt_norm,
                    xnxt_p_ynyt_norm > 0),
                  (n_rand_v[0, 0], xnxt_p_ynyt_norm <= 0))

_xu_2 = Piecewise(((xn * xt + yn * yt)[1, 0] / xnxt_p_ynyt_norm,
                   xnxt_p_ynyt_norm > 0),
                  (n_rand_v[1, 0], xnxt_p_ynyt_norm <= 0))



u_1 = 0.5 * Matrix([1,
                    -_xu_1,
                    -_xu_2])

u_2 = 0.5 * Matrix([1,
                    _xu_1,
                    _xu_2])


lambda_1 = x_norm**2 + y_norm**2 - 2 * xnxt_p_ynyt_norm
lambda_2 = x_norm**2 + y_norm**2 + 2 * xnxt_p_ynyt_norm

phi_fb = x + y - (sqrt(lambda_1) * u_1 + sqrt(lambda_2) * u_2)

FAC = phi_fb

# ut_norm == 0
# xnxt_p_ynyt_norm == 0
# x_norm == 0
# y_norm == 0
# lambda_1 == 0
# lambda_2 == 0

A_ = phi_fb.jacobian(u)
B_ = phi_fb.jacobian(r)

t = Symbol('t', positive=True, real=True)

from sympy.functions.elementary.piecewise import ExprCondPair


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

    maple.set('Digits', 36)

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

maple.assume('un>=0')
maple.assume('rn>=0')
maple.assume('t>0')

#print mlimit(A_[2,2].subs(ut1, t).subs(ut2, t).subs(rt1, t).subs(rt2, t), t, 0)

#print cse(mlimit(A_[2,2].subs(ut1, t).subs(ut2, t).subs(rt1, t).subs(rt2, t), t, 0))

#exit(0)

def utzero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t)
        return mlimit(xe, t, 0)
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def uzero(e):


    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(un, t)
        return mlimit(xe, t, 0)
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def rzero(e):


    if hasattr(e, 'subs'):
        xe = e.subs(rt1, t).subs(rt2, t).subs(rn, t)
        return mlimit(xe, t, 0)
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def xnxt_p_ynyt_zero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(rt1, t).subs(rt2, t)
        #return limit(xe._sage_(), t=0, dir='+')._sympy_()
        return mlimit(xe, t, 0)

    else:
        return e

def lambda_1_zero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(un, t).subs(rt1, t).subs(rt2, t).subs(rn, t)

        return mlimit(xe, t, 0)


def lambda_2_zero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(un, t).subs(rt1, t).subs(rt2, t).subs(rn, t)

        return mlimit(xe, t, 0)

A = Matrix(3, 3, lambda i, j:
           Piecewise(
                     (uzero(A_[i, j]),  x_norm <= ZERO),
                     (xnxt_p_ynyt_zero(A_[i, j]), xnxt_p_ynyt_norm<=ZERO),
                     (lambda_1_zero(A_[i, j]), Abs(lambda_1)<=ZERO),
                     (lambda_2_zero(A_[i, j]), Abs(lambda_2)<=ZERO),
                     (utzero(A_[i, j]), ut.norm()<=ZERO),
                     (A_[i, j], And(ut.norm()>ZERO, x_norm>ZERO, xnxt_p_ynyt_norm>ZERO, Abs(lambda_1)>ZERO, Abs(lambda_2)>ZERO))))


B = Matrix(3, 3, lambda i, j:
           Piecewise((rzero(B_[i, j]),  y_norm<=ZERO),
                     (xnxt_p_ynyt_zero(B_[i, j]), xnxt_p_ynyt_norm<=ZERO),
                     (lambda_1_zero(B_[i, j]), Abs(lambda_1)<=ZERO),
                     (lambda_2_zero(B_[i, j]), Abs(lambda_2)<=ZERO),
                     (B_[i, j], And(y_norm>ZERO, xnxt_p_ynyt_norm>ZERO, Abs(lambda_1)>ZERO, Abs(lambda_2)>ZERO))))


with open('A', 'w') as A_file:
    pickle.dump(A, A_file)

with open('B', 'w') as B_file:
    pickle.dump(B, B_file)

with open('FAC', 'w') as FAC_file:
    pickle.dump(FAC, FAC_file)

Rnow = FAC.row_join(A).row_join(B)

with open('Rnow', 'w') as Rnow_file:
    pickle.dump(Rnow, Rnow_file)
