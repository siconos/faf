#!/usr/bin/env sage

#
# print the Natural Map function and its derivatives in Fnat A B Rnow pickle files
#

# Usage: ./nm1.py

from sympy import *
from sympy.core.numbers import NaN
from sympy.core.relational import StrictGreaterThan, GreaterThan, StrictLessThan, LessThan
import pickle

import sage
from sage.all import Maple
import numpy as np

init_printing()

maple = Maple(server="bizet.inria.fr")

# ZERO = Symbol('ZERO', real=True)
EPSILON = np.finfo(float).eps


mu = Symbol('mu', negative=False, real=True)
rn = Symbol('rn', real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)

xn = un + mu * sqrt(ut1 ** 2 + ut2 ** 2)
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

xmyt = Matrix([xt1 - yt1, xt2 - yt2])

assert x.shape == (3, 1)
assert y.shape == (3, 1)

# with spectral decomposition
import JordanAlgebra

Fnat = y - JordanAlgebra.projection(y - x)


FAC = Fnat

A_ = Fnat.jacobian(u)
B_ = Fnat.jacobian(r)


t = Symbol('t', positive=True, real=True)

maple_assumes = 'mu >= 0 and mu < 1 and rn >=0 and t > 0'

maple.assume(maple_assumes)


class Memoize():

    def __init__(self, fun):
        self._fun = fun
        self._done = dict()

    def __call__(self, *args):
        if args in self._done:
            return self._done[args]
        else:
            r = self._fun(*args)
            self._done[args] = r
            return r


def _to_maple(expr):

    p_expr = piecewise_fold(expr)
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        args = []
        for e, c in p_expr.args:
            args += [c, to_maple(e)]
        return maple('piecewise({0})'.format(', '.join([str(e) for e in args])))
    else:
        return maple(str(p_expr))

to_maple = Memoize(_to_maple)


def level1_args(istr):
    """
    args at first level
    """
    p = 0
    current_arg = ''
    state = 0
    while p < len(istr):
        current_char = istr[p]
        p += 1

        if current_char == '(':
            if state > 0:
                current_arg += current_char
            state += 1

        elif current_char == ')':
            if state > 1:
                current_arg += current_char
            state -= 1
            if state == 0:
                yield current_arg

        elif state == 1:
            if current_char == ',':
                return_arg = current_arg
                current_arg = ''
                yield return_arg
            else:
                current_arg += current_char

        else:
            if state > 0:
                current_arg += current_char


def _to_sympy(mexpr):
    if type(mexpr) == sage.interfaces.maple.MapleElement:
        try:
            return mexpr._sage_()._sympy_()
        except:
            mstr = str(maple.simplify(mexpr))
            if mstr[0:9] == 'piecewise':
                args = list(level1_args(mstr))
                l = [(_to_sympy(maple.simplify(e)), _to_sympy(maple.simplify(c)))
                     for e, c in zip(args[1::2], args[::2])]
                return Piecewise(*l)
            else:
                result = None
                for comp, tcomp in [(
                    '>=', GreaterThan),
                    ('<=', LessThan),
                    ('>', StrictGreaterThan),
                        ('<', StrictLessThan)]:
                    if comp in mstr:
                        args = mstr.split(comp)
                        result = tcomp(
                            _to_sympy(maple.simplify(args[0])), _to_sympy(maple.simplify(args[1])))
                        break
                if result is not None:
                    return result
                else:
                    print type(mstr), mstr
                    assert False

    else:
        print type(mexpr), mexpr
        assert False

to_sympy = Memoize(_to_sympy)

# P = Piecewise((ut1, ut1 * ut1 <= 0), (1, ut1 * ut1 > 0))

# print to_maple(P)
# print to_sympy(to_maple(P))
# exit(0)


def _maple_simplify(mexpr):
    return maple.simplify(mexpr)

maple_simplify = Memoize(_maple_simplify)


def piecewise_simplify(expr):

    p_expr = piecewise_fold(expr)
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        return Piecewise(*[(e, c) for (e, c) in p_expr.args])
    else:
        return p_expr


def _mlimit(expr, var, lim, dir=None):

    if dir is None:
        dir = maple.right

    maple.set('Digits', 256)

    mexpr = to_maple(expr)

    # failure if the limit is piecewise!
    return maple('limit({0}, {1}={2}, {3})'.format(mexpr,
                                                   var,
                                                   lim, dir))._sage_()._sympy_()


mlimit = Memoize(_mlimit)


def msimplify(expr):
    xexpr = expr
    if hasattr(expr, 'expand'):
        xexpr = expand(xexpr)
    if hasattr(expr, 'factor'):
        xexpr = factor(xexpr)
#    if hasattr(expr, 'collect'):
#        xexpr = collect(xexpr)
    if hasattr(expr, 'cancel'):
        xexpr = cancel(xexpr)
    return xexpr

from codegen.funcodegen import flatten_piecewise


def fix(expr):

    if hasattr(expr, 'subs') and not hasattr(expr, '_maple_'):
        xexpr = expr
        for subexpr in list(postorder_traversal(expr)):
            xexpr = xexpr.subs(
                subexpr, factor(subexpr))

        return xexpr
    else:
        return expr


def msubs(e, syms):

    me = to_maple(e)

    for sym1, sym2 in syms:

        me = maple('subs({0}={1}, {2})'.format(sym1, sym2, me))

    return me


# find a limit of e when ||ut|| -> 0
def utzero(e):

    if hasattr(e, 'subs'):
#        xe = msubs(e, [(ut1, t), (ut2, t)])
        xe = e.subs(ut1, t).subs(ut2, t)
        return fix(mlimit(xe, t, 0))
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

# find a limit of e when ||rt|| -> 0


def rtzero(e):

    if hasattr(e, 'subs'):
#        xe = msubs(e, [(ut1, t), (ut2, t)])
        xe = e.subs(rt1, 0).subs(rt2, 0)
        return fix(xe)
#        return mlimit(xe, t, 0)
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

# find a limit of e when (mu*ut1 - rt1)^2  + (mu*ut2 - rt2)^2 -> 0


def projzero(e):

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t * rt1 / mu).subs(ut2, t * rt2 / mu)
        return fix(mlimit(xe, t, 0))
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

# ||u|| -> 0 ||r|| -> 0


def utrtzero(e):
    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(rt1, 0).subs(rt2, 0)
        return fix(mlimit(xe, t, 0))
#        return limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e


for i in range(0, 3):
    for j in range(0, 3):
        print '======================='
        print i, j
        print 'utrtzero(A_[i,j])=', utrtzero(A_[i, j])
        print 'utzero(A_[i,j])=',  utzero(A_[i, j])
        print 'rtzero(A_[i,j])=', rtzero(A_[i, j])

        print 'utrtzero(B_[i,j])=', utrtzero(B_[i, j])
        print 'projzero(B_[i,j])=', projzero(B_[i, j])
        print 'utzero(B_[i,j])=',  utzero(B_[i, j])
        print 'rtzero(B_[i,j])=',  rtzero(B_[i, j])

A = Matrix(3, 3, lambda i, j:
           Piecewise(
               (utrtzero(A_[i, j]), And(
                   ut.norm() <= EPSILON, rt.norm() <= EPSILON)),
               (rtzero(A_[i, j]),   And(
                   ut.norm() > EPSILON,  rt.norm() <= EPSILON)),
               (utzero(A_[i, j]),   And(
                   ut.norm() <= EPSILON, rt.norm() > EPSILON)),
               (projzero(A_[i, j]), And(xmyt.norm() <= EPSILON,
                                        ut.norm() > EPSILON,
                                        rt.norm() > EPSILON)),
               (A_[i, j],           And(xmyt.norm() > EPSILON,
                                        ut.norm() > EPSILON,
                rt.norm() > EPSILON))))

B = Matrix(3, 3, lambda i, j:
           Piecewise(
               (utrtzero(B_[i, j]), And(
                   ut.norm() <= EPSILON, rt.norm() <= EPSILON)),
               (rtzero(B_[i, j]),   And(
                   ut.norm() > EPSILON,  rt.norm() <= EPSILON)),
               (utzero(B_[i, j]),   And(
                   ut.norm() <= EPSILON, rt.norm() > EPSILON)),
               (projzero(B_[i, j]), And(xmyt.norm() <= EPSILON,
                                        rt.norm() > EPSILON,
                ut.norm() > EPSILON)),
               (B_[i, j], And(
                xmyt.norm() > EPSILON, rt.norm() > EPSILON,
                ut.norm() > EPSILON))))


def isnan(x):
    return isinstance(x, NaN)


print isnan(A_[2, 2].subs(ut1, 0).subs(ut2, 0))
print isnan(B_[2, 2].subs(ut1, rt1 / mu).subs(ut2, rt2 / mu))


# bad substitutions in sympy!!!
for i in range(3):
    for j in range(3):
        print i, j
        assert not isnan(A[i, j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(A[i, j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(A[i, j].subs('ut1', 0).subs('ut2', 0).subs(
            'rt1', 0).subs('rt2', 0).subs('rn', 0).subs('un', 0))
        assert not isnan(A[i, j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))
        assert not isnan(B[i, j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(B[i, j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(B[i, j].subs('ut1', 0).subs('ut2', 0).subs(
            'rt1', 0).subs('rt2', 0).subs('rn', 0).subs('un', 0))
        assert not isnan(B[i, j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))


with open('Fnat_A', 'w') as A_file:
    pickle.dump(A, A_file)

with open('Fnat_B', 'w') as B_file:
    pickle.dump(B, B_file)

with open('Fnat', 'w') as FAC_file:
    pickle.dump(Fnat, FAC_file)

Rnow = FAC.row_join(A).row_join(B)

with open('Fnat_Rnow', 'w') as Rnow_file:
    pickle.dump(Rnow, Rnow_file)

theta_Fnat = .5 * (Fnat[0] ** 2 + Fnat[1] ** 2 + Fnat[2] ** 2)

with open('theta_Fnat', 'w') as theta_Fnat_file:
    pickle.dump(theta_Fnat, theta_Fnat_file)

grad_theta_Fnat = Matrix(
    [[theta_Fnat.diff(rn), theta_Fnat.diff(rt1), theta_Fnat.diff(rt2)]])

with open('grad_theta_Fnat', 'w') as grad_theta_Fnat_file:
    pickle.dump(grad_theta_Fnat, grad_theta_Fnat_file)


# intervals = {mu : Interval(0, 1), rn : Interval(0, 1e6)}

# from codegen.funcodegen import funcodegen

# print funcodegen('natural_mapA', A, intervals=intervals,
# array_format='Fortran', epsilon_inf=np.finfo(float).eps,
# epsilon_power=3, assertions=True, main_check=True, contracts=False)

# ut_norm == 0
# xnxt_p_ynyt_norm == 0
# x_norm == 0
# y_norm == 0
# lambda_1 == 0
# lambda_2 == 0
