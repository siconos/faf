#!/usr/bin/env python

# print the Alart Curnier function (Standard form or Moreau &
# Jean) and its derivatives in latex or c code
#

#
# Usage ./fac.py [--JeanMoreau] [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#

import sys
sys.path.append('..')

from sympy import *
from funcodegen import funcodegen
import numpy as np
EPSILON=np.finfo(float).eps

mu = Symbol('mu', positive=True, real=True)
rn = Symbol('rn', real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)

rhon = Symbol('rhon', real=True)
rhot1 = Symbol('rhot1', real=True)
rhot2 = Symbol('rhot2', real=True)

D = Symbol('D', real=True)

D0 = rn-rhon*un
D1 = rt1-rhot1*ut1
D2 = rt2-rhot2*ut2

x = Wild('x')
y = Wild('y')

# max(0,x)
#_sup0 = Lambda(x, Piecewise((0, x<=0),(x, x>0)))
_sup0 = Lambda(x, Max(0, x))


# not the same derivative in 0
#_sup0 = Lambda(x, Rational(1,2)*(x+abs(x)))


# disk projection
px = Piecewise((x,Matrix([x,y]).norm()<=D),(D*x/Matrix([x,y]).norm(),Matrix([x,y]).norm()>D))
py = Piecewise((y,Matrix([x,y]).norm()<=D),(D*y/Matrix([x,y]).norm(),Matrix([x,y]).norm()>D))
max0D0 = _sup0(D0)

#
#  The function definition
#

# cf Dynamic in the presence of unilateral contacts and dry friction :
# a numerical approach. M Jean, J.J. Moreau. Unilateral Problems in
# Structural Analysis - 2 International Centre for Mechanical Sciences
# Volume 304, 1987, pp 151-196
if sys.argv[1:] == 1:
    Radius = mu*_sup0(rn) # JeanMoreau
else:
    Radius = mu*max0D0 # max(0,D0)

# phi
phi2 =  rn - max0D0

phi31 = rt1 - px.subs(x,D1).subs(y,D2).subs(D,Radius)

phi32 = rt2 - py.subs(x,D1).subs(y,D2).subs(D,Radius)

u = Matrix([[un],[ut1],[ut2]])

r = Matrix([[rn],[rt1],[rt2]])

FAC_ = Matrix([[phi2],
              [phi31],
              [phi32]])

FAC = Matrix(FAC_.shape[0], FAC_.shape[1],
             lambda i, j: Piecewise((FAC_[i, j], Radius > EPSILON),
                                    (FAC_[i, j].subs(Radius, 0), Radius <= EPSILON)))


A = FAC.jacobian(u)
B = FAC.jacobian(r)

if sys.argv[1] == 1:
    print funcodegen('JeanMoreau', FAC.row_join(A).row_join(B),
                     array_format='Fortran',
                     epsilon_inf=EPSILON,
                     assertions=True, main_check=True)
else:
    print funcodegen('AlartCurnier', FAC.row_join(A).row_join(B),
                     array_format='Fortran',
                     epsilon_inf=EPSILON,
                     assertions=True, main_check=True)
