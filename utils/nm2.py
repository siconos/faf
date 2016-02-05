#!/usr/bin/env python

#
# print the Fischer Burmeister function and its derivatives in latex or c code
#

#
# Usage ./fb.py [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#

import JordanAlgebra
from sympy import *
from sympy.core.numbers import  NaN
from codegen import localccode
from locallatex import print_latex_by_conditions
import pickle
import sys

init_printing()

# not in python 2.7 -> see argparse
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--latex", action="store_true")
parser.add_option("--ccode", action="store_true")
parser.add_option("--ccodefac", action="store_true")
parser.add_option("--ccodeAB", action="store_true")
parser.add_option("--wrapper", action="store_true")
parser.add_option("--merit", action="store_true")
(options, args) = parser.parse_args()

mu = Symbol('mu', positive=True, real=True)
rn = Symbol('rn', nonnegative=True, real=True)
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

Fnat = y - JordanAlgebra.projection(y - x)
A_ = Fnat.jacobian(u)
B_ = Fnat.jacobian(r)

def load(v):

    with open(v) as r_file:
        return pickle.load(r_file)

A = load('Fnat_A')
B = load('Fnat_B')
Fnat = load('Fnat')
Rnow = load('Fnat_Rnow')
#theta_phi_fb = load('theta_phi_fb')
#grad_theta_phi_fb = load('grad_theta_phi_fb')

def isnan(x):
    return isinstance(x, NaN)

for i in range(3):
    for j in range(3):
        assert not isnan(A[i,j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(A[i,j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(A[i,j].subs('ut1', 0).subs('ut2', 0).subs('rt1', 0).subs('rt2', 0))
        assert not isnan(A[i,j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))
        assert not isnan(B[i,j].subs('ut1', 0).subs('ut2', 0))
        assert not isnan(B[i,j].subs('rt1', 0).subs('rt2', 0))
        assert not isnan(B[i,j].subs('ut1', 0).subs('ut2', 0).subs('rt1', 0).subs('rt2', 0))
        assert not isnan(B[i,j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))

def def_fun(name):
    print("void " + name + r"""(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{""")

def end_fun():
    print(r"}")


ac_name = "fc3d_NaturalMap"


if options.ccode:
    def_fun(ac_name + "FABGenerated")

#    nrows, ncols = Rnow.shape

#    for i in range(0, nrows):
#        for j in range(0, ncols):
#            print '{'
#            print localccode(Rnow[i, j], assign_to='result[{0}]'.format(i+j*nrows), level=1)
#            print '}'

    print localccode(Rnow, assign_to='result', array_format='Fortran', epsilon_inf=1e-10, assertions=True)

    end_fun()

if options.ccodefac:
    def_fun(ac_name + "FGenerated")

    print localccode(Fnat, assign_to='result', array_format='Fortran')

    end_fun()


if options.ccodeAB:
    def_fun(ac_name + "ABGenerated")

    print localccode(A.row_join(B), assign_to='result', array_format='Fortran')

    end_fun()

if options.merit:

    def_fun(ac_name + "FMeritGenerated")

#    print localccode(theta_phi_fb, assign_to='result', array_format='Fortran')

    end_fun()

    def_fun(ac_name + "GradFMeritGenerated")

#    print localccode(grad_theta_phi_fb, assign_to='result', array_format='Fortran')

    end_fun()


if options.wrapper:

    print(
        """
void {0}FunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{{
  double result[21];

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if (f && A && B)
  {{

    {0}FABGenerated(
      *reaction0, *reaction1, *reaction2,
      *velocity0, *velocity1, *velocity2,
      mu,
      *rho0, *rho1, *rho2,
      result);
    cpy3(result, f);
    cpy3x3(result + 3, A);
    cpy3x3(result + 12, B);
  }}

  else
  {{
    if (f)
    {{
      {0}FGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }}

    if (A && B)
    {{
      {0}ABGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3x3(result, A);
      cpy3x3(result + 9, B);
    }}
  }}
}}
        """.format(ac_name))

