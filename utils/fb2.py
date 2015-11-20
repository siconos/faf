#!/usr/bin/env python

#
# print the Fischer Burmeister function and its derivatives in latex or c code
#

#
# Usage ./fb.py [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#


from sympy import *
from sympy.core.numbers import  NaN
from localcodegen import localccode
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

def load(v):

    with open(v) as r_file:
        return pickle.load(r_file)

A = load('FB_A')
B = load('FB_B')
FB = load('FB')
Rnow = load('FB_AB')
theta_phi_fb = load('theta_phi_fb')
grad_theta_phi_fb = load('grad_theta_phi_fb')

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


ac_name = "fc3d_FischerBurmeister"

if options.ccode:
    def_fun(ac_name + "FABGenerated")

    print localccode(Rnow, assign_to='result', array_format='Fortran')

    end_fun()

if options.ccodefac:
    def_fun(ac_name + "FGenerated")

    print localccode(FB, assign_to='result', array_format='Fortran')

    end_fun()


if options.ccodeAB:
    def_fun(ac_name + "ABGenerated")

    print localccode(A.row_join(B), assign_to='result', array_format='Fortran')

    end_fun()

if options.merit:

    def_fun(ac_name + "FMeritGenerated")

    print localccode(theta_phi_fb, assign_to='result', array_format='Fortran')

    end_fun()

    def_fun(ac_name + "GradFMeritGenerated")

    print localccode(grad_theta_phi_fb, assign_to='result', array_format='Fortran')

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



