#!/usr/bin/env python

# print the Alart Curnier function (Standard form or Moreau &
# Jean) and its derivatives in latex or c code
#

#
# Usage ./fac.py [--JeanMoreau] [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#


from sympy import *
from sympy.core.numbers import  NaN
from localcodegen import localccode
from locallatex import print_latex_by_conditions

# not in python 2.7 -> see argparse
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--JeanMoreau", action="store_true")
parser.add_option("--latex", action="store_true")
parser.add_option("--ccode", action="store_true")
parser.add_option("--ccodefac", action="store_true")
parser.add_option("--ccodeAB", action="store_true")
parser.add_option("--wrapper", action="store_true")
(options, args) = parser.parse_args()



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


var('R_N R_T1 R_T2 U_N U_T1 U_T2 rho_N rho_T1 rho_T2 D_N D_T1 D_T2')

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
if (options.JeanMoreau):
    Radius = mu*_sup0(rn)
else:
    Radius = mu*max0D0 # max(0,D0)

# phi
phi2 =  rn - max0D0

phi31 = rt1 - px.subs(x,D1).subs(y,D2).subs(D,Radius)

phi32 = rt2 - py.subs(x,D1).subs(y,D2).subs(D,Radius)

u = Matrix([[un],[ut1],[ut2]])

r = Matrix([[rn],[rt1],[rt2]])

FAC = Matrix([[phi2],
              [phi31],
              [phi32]])

A = Matrix([[diff(phi2,un), diff(phi2,ut1), diff(phi2, ut2)],
            [diff(phi31,un), diff(phi31,ut1), diff(phi31, ut2)],
            [diff(phi32,un), diff(phi32,ut1), diff(phi32, ut2)]])

B = Matrix([[diff(phi2,rn), diff(phi2,rt1), diff(phi2, rt2)],
            [diff(phi31,rn), diff(phi31,rt1), diff(phi31, rt2)],
            [diff(phi32,rn), diff(phi32,rt1), diff(phi32, rt2)]])


#
# Toward Latex
#


phi_2 = Symbol(r'\phi_2')
phi_31 = Symbol(r'\phi_{3,T1}')
phi_32 = Symbol(r'\phi_{3,T2}')


def map_matrix(fun,M):
    m=M.shape[0]
    n=M.shape[1]
    return Matrix(m,n, lambda i,j: fun(M[i,j]))


def with_latex_symbols(expr):
    return expr.subs(D0,D_N).subs(D1,D_T1).subs(D2,D_T2).subs(rn,R_N).subs(rt1,R_T1).subs(rt2,R_T2).subs(un,U_N).subs(ut1,U_T1).subs(ut2,U_T2).subs(rhon,rho_N).subs(rhot1,rho_T1).subs(rhot2,rho_T2)


def print_latex_value(name,value):
    if (hasattr(name,'args')):
        lname = latex(name)
    else:
        lname = name
    print latex(r'\subsubsection{$' + lname + r'$}')
    print latex(r'\begin{align*}')
    print latex(with_latex_symbols(value))
    print latex(r'\end{align*}')
    print r'\paragraph{}'


textU = with_latex_symbols(u)
textR = with_latex_symbols(r)

textA = Matrix([[diff(phi_2(U_N),U_N), diff(phi_2(U_T1),U_T1), diff(phi_2(U_T2),U_T2)],
                [diff(phi_31(U_N),U_N), diff(phi_31(U_T1),U_T1), diff(phi_31(U_T2),U_T2)],
                [diff(phi_32(U_N),U_N), diff(phi_32(U_T1),U_T1), diff(phi_32(U_T2),U_T2)]])

textB = Matrix([[diff(phi_2(R_N),R_N), diff(phi_2(R_T1),R_T1), diff(phi_2(R_T2),R_T2)],
                [diff(phi_31(R_N),R_N), diff(phi_31(R_T1),R_T1), diff(phi_31(R_T2),R_T2)],
                [diff(phi_32(R_N),R_N), diff(phi_32(R_T1),R_T1), diff(phi_32(R_T2),R_T2)]])


def resultFAC(i,j):
    if (j==0):
        if(i==0):
            return r"\phi_2"
        if(i==1):
            return r"\phi_{3,T1}"
        if(i==2):
            return r"\phi_{3,T2}"
    if (j>=1 and j<=3):
            return r"A[{0},{1}]={2}".format(i,j-1,latex(textA[i,j-1]))
    if (j>3):
            return r"B[{0},{1}]={2}".format(i,j-4,latex(textB[i,j-4]))

def dump_latex():

    print (r"""% this file has been generated
\documentclass[10pt]{article}
\usepackage{amssymb,amsmath}
\begin{document}""")
    if options.JeanMoreau:
        print(r"""
\title{Alart \& Curnier function (Christensen \& Pang variant)}""")
    else:
        print(r"""
\title{Alart \& Curnier function}""")
    print(r"""
\date{\today}
\maketitle
""")

    print latex(r'\section{Symbols}')

    print r'\begin{itemize}'
    print r'\item ' + latex(R_N, mode='inline') + r': \text{normal reaction}'
    print r'\item ' + latex(R_T1, mode='inline') + r': \text{first tangential reaction}'
    print r'\item ' + latex(R_T2, mode='inline') + r': \text{second tangential reaction}'

    print r'\item ' + latex(U_N, mode='inline') + r': \text{normal velocity}'
    print r'\item ' + latex(U_T1, mode='inline') + r': \text{first tangential velocity}'
    print r'\item ' + latex(U_T2, mode='inline') + r': \text{second tangential velocity}'

    print r'\item ' + latex(rho_N, mode='inline') + r': \text{normal component of} $\rho$'
    print r'\item ' + latex(rho_T1, mode='inline') + r': \text{first tangential component of} $\rho$'
    print r'\item ' + latex(rho_T2, mode='inline') + r': \text{second tangential component of} $\rho$'
    print r'\end{itemize}'

    print r'\paragraph{basic relations}'

    print r'\begin{align*}'
    print latex(D_N) + '&=' + latex(R_N-U_N*rho_N) +r'\\'
    print latex(D_T1) + '&=' + latex(R_T1-U_T1*rho_T1) +r'\\'
    print latex(D_T2) + '&=' + latex(R_T2-U_T2*rho_T2) +r'\\'
    print r'\end{align*}'



    print(r'\section{$\phi_2,\phi_{3,T1},\phi_{3,T2}$,A,B by values}')


    print latex(r'\subsection{$\phi$}')
    print_latex_value(phi_2,phi2)
    print_latex_value(phi_31,FAC[1])
    print_latex_value(phi_32,FAC[2])

    print latex(r'\subsection{A}')

    for i in range(0,3):
        for j in range(0,3):
            print_latex_value("A[{0},{1}]={2}".format(i,j,latex(textA[i,j])),A[i,j])

    print latex('\subsection{B}')

    for i in range(0,3):
        for j in range(0,3):
            print_latex_value("B[{0},{1}]={2}".format(i,j,latex(textB[i,j])),B[i,j])

    R = FAC.row_join(A).row_join(B)

    print(r'\section{$\phi_2,\phi_{3,T1},\phi_{3,T2}$,A,B by conditions}')

    print_latex_by_conditions(Matrix(with_latex_symbols(R)),resultFAC)

    print(r"\end{document}")


if options.latex:
    dump_latex()

if options.JeanMoreau:
    ac_name = "frictionContact3D_AlartCurnierJeanMoreau"
else:
    ac_name = "frictionContact3D_AlartCurnier"

if options.ccode:
    def_fun(ac_name + "FABGenerated")
    Rnow = FAC.row_join(A).row_join(B)

    print localccode(Rnow, assign_to='result', array_format='Fortran')
    end_fun()

if options.ccodefac:
    def_fun(ac_name + "FGenerated")
    print localccode(FAC, assign_to='result', array_format='Fortran')
    end_fun()


if options.ccodeAB:
    def_fun(ac_name + "ABGenerated")
    print localccode(A.row_join(B), assign_to='result', array_format='Fortran')
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


