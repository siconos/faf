#!/usr/bin/env python

#
# print the Fischer Burmeister function and its derivatives in latex or c code
#

#
# Usage ./fb.py [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#


from sympy import *
from localcodegen import dump_ccode
from locallatex import print_latex_by_conditions

# not in python 2.7 -> see argparse
from optparse import OptionParser
parser = OptionParser()
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

r = Matrix([rn, rt1, rt2]).transpose()

u = Matrix([un, ut1, ut2]).transpose()

rho = Matrix([rhon, rhot1, rhot2]).transpose()


uun = un + mu * sqrt(ut1*ut1+ut2*ut2)
uut1 = ut1
uut2 = ut2



yn = mu * rn
yt1 = rt1
yt2 = rt2

xn = uun
xt1 = mu * ut1
xt2 = mu * ut2

x = Matrix([[xn], [xt1], [xt2]])
y = Matrix([[yn], [yt1], [yt2]])

assert x.shape == (3,1)
assert y.shape == (3,1)

#1. direct computation 
x_o_y = Matrix((x.transpose() * y)[:] + (y[0,0] * Matrix(x[1:]) + x[0,0] * Matrix(y[1:]))[:])

assert x_o_y.shape == (3,1)

#print x_o_y

x_o_x = x_o_y.subs(y,x)

y_o_y = x_o_y.subs(x,y)

s = sqrt(.5*(x[0,0] + sqrt(x[0,0]**2 - (x[1,0]**2+x[2,0]**2))))

sqrt_x = Matrix([[Piecewise((0, s<=0),(s, s>0))],
                 [Piecewise((0, s<=0),(x[1,0]/ (2.* s), s>0))],
                 [Piecewise((0, s<=0),(x[2,0]/ (2.* s), s>0))]])

#print x
#print y
#print x_o_y



phi_fb1 = x + y + sqrt_x.subs(x, x_o_x + y_o_y)


#2. with spectral decomposition



xt = Matrix(x[1:])
yt = Matrix(y[1:])

assert xt.shape == (2,1)
assert yt.shape == (2,1)

Rand = Function('Rand')

xnxt_p_ynyt = xn * xt + yn * yt
xnxt_p_ynyt_norm = sqrt(xnxt_p_ynyt[0,0]**2 + xnxt_p_ynyt[1,0]**2)

rand_v = Matrix([Rand(1), Rand(2)])

n_rand_v = rand_v / rand_v.norm()

_xu_1 = Piecewise((xnxt_p_ynyt[0,0] / xnxt_p_ynyt_norm, xnxt_p_ynyt_norm>0),
                  (n_rand_v[0,0], xnxt_p_ynyt_norm <= 0))

_xu_2 = Piecewise((xnxt_p_ynyt[1,0] / xnxt_p_ynyt_norm, xnxt_p_ynyt_norm>0),
                  (n_rand_v[1,0], xnxt_p_ynyt_norm <= 0))


# Abs in x.norm() even with assumptions.

x_norm = sqrt(x[0,0]**2 + x[1,0]**2 + x[2,0]**2)
y_norm = sqrt(y[0,0]**2 + y[1,0]**2 + y[2,0]**2)


lambda_1 = x_norm ** 2 + y_norm ** 2 - 2 * xnxt_p_ynyt_norm
lambda_2 = x_norm ** 2 + y_norm ** 2 + 2 * xnxt_p_ynyt_norm


u_1 = 0.5 * Matrix([1, - _xu_1, - _xu_2])
u_2 = 0.5 * Matrix([1, _xu_1, _xu_2])

phi_fb2 = x + y - (sqrt(Max(0,lambda_1)) * u_1 + sqrt(Max(0,lambda_2)) * u_2)

#print phi_fb1
#print phi_fb2

phi_fb = phi_fb2

phi2 =  phi_fb[0,0]

phi31 = phi_fb[1,0]

phi32 = phi_fb[2,0]

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
    return expr


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
    print(r"""
\title{Fischer - Burmeister function}""")
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

ac_name = "frictionContact3D_FischerBurmeister"

if options.ccode:
    def_fun(ac_name + "FABGenerated")
    Rnow = FAC.row_join(A).row_join(B)
    dump_ccode(Rnow, array_format='fortran')
    end_fun()

if options.ccodefac:
    def_fun(ac_name + "FGenerated")
    dump_ccode(FAC, array_format='fortran')
    end_fun()


if options.ccodeAB:
    def_fun(ac_name + "ABGenerated")
    dump_ccode(A.row_join(B), array_format='fortran')
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
