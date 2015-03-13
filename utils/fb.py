#!/usr/bin/env sage

#
# print the Fischer Burmeister function and its derivatives in latex or c code
#

#
# Usage ./fb.py [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#


from sympy import *
from sympy.core.numbers import  NaN
from localcodegen import dump_ccode, dump_var
from locallatex import print_latex_by_conditions
import pickle

init_printing()

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


mu = Symbol('mu', positive=True, real=True)
rn = Symbol('rn', real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', real=True)
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



def mlimit(expr, var, lim, dir=None):
    import sage 

    if dir is None:
        dir = sage.all.maple.right
    return sage.all.maple('limit({0}, {1}={2}, {3})'.format(expr._sage_(),
                                                            var, 
                                                            lim, dir))._sage_()._sympy_()

def piecewise_f(expr, f, i, j):

    print i, j, f

    p_expr = piecewise_fold(expr)

    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        args = [(piecewise_f(e, f, i, j), c) for e, c in p_expr.args]

        return Piecewise(*args)
    else:

        return f(p_expr)



def utzero(e):
    import sage

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t)
        return mlimit(xe, t, 0)
#        return sage.all.limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def uzero(e):
    import sage

    if hasattr(e, 'subs'):
        xe = e.subs(ut1, t).subs(ut2, t).subs(un, t)
        return mlimit(xe, t, 0)
#        return sage.all.limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def rzero(e):
    import sage

    if hasattr(e, 'subs'):
        xe = e.subs(rt1, t).subs(rt2, t).subs(rn, t)
        return mlimit(xe, t, 0)
#        return sage.all.limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e

def xnxt_p_ynyt_zero(e):
    import sage

    if hasattr(e, 'subs'):
        xe = e.subs(un, t).subs(rn, t).subs(ut1, t).subs(ut2, t).subs(rt1, t).subs(rt2, t)
        
        xxe = sage.all.maple.simplify(xe)

        print xxe

        return mlimit(xxe, t, 0)
#        return sage.all.limit(xe._sage_(), t=0, dir='+')._sympy_()
    else:
        return e


print piecewise_f(A_[1, 0], uzero, 1, 0)

for i in range(0, 3):
    for j in range(0, 3):
        print i, j
        piecewise_f(A_[i, j], utzero, i, j)
        piecewise_f(A_[i, j], uzero, i, j)
        piecewise_f(A_[i, j], xnxt_p_ynyt_zero, i, j)



A = Matrix(3, 3, lambda i, j: 
           Piecewise(
                     (piecewise_f(A_[i, j], uzero, i, j),  x_norm<=0),
                     (piecewise_f(A_[i, j], xnxt_p_ynyt_zero, i, j), And(xnxt_p_ynyt_norm<=0, x_norm>0)),
                     (piecewise_f(A_[i, j], xnxt_p_ynyt_zero, i, j), And(lambda_1<=0, xnxt_p_ynyt_norm>0, x_norm>0)),
                     (piecewise_f(A_[i, j], xnxt_p_ynyt_zero, i, j), And(lambda_2<=0, lambda_1>0, xnxt_p_ynyt_norm>0, x_norm>0)),
                     (piecewise_f(A_[i, j], utzero, i, j), And(ut.norm()<=0,lambda_2>0, lambda_1>0, xnxt_p_ynyt_norm>0, x_norm>0)),
                     (A_[i, j], And(ut.norm()>0, x_norm>0, xnxt_p_ynyt_norm>0, lambda_1>0, lambda_2>0))))


B = Matrix(3, 3, lambda i, j: 
           Piecewise((piecewise_f(B_[i, j], rzero, i, j),  y_norm<=0),
                     (piecewise_f(B_[i, j], xnxt_p_ynyt_zero, i, j), xnxt_p_ynyt_norm<=0),
                     (piecewise_f(B_[i, j], xnxt_p_ynyt_zero, i, j), lambda_1<=0),
                     (piecewise_f(B_[i, j], xnxt_p_ynyt_zero, i, j), lambda_2<=0),
                     (B_[i, j], And(y_norm>0, xnxt_p_ynyt_norm>0, lambda_1>0, lambda_2>0))))


with open('A', 'w') as A_file:
    pickle.dump(A, A_file)

with open('B', 'w') as B_file:
    pickle.dump(B, B_file)

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

def change_var():
    pass

if options.ccode:
    def_fun(ac_name + "FABGenerated")
    Rnow = FAC.row_join(A).row_join(B)

    nrows, ncols = Rnow.shape

    for i in range(0, nrows):
        for j in range(0, ncols):
            print '{'
            dump_ccode(Rnow[i, j], array_format='fortran', offset=(i, j*nrows))
            print '}'
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



