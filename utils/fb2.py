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
(options, args) = parser.parse_args()

def load(v):

    with open(v) as r_file:
        return pickle.load(r_file)

A = load('A')
B = load('B')
FAC = load('FAC')
Rnow = load('Rnow')

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

#    nrows, ncols = Rnow.shape

#    for i in range(0, nrows):
#        for j in range(0, ncols):
#            print '{'
#            print localccode(Rnow[i, j], assign_to='result[{0}]'.format(i+j*nrows), level=1)
#            print '}'

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



