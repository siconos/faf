#!/usr/bin/env sage

#
# print the Fischer Burmeister function and its derivatives in latex or c code
#

#
# Usage ./fb.py [--latex] [--ccode] [--ccodefac] [--ccodeAB]
#

import codegen.JordanAlgebra
from sympy import *
from sympy.core.numbers import NaN
from codegen.funcodegen import funcodegen
from codegen.maple import denoms, limzero, Maple, set_maple, mlimit
from locallatex import print_latex_by_conditions
import pickle
import sys
import numpy as np

#from codegen.maple import limzero, Maple, set_maple
#maple = Maple(server='bizet.inria.fr')
#set_maple(maple)

#maple = Maple(server='bizet.inria.fr')
#set_maple(maple)
#t = Symbol('t', positive=True, real=True)
#maple_assumes = 'mu >= 0 and mu < 1 and rn >=0 and t > 0 and epsilon > 0'
#maple.assume(maple_assumes)

init_printing()
from codegen.friction import mu, rn, rt1, rt2, un, ut1, ut2, xn, xt1, xt2, yn, yt1, yt2, u, r, x, y, ut, epsilon

Fnat = y - codegen.JordanAlgebra.projection(y - x)
A_ = Fnat.jacobian(u)
B_ = Fnat.jacobian(r)

from sympy.core.numbers import NaN, Infinity, ComplexInfinity

def isfinite(x):
    return not(isinstance(x, NaN) or isinstance(x, Infinity) or isinstance(x, ComplexInfinity))

def load(v):

    with open(v) as r_file:
        return pickle.load(r_file)

A = load('Fnat_A')
B = load('Fnat_B')
Fnat = load('Fnat')
#Rnow = load('Fnat_Rnow')
# theta_phi_fb = load('theta_phi_fb')
# grad_theta_phi_fb = load('grad_theta_phi_fb')


def isnan(x):
    return isinstance(x, NaN)

#for i in range(3):
#    for j in range(3):
#        assert not isnan(A[i, j].subs('ut1', 0).subs('ut2', 0))
#        assert not isnan(A[i, j].subs('rt1', 0).subs('rt2', 0))
#        assert not isnan(A[i, j].subs('ut1', 0).subs(
#            'ut2', 0).subs('rt1', 0).subs('rt2', 0))
#        assert not isnan(A[i, j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))
#        assert not isnan(B[i, j].subs('ut1', 0).subs('ut2', 0))
#        assert not isnan(B[i, j].subs('rt1', 0).subs('rt2', 0))
#        assert not isnan(B[i, j].subs('ut1', 0).subs(
#            'ut2', 0).subs('rt1', 0).subs('rt2', 0))
#        assert not isnan(B[i, j].subs('ut1', 'rt1/mu').subs('ut2', 'rt2/mu'))

intervals = {mu: Interval(0, 1), rn: Interval(0, 1e6)}

deriv = {'A': A, 'B': B}

# vv = sys.argv[1]
# ii = int(sys.argv[2])
# jj = int(sys.argv[3])

function_name = sys.argv[1]
target_name = sys.argv[2]

# print funcodegen('natural_map{0}{1}{2}'.format(vv, ii, jj), deriv[vv][ii, jj], intervals=intervals, array_format='Fortran', epsilon_inf=np.finfo(float).eps, epsilon_power=3, assertions=True, main_check=True, contracts=False)
# print funcodegen('natural_mapA', A, intervals=intervals,
# array_format='Fortran', epsilon_inf=np.finfo(float).eps,
# epsilon_power=3, assertions=True, main_check=True, contracts=False)

EPSILON = np.finfo(float).eps

expression = {
    'fc3d_{0}FABGenerated'.format(function_name): Fnat.row_join(A).row_join(B),
    'fc3d_{0}FGenerated'.format(function_name):   Fnat,
    'fc3d_{0}ABGenerated'.format(function_name):  A.row_join(B)}

C_function_name = 'fc3d_{0}{1}Generated'.format(
    function_name, target_name)

var('rhon rhot1 rhot2')

print funcodegen(C_function_name, expression[C_function_name],
                variables=[rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1,
                           rhot2],
                intervals=intervals,
                array_format='Fortran',
                epsilon_inf=epsilon,
                epsilon_power=1,
                with_files_generation=True,
                user_exprs=[(2*mu*mu-2*mu+1,[]), (epsilon*(mu+1),['>= epsilon'])],
                assertions=True, contracts=False, main_check=True)

#for i in range(3):
#    for j in range(3):
#        print funcodegen(
#            'tA{0}{1}'.format(i, j), A[i, j],
#            variables=[
#                rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1,
#                rhot2],
#            intervals=intervals,
#            array_format='Fortran',
#            epsilon_inf=epsilon,
#            epsilon_power=1,
#            do_cse=True,
#            with_files_generation=True,
#            user_exprs=[(2*mu*mu-2*mu+1,[]), (epsilon*(mu+1),['>= epsilon'])],
#            assertions=True, contracts=False, main_check=True)
