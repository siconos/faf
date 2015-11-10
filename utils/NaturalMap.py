#!/usr/bin/env python

# failure!

# Joli & Feng p 326  missing ||xt|| -> 0 ?

# cf nm1.py + nm2.py


from sympy import var, Matrix, init_printing, Symbol, sqrt, Piecewise, piecewise_fold, And
from sympy.core.numbers import  NaN
from localcodegen import localccode

init_printing()

import JordanAlgebra

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
rn = Symbol('rn', positive=True, real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', positive=True, real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)


u = Matrix([un, ut1, ut2])
r = Matrix([rn, rt1, rt2])





def piecewise_matrix(p):

    p_f = piecewise_fold(p)

    if p_f.is_Piecewise:
        nbr = p_f.args[0][0].shape[0]
        nbc = p_f.args[0][0].shape[1]

        return Matrix(nbr, nbc, lambda i, j:
                      Piecewise(*[(e[i,j],c) for e, c in p_f.args]))

    else:
        return p


def cone_proj(mu, v):
    vn = v[0]
    vt = v[1:, 0]
    vt1 = vt[0]
    vt2 = vt[1]

    vt_norm = vt.norm()
    zero3 = Matrix([0., 0., 0.])

    return Piecewise((zero3, mu * vt_norm <= - vn),
                     (v, vt_norm <= mu * vn),
                     (v - ((vt_norm - mu * vn)/(1 + mu*mu)) * (Matrix([-mu]).col_join(vt/vt_norm)), vt_norm > mu * vn))

# x_star
xs = Matrix([
    un + mu * sqrt(ut1**2 + ut2**2),
    ut1, 
    ut2])

# r_star
rs = r - xs

rsn = rs[0]
rst = rs[1:, 0]
rst_norm = rst.norm()


xn =  un + mu * sqrt(ut1**2 + ut2**2)
xt1 = mu * ut1
xt2 = mu * ut2

yn = mu * rn
yt1 = rt1
yt2 = rt2

x = Matrix([xn, xt1, xt2])
y = Matrix([yn, yt1, yt2])


Fnat = y - JordanAlgebra.projection(y - x)
#Fnat = piecewise_matrix(Piecewise(
#    (xs + ((rst_norm-mu*rsn)/(1+mu*mu))*(Matrix([-mu]).col_join(rst/rst_norm)),
#     And(mu*rst_norm >= -rsn, rst_norm>=mu*rsn)),
#    (xs, And(mu*rst_norm>= -rsn, rst_norm<mu*rsn))))

theta_Fnat = .5*Fnat.norm()**2

#print Fnat

A = Fnat.jacobian(u)
B = Fnat.jacobian(r)

#print A_
#print B_

Rnow = Fnat.row_join(A).row_join(B)


ac_name = "fc3d_NaturalMap"

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

    print localccode(Fnat, assign_to='result', array_format='Fortran')

    end_fun()


if options.ccodeAB:
    def_fun(ac_name + "ABGenerated")

    print localccode(A.row_join(B), assign_to='result', array_format='Fortran')

    end_fun()

if options.merit:

    def_fun(ac_name + "FMeritGenerated")

    print localccode(theta_Fnat, assign_to='result', array_format='Fortran')

    end_fun()

    def_fun(ac_name + "GradFMeritGenerated")

    print localccode(grad_theta_Fnat, assign_to='result', array_format='Fortran')

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

