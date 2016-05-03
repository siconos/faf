from sympy import Symbol, Function, Matrix, sqrt

mu = Symbol('mu', negative=False, real=True, finite=True)
rn = Symbol('rn', negative=False, real=True)
rt1 = Symbol('rt1', real=True)
rt2 = Symbol('rt2', real=True)

un = Symbol('un', real=True)
ut1 = Symbol('ut1', real=True)
ut2 = Symbol('ut2', real=True)

#class N(Function):
#      nargs = (1, 2)

def N(x, y):
    return sqrt(x*x+y*y)
      
xn = un + mu * N(ut1, ut2)
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
