from sympy import Matrix, Piecewise, Max, piecewise_fold, Symbol
from sympy import sqrt as std_sqrt
import numpy as np
import friction

epsilon = Symbol('epsilon', real=True, positive=True)

def piecewise_matrix(p):

    p_f = piecewise_fold(p)

    if p_f.is_Piecewise:
        nbr = p_f.args[0][0].shape[0]
        nbc = p_f.args[0][0].shape[1]

        return Matrix(nbr, nbc, lambda i, j:
                      Piecewise(*[(e[i,j],c) for e, c in p_f.args]))

    else:
        return p

def prod(x, y):
    """
    Jordan product.
    """

    return (x.transpose() * y).col_join(
        y[0] * x[1:, 0] + x[0] * y[1:, 0])

def spectral_decomposition(x, epsilon=epsilon, norm=friction.N):
    """
    Jordan spectral decomposition.
    """
    xn = x[0]
    xt = x[1:, 0]
    xt_norm = norm(xt[0], xt[1])

    omega = Matrix([0, 1])

    lambda1 = xn - xt_norm
    lambda2 = xn + xt_norm

    u1 = (1./2.)*Matrix([
        1]).col_join(
            piecewise_matrix(Piecewise((-xt/xt_norm, xt_norm**2 > epsilon),
                                       (- omega, xt_norm**2 <= epsilon))))

    u2 = (1./2.)*Matrix([
        1]).col_join(
            piecewise_matrix(Piecewise((xt/xt_norm, xt_norm**2 > epsilon),
                                       (omega, xt_norm**2 <= epsilon))))

    return (lambda1, lambda2, u1, u2)

def sqrt(x, epsilon=epsilon, norm=friction.N):
    """
    Jordan Algebra sqrt
    """
    lambda1, lambda2, u_1, u_2 = spectral_decomposition(x, epsilon, norm)
    return std_sqrt(lambda1) * u_1 + std_sqrt(lambda2) * u_2

def projection(x, epsilon=epsilon, norm=friction.N):

    lambda1, lambda2, u_1, u_2 = spectral_decomposition(x, epsilon, norm)

    return Max(0, lambda1) * u_1 + Max(0, lambda2) * u_2

