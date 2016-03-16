from sympy import Matrix, Piecewise, Max, piecewise_fold


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

import numpy as np

def spectral_decomposition(x, epsilon=np.finfo(float).eps):
    """
    Jordan spectral decomposition.
    """
    xn = x[0]
    xt = x[1:, 0]
    xt_norm = xt.norm()

    omega = Matrix([0, 1])

    lambda1 = xn - xt_norm
    lambda2 = xn + xt_norm

    u1 = (1./2.)*Matrix([
        1]).col_join(
            piecewise_matrix(Piecewise((-xt/xt_norm, xt_norm > epsilon),
                                       (- omega, xt_norm <= epsilon))))

    u2 = (1./2.)*Matrix([
        1]).col_join(
            piecewise_matrix(Piecewise((xt/xt_norm, xt_norm > epsilon),
                                       (omega, xt_norm <= epsilon))))

    return (lambda1, lambda2, u1, u2)

def projection(x):

    lambda1, lambda2, u_1, u_2 = spectral_decomposition(x)

    return Max(0, lambda1) * u_1 + Max(0, lambda2) * u_2

