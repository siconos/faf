import sage
from sage.all import Maple
from funcodegen import Memoize
from sympy import piecewise_fold, Piecewise, Matrix, sympify, simplify,\
    GreaterThan, StrictGreaterThan, LessThan, StrictLessThan

maple = None
assumes = None


def set_maple(umaple):
    global maple
    maple = umaple


def _to_maple(expr):

    p_expr = piecewise_fold(expr)
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        args = []
        for e, c in p_expr.args:
            args += [c, to_maple(e)]
        return maple('piecewise({0})'.format(', '.join([str(e) for e in args])))
    else:
        return maple(str(p_expr))

to_maple = Memoize(_to_maple)


def _to_sympy(mexpr):
    if type(mexpr) == sage.interfaces.maple.MapleElement:
        try:
            return mexpr._sage_()._sympy_()
        except:
            mstr = str(maple.simplify(mexpr))
            mstr = mstr.replace('D[1](Max)(0,','Heaviside(')
            mstr = mstr.replace('D[2](Max)(0,','DiracDelta(')
            if mstr[0:9] == 'piecewise':
                args = list(level1_args(mstr))
                l = [(_to_sympy(maple.simplify(e)), _to_sympy(maple.simplify(c)))
                     for e, c in zip(args[1::2], args[::2])]
                return Piecewise(*l)
            else:
                result = None
                for comp, tcomp in [(
                    '>=', GreaterThan),
                    ('<=', LessThan),
                    ('>', StrictGreaterThan),
                        ('<', StrictLessThan)]:
                    if comp in mstr:
                        args = mstr.split(comp)
                        result = tcomp(
                            _to_sympy(maple.simplify(args[0])), _to_sympy(maple.simplify(args[1])))
                        break
                if result is not None:
                    return result
                else:
                    print mstr
                    return sympify(mstr)

    else:
        if hasattr(mexpr, 'sage'):
            return mexpr._sage_()._sympy_()
        else:
            if hasattr(mexpr, '_sympy_'):
                return mexpr._sympy_()
            else:
                return mexpr

to_sympy = Memoize(_to_sympy)


def _maple_simplify(mexpr):
    return maple.simplify(mexpr)

maple_simplify = Memoize(_maple_simplify)


def piecewise_simplify(expr):

    p_expr = piecewise_fold(expr)
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        return Piecewise(*[(e, c) for (e, c) in p_expr.args])
    else:
        return p_expr


def _mlimit(expr, var, lim, dir=None):

    if dir is None:
        dir = maple.right

    maple.set('Digits', 64)

    mexpr = to_maple(expr)

    # failure if the limit is piecewise!
    return _to_sympy(maple('limit({0}, {1}={2}, {3})'.format(mexpr,
                                                    var,
                                                    lim, dir)))


mlimit = Memoize(_mlimit)


def _x_simplify(expr):
    if hasattr(expr, 'is_Matrix') and expr.is_Matrix:
        return Matrix(expr.shape[0], expr.shape[1], lambda i, j: x_simplify(expr[i, j]))

    p_expr = piecewise_fold(expr)
    try:
        r0 = cse([p_expr])
        if is_iterable(r0) and is_iterable(r0[0]):
            len0 = len(r0[0])
        else:
            len0 = 0
    except:
        len0 = 0
    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:
        return Piecewise(*[(x_simplify(e), x_simplify(c)) for (e, c) in p_expr.args])
    elif hasattr(p_expr, 'is_Relational') and p_expr.is_Relational:
        return type(p_expr)(x_simplify(p_expr.args[0]), x_simplify(p_expr.args[1]))
    else:
        mexpr = maple.simplify(p_expr)
        r = to_sympy(mexpr)
        try:
            r1 = cse([r])
            if is_iterable(r1) and is_iterable(r1[0]):
                len1 = len(r1[0])
            else:
                len1 = 0
        except:
            len1 = 0
        if len1 < len0:
            return r
        else:
            return expr

x_simplify = Memoize(_x_simplify)


def fix(expr):

    if hasattr(expr, 'subs') and not hasattr(expr, '_maple_'):
        xexpr = expr
#        for subexpr in list(postorder_traversal(expr)):
#            xexpr = xexpr.subs(
#                subexpr, x_simplify(subexpr))

        return x_simplify(xexpr)
    else:
        return expr


def limzero(e, lv):
    if hasattr(e, 'subs'):
        xe = e
        limvs = []
        for v, vsub in lv:
            xe = xe.subs(v, vsub)
            if hasattr(vsub, 'is_Symbol') and vsub.is_Symbol:
                limvs.append(vsub)
        result = xe
        for limv in limvs:
            result = mlimit(result, limv, 0)
        return fix(result)
    else:
        return e

from sympy.core.numbers import Infinity, NegativeInfinity


def minimize(e, intervals={}):
    intervkeys = dict((str(k), k) for k in intervals.keys())
    minterv = []
    for sym in e.free_symbols:
        if str(sym) in intervkeys:
            bmin = intervals[intervkeys[str(sym)]].args[0]
            bmax = intervals[intervkeys[str(sym)]].args[1]
            minterv += ['{0}={1}..{2}'.format(sym, bmin, bmax)]
    if len(minterv) > 0:
        try:
            return to_sympy(maple(
                'minimize({0}, {1})'.format(e, ','.join(minterv))))
        except:
            return NegativeInfinity
    else:
        return to_sympy(maple.minimize(e))

# minimize = Memoize(_minimize)


def maximize(e, intervals={}):
    intervkeys = dict((str(k), k) for k in intervals.keys())
    minterv = []
    for sym in e.free_symbols:
        if str(sym) in intervkeys:
            bmin = intervals[intervkeys[str(sym)]].args[0]
            bmax = intervals[intervkeys[str(sym)]].args[1]
            minterv += ['{0}={1}..{2}'.format(sym, bmin, bmax)]
    if len(minterv) > 0:
        try:
            return to_sympy(maple(
                'maximize({0}, {1})'.format(e, ','.join(minterv))))
        except:
            return Infinity
    else:
        return to_sympy(maple.maximize(e))

# maximize = Memoize(_maximize)

from sympy import postorder_traversal, Pow


def denoms(expr):
    return set(
        to_sympy(
            e.args[0]) for e in filter(lambda e: type(e) == Pow and e.args[1] < 0,
                                       postorder_traversal(expr)))


def all_denoms_bounds(expr, intervals={}):
    return (((e, (minimize(e, intervals=intervals), maximize(e, intervals=intervals))) for e in denoms(expr)))

from sympy import And

def is_finite(expr):
    if hasattr(expr, 'is_finite'):
        return expr.is_finite
    else:
        if type(expr) == sage.rings.real_double.RealDoubleElement:
            return True

def finite_denoms_bounds(expr, intervals={}):
    return dict(filter(lambda b: is_finite(b[1][0]) and is_finite(b[1][1]),
                       all_denoms_bounds(expr, intervals=intervals)))
