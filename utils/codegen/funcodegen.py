#!/usr/bin/env python

from sympy import *
from sympy.core.mul import _keep_coeff
from sympy.core.logic import Logic
from sympy.utilities.iterables import postorder_traversal
from sympy.core.compatibility import string_types
from sympy.printing.codeprinter import CodePrinter, Assignment
from sympy.printing.ccode import ccode, CCodePrinter
from sympy.printing.precedence import precedence
from sympy.printing.str import StrPrinter


class Memoize():

    def __init__(self, fun):
        self._fun = fun
        self._done = dict()

    def __call__(self, *args):
        if args in self._done:
            return self._done[args]
        else:
            r = self._fun(*args)
            self._done[args] = r
            return r

# dummy sqrt for cse and to avoid sympy transforms on Pow(_x, Rational(1, 2))


class fsqrt(Function):

    def _eval_is_negative(self):
        return False

    def _eval_is_positive(self):
        return self._args[0].is_positive

# dummy mul function to avoid x*x*x*x* ... *x -> Pow(x, n)


class fmul(Function):

    def _eval_is_negative(self):
        return Mul(self._args[0], self._args[1]).is_negative


def print_float(x):
    return StrPrinter({'full_prec': False}).doprint(Float(x, 3))


def append(d, k, v):
    if k in d:
        if not v in d[k]:
            d[k].append(v)
    else:
        d[k] = list()
        d[k].append(v)


def is_not_zero(a, epsilon):
    if epsilon == 0.:
        return '{0} != 0'.format(a)
    else:
        return '{0} < -{1} || {0} > {1}'.format(a, epsilon)


def is_positive_or_zero(a, epsilon):
    if epsilon == 0.:
        return '{0} >= 0'.format(a)
    else:
        return '{0} >= -{1}'.format(a, epsilon)


def is_strict_positive(a, epsilon):
    if epsilon == 0.:
        return '{0} > 0'.format(a)
    else:
        return '{0} > -{1}'.format(a, epsilon)

    
def assert_is_not_zero(a, epsilon):

    return '/*@ assert {0}; */'.format(is_not_zero(a, epsilon))


def assert_is_positive_or_zero(a, epsilon):

    return '/*@ assert {0}; */'.format(is_positive_or_zero(a, epsilon))


def asserts(l):
    if len(l) > 0:
        return '\n'.join(['/*@ assert {0}; */'.format(a) for a in l])
    else:
        return ''


class Is():

    """
    check type of an object x:
    Is(x).Matrix <=> hasattr(x, 'is_Matrix') and x.is_Matrix
    """

    def __init__(self, x):
        self._x = x

    def __getattr__(self, name):
        if hasattr(self._x, 'is_{0}'.format(name)):
            return getattr(self._x, 'is_{0}'.format(name))
        else:
            return False


class SomeVars():

    def __init__(self, name='x'):
        self._name = name
        self._counter = 0

    def __iter__(self):
        while(True):
            self._counter += 1
            yield Symbol('{0}{1}'.format(self._name, self._counter), real=True)


def flatten_piecewise(expr, upper_conds=None, conds=None, lconds=None):

    if hasattr(expr, 'is_Matrix') and expr.is_Matrix:
        return Matrix(expr.shape[0], expr.shape[1],
                      lambda i, j: flatten_piecewise(expr[i, j]))

    if upper_conds is None:
        upper_conds = []

    if conds is None:
        conds = dict()

    if lconds is None:
        lconds = []

    p_expr = piecewise_fold(expr)

    if hasattr(p_expr, 'is_Piecewise') and p_expr.is_Piecewise:

        for e, c in p_expr.args:
            # O(n)
            if c not in upper_conds:
                current_conds = upper_conds + [c]

                t_cc = tuple(current_conds)

                if t_cc in conds:
                    print p_expr
                    print t_cc, e
                    assert(False)

                else:
                    if not piecewise_fold(e).is_Piecewise:
                        conds[t_cc] = e
                        lconds.append(t_cc)
                    else:
                        flatten_piecewise(e, current_conds, conds, lconds)

    if len(upper_conds) == 0:

        if len(conds) > 0:
            return Piecewise(*[(conds[k], And(*(list(reversed(k))))) for k in filter(lambda k_: satisfiable(k_), lconds)])
        else:
            return p_expr

    else:
        return None

class LocalCCodePrinter(CCodePrinter):

    def __init__(self, settings={}, tab='    ', level=0, array_format='C',
                 epsilon_nan=0, epsilon_inf=0, epsilon_power=1, assertions=False,
                 contracts=False,
                 postcheck_hooks=False, do_cse=True):
        CCodePrinter.__init__(self, settings)
        self._some_vars = SomeVars()
        self._value_type = 'double'
        self._relational_type = 'int'
        self._assignments_values = dict()
        self._decls = dict()
        self._varid = dict()
        self._declared = dict()
        self._affcts = dict()
        self._cond_affcts = dict()
        self._tab = tab
        self._level = level
        self._do_cse = do_cse
        self._array_format = array_format
        self._epsilon_nan = epsilon_nan
        self._epsilon_inf = epsilon_inf
        self._epsilon_power = epsilon_power
        self._assertions = assertions
        self._contracts = contracts
        self._postcheck_hooks = postcheck_hooks
        self._current_condition = None
        self._sign = dict()

    def implied_positiveness(self, a, expr):
        if self._epsilon_inf == 0.:
            return '{0} > 0 ==> {1} > 0'.format(a, expr.args[0])
        else:
            return '{0} > {1} ==> {2} > {3}'.format(
                a,
                self._print(self._epsilon_inf),
                self._print(expr.args[0]),
                self._print(expr.subs(expr.args[0], self._epsilon_inf)))

    def recurs_mul(self, expr, nn):

        def raw_mul(_expr, _nn):
            def irm(x, n, ac):

                if n == 1:
                    return ac
                else:
                    return irm(x, n - 1, fmul(x, ac, evaluate=False))

            if nn > 0:
                return irm(_expr, _nn, _expr)
            elif nn == 0.:
                return 1.
            else:
                return irm(1. / _expr, -_nn, 1. / _expr)

        result = raw_mul(expr, nn)

        return result

    def _preconditions(self, var, expr):

        if not self._assertions:
            return ''

        if type(expr) == fsqrt:
            return self._preconditions(var, Pow(expr.args[0], Rational(1, 2)))

        def parenthesize(e):
            return self.parenthesize(e, precedence(e))

        sa0 = ['']

        sa = ['', '/*@']

        requires = []

        ensures = []
        assigns = ['assigns {0};'.format(var)]

        if expr.is_negative is not None and not expr.is_negative:
            ensures.append(is_positive_or_zero(var, 0.))

        if Is(expr).Relational or Is(expr).Boolean:
            ensures.append('{0} <==> ({1})'.format(var, self._print(expr)))
        else:
            assigns.append('assigns {0};'.format(var))
            ensures.append(
                'same_sign({0}, {1})'.format(var, self._print(expr)))
            ensures.append(
                r'\is_finite(({0}) {1})'.format(self._value_type, var))

        if Is(expr).Pow:

            sestr = self._print(expr.base)

            requires.append(r'\is_finite(({0}) ({1}))'.format(self._value_type,
                                                              sestr))

            if abs(float(expr.exp) - int(expr.exp)) > 0:
                requires.append(is_positive_or_zero(sestr, self._epsilon_nan))

            if expr.exp < 0:
                requires.append(is_not_zero(sestr, pow(self._epsilon_inf,
                                                       self._epsilon_power)))

            if expr.exp.is_even:
                ensures.append(is_positive_or_zero(var, self._epsilon_nan))

            if expr.exp.is_Rational:

                if expr.exp == Rational(1, 2):
                    ensures.append(is_positive_or_zero(var, self._epsilon_nan))

                if expr.exp == -Rational(1, 2):
                    ensures.append(is_strict_positive(var, self._epsilon_inf))

        if len(requires) > 0:
            sa0 += ['/*@ assert {0}; */'.format(r)
                    for r in map(parenthesize, requires)]
            sa += ['requires {0};'.format(r)
                   for r in map(parenthesize, requires)]

        sa += assigns

        pc = self._postconditions(var, expr)
        if len(pc) > 0:
            ensures += pc

        if len(ensures) > 0:
            sa += ['ensures {0};'.format(e) for e in map(parenthesize,
                                                         ensures)]
        if self._contracts:
            if len(requires) > 0:
                return sa0 + sa + ['*/']
            else:
                return sa + ['*/']

        else:
            if len(requires) > 0:
                return sa0
            else:
                return []

    def _postconditions(self, lhs, rhs, on_lhs=True):

        if not self._assertions:
            return ''

        sa = []

        if on_lhs:
            sestr = '({0})'.format(self._print(lhs))
        else:
            sestr = '({0})'.format(self._print(rhs))

        if Is(rhs).Relational or Is(rhs).Boolean:
            sa.append('{0} <==> ({1})'.format(
                self._print(lhs), self._print(rhs)))

        if not Is(rhs).Relational and not isinstance(rhs,
                                                     (Logic, And, Or, Not)):

            sa.append(r'\is_finite(({0}) {1})'.format(
                self._value_type, sestr))

            if not Is(rhs).Piecewise:

                if Is(rhs).positive:
                    sa.append(is_positive_or_zero(sestr, 0.))
                    sa.append(is_not_zero(sestr, 0.))
#                    sa.append(self.implied_strict_positive(sestr, rhs))
                    
                else:
                    if Is(rhs).negative is not None and not Is(rhs).negative:
                        sa.append(is_positive_or_zero(sestr, 0.))
                        if type(rhs) == fsqrt:
                            sa.append(self.implied_positiveness(sestr, Pow(rhs.args[0], 2)))
                    else:
                        if Is(rhs).negative is not None and Is(rhs).negative:
                            sa.append(is_positive_or_zero(self.parenthesize(rhs, precedence(
                                - rhs)), 0.))
                        else:
                            if Is(rhs).nonzero:
                                sa.append(is_not_zero(sestr, 0.))
        return sa

    def _print_MatrixElement(self, expr):
        if self._array_format == 'C':
            return "{0}[{1}]".format(expr.parent, expr.j +
                                     expr.i * expr.parent.shape[1])

        elif self._array_format == 'Fortran':
            return "{0}[{1}]".format(expr.parent, expr.i +
                                     expr.j * expr.parent.shape[0])

    def _traverse_matrix_indices(self, mat):

        rows, cols = mat.shape
        if self._array_format == 'C':
            return ((i, j) for i in range(rows) for j in range(cols))
        elif self._array_format == 'Fortran':
            return ((i, j) for j in range(cols) for i in range(rows))

    def _fix_integer_power(self, expr):
        subs = dict()

        for subexpr in list(postorder_traversal(expr)):
            if Is(subexpr).Pow:
                if Is(subexpr.args[1]).Integer and subexpr.args[1] > 0:
                    expr = expr.subs(
                        subexpr, Mul(*([subexpr.args[0]] * subexpr.args[1])))

        return expr

    def _needed_symbols(self, expr):

        l = set()
        symbols = set()

        for subexpr in postorder_traversal(expr):
            if hasattr(subexpr, 'free_symbols'):
                symbols = symbols.union(subexpr.free_symbols)

        for symb in symbols:
            if symb in self._decls:
                l.add(symb)
                l = l.union(self._needed_symbols(self._decls[symb]))
        return l

    def _cse(self, expr):

        l = set()
        y = Wild('y')
        p = Wild('p')

        exprs = set()

        expr_n = expr.\
            replace(sqrt(2),
                    Float(sqrt(2).evalf(n=128), 128)). \
            replace(Pow(y, -Rational(3, 2)),
                    lambda y: 1. / Mul(y, fsqrt(y),
                                       evaluate=False)).\
            replace(Pow(y, Rational(3, 2)),
                    lambda y: Mul(y, fsqrt(y),
                                  evaluate=False)).\
            replace(Pow(y, Rational(5, 2)),
                    lambda y: Mul(y, y, fsqrt(y),
                                  evaluate=False)).\
            replace(Pow(y, -Rational(5, 2)),
                    lambda y: 1. / Mul(y, y, fsqrt(y),
                                       evaluate=False)).\
            replace(Pow(y, Rational(1, 2)),
                    lambda y: fsqrt(y)).\
            replace(Pow(y, -Rational(1, 2)),
                    lambda y: 1. / fsqrt(y))
#                 replace(lambda expr: expr.is_Pow and expr.args[1].is_Integer and expr.args[1]>2,
#                         lambda z: self.recurs_mul(z.args[0], z.args[1])).\
#                 replace(lambda expr: expr.is_Pow and expr.args[1].is_Integer and expr.args[1]<-2,
# lambda z: 1./(self.recurs_mul(z.args[0], -z.args[1])))

        for subexpr in postorder_traversal(expr_n):

# AC & JM failure
#            if Is(subexpr).Piecewise:
#                for (e, c) in subexpr.args:
#                    exprs.add(e)
#                    exprs.add(c)

            if Is(subexpr).Function and not subexpr.is_Piecewise:
                exprs.add(subexpr)
            elif Is(subexpr).Pow:
                exprs.add(subexpr)
                exprs.add(subexpr.args[0])

                    
# AC & JM failure
            if Is(subexpr).Mul or Is(subexpr).Add:
                for e in subexpr.args:
                    exprs.add(e)

        return cse([expr_n] + list(exprs), self._some_vars)

    def _print_declarations(self, lexpr):

        l = set()
        decls = []
        for expr in lexpr:
            l = l.union(self._needed_symbols(expr))

        for var, _ in sorted([(s, self._varid[s]) for s in l], key=lambda t: t[1]):
            if var in self._decls and not var in self._declared:
                if self._decls[var].is_Relational or isinstance(self._decls[var], (Logic, And, Or, Not)):
                    decls.append(
                        '{0} {1} = 0;'.format(self._relational_type, var))
                else:
                    decls.append('{0} {1} = 0.;'.format(self._value_type, var))

                self._declared[var] = True

        if len(decls) > 0:

            return '\n'.join(decls) + '\n'

        else:
            return ''

    def _print_affectations(self, lexpr, condition=None):

        l = set()
        affcts = []
        for expr in lexpr:
            l = l.union(self._needed_symbols(expr))

        def tag(var):
            if condition is not None:
# affcts.append('// cond tag({0})={1}'.format(var,var in self._affcts or
# (var in self._cond_affcts and condition in self._cond_affcts[var])))
                return var in self._affcts or (var in self._cond_affcts and condition in self._cond_affcts[var])
            else:
# affcts.append('// tag({0})={1}'.format(var,var in self._affcts))
                if self._current_condition is None:
                    return var in self._affcts
                else:
                    return var in self._affcts or (var in self._cond_affcts and self._current_condition in self._cond_affcts[var])

        for var, _ in sorted([(s, self._varid[s]) for s in l], key=lambda t: t[1]):
            if var in self._decls and not tag(var):
                for sa in self._preconditions(var, self._decls[var]):
                    affcts.append(sa)
                postcheck = None
                if self._postcheck_hooks:
                    if isinstance(self._decls[var], Add):
                        postcheck = 'POST_CHECK_ADD'
                    elif isinstance(self._decls[var], Mul):
                        postcheck = 'POST_CHECK_MUL'
                    elif isinstance(self._decls[var], Pow):
                        postcheck = 'POST_CHECK_POW'
                    else:
                        postcheck = 'POST_CHECK'
#                affcts.append('/*@ assigns {0}; */'.format(var))
                if postcheck is not None:
                    affcts.append('{0} = {1}; {2}({0});'.format(
                        var, super(LocalCCodePrinter, self)._print(self._decls[var]), postcheck))
                else:
                    affcts.append(
                        '{0} = {1};'.format(var, super(LocalCCodePrinter, self)._print(self._decls[var])))

                affcts.append(
                    asserts(self._postconditions(var, self._decls[var])))

                if condition is not None:
                    if var in self._cond_affcts:
                        assert type(self._cond_affcts[var]) == list
                        self._cond_affcts[var].append(condition)
# affcts.append('// cond_affcts[{0}].append({1})'.format(var, condition))
                    else:
                        self._cond_affcts[var] = [condition]
# affcts.append('// cond_affcts[{0}] = {1}'.format(var, condition))

                else:
                    if self._current_condition is None:
                        self._affcts[var] = True

        if len(affcts) > 0:

            return '\n'.join(affcts) + '\n'

        else:
            return ''

    def _piecewise_matrix_fold(self, expr):
        assert Is(expr).Matrix

        nrows = expr.shape[0]
        ncols = expr.shape[1]

        exprs = dict()

        expr = Matrix(expr.shape[0], expr.shape[1],
                      lambda i, j: piecewise_fold(expr[i, j]))

        for i in range(nrows):
            for j in range(ncols):
                if expr[i, j].is_Piecewise:
                    for e, c in expr[i, j].args:
                        if c in exprs:
                            exprs[c][(i, j)] = e
                        else:
                            exprs[c] = dict()
                            exprs[c][(i, j)] = e

        if len(exprs) > 0:
            m = dict()

            for c in exprs:
                m[c] = Matrix(
                    nrows, ncols,
                    lambda i, j: exprs[c][(i, j)]
                    if (i, j) in exprs[c] else expr[i, j])

            return Piecewise(*((m[c], c) for c in exprs))

        else:
            return expr

    def doprint(self, expr, assign_to=None):

        if self._assertions:
            pre_asserts = '\n'.join(['/*@ assert {0}; */'.format(a) for a in
                                     [r'\is_finite(({0}) {1})'.format(self._value_type, fsym)
                                      for fsym in expr.free_symbols]]) + '\n'
        else:
            pre_asserts = ''

        # avoid bad (cond) ? (...) : (...) sequences
        expr = self._fix_integer_power(expr)

        if expr.is_Matrix:

            expr = Matrix(expr.shape[0], expr.shape[1],
                          lambda i, j: flatten_piecewise(expr[i, j]))

        elif expr.is_Piecewise:
            expr = flatten_piecewise(expr)

        # cse is meaningful only if assignment is specified
        if assign_to is not None:

            if self._do_cse:
                (assignments, substitued_expr) = self._cse(expr)
            else:
                (assignments, substitued_expr) = ([], [expr])

            subs = []
            reworked_expr = substitued_expr[0]

            for (i, (var, val)) in enumerate(assignments):

                # cse bug ?
                # sympy 7.6
                # cse([Matrix([[x**2 + 1.0/fsqrt(x**2 + y**2)]])])
                # => ([], [Matrix([[x**2 + 1.0/fsqrt(x**2 + y**2)]])])
                # cse(Matrix([[x**2 + 1.0/fsqrt(x**2 + y**2)]]))
                # => ([(x0, x**2)], [Matrix([[x0 + 1.0/fsqrt(x0 + y**2)]])])
                reworked_expr = reworked_expr.subs(val, var)

                self._varid[var] = i

                self._decls[var] = val

#                var._assumptions = val._assumptions
                var._assumptions['real'] = True

                if val in self._assignments_values:
                    subs += (var, self._assignments_values[val])
                else:
                    self._assignments_values[val] = var

            if subs != []:
                for var in self._decls:
                    if not self._decls[var].is_Symbol:

                        try:
                            self._decls[var] = self._decls[var].subs(subs)
                        except:
                            pass

        else:
            reworked_expr = expr

        decls = ''
        affcts = ''

        return self.indent_code(pre_asserts + decls + affcts + super(LocalCCodePrinter, self).doprint(
            piecewise_fold(reworked_expr), assign_to))

    def _print_Assignment(self, expr):

        if expr.rhs.is_Matrix:

            clines = dict()
            conds = []

            for i in range(expr.rhs.shape[0]):
                for j in range(expr.rhs.shape[1]):
                    elem = piecewise_fold(expr.rhs[i, j])

                    if elem.is_Piecewise:
                        for e, c in elem.args:
                            if c not in clines:
                                clines[c] = []
                                conds.append(c)
                                clines[c].append(e)
                            else:
                                clines[c].append(e)

            lines = []
            if len(conds) > 0:

                for c in conds:
                    lines.append(self._print_declarations([c]))
                    lines.append(self._print_declarations(clines[c]))
                    lines.append(self._print_affectations([c]))

                if_s = ['if', 'else if', 'else']
                if_s_index = 0

                def affects(cond):
                    return len(set().union(*[self._needed_symbols(e)
                                             for e in clines[cond]])) > 0

                affects_cond = filter(affects, conds)

                laffects_cond = len(affects_cond)

#                if self._assertions:
#                    lines.append('/*@ assert {0}; */'.format(
#                        ' || '.join([self._print(c) for c in conds])))

                for c in affects_cond:
                    lines.append(
                        '{0} ({1})'.format(if_s[if_s_index], self._print(c)))
#                    if_s_index = 1
                    lines.append('{')
                    lines.append(
                        self._print_affectations(clines[c], condition=c))
                    lines.append('}')

                lines.append('')

            return '\n'.join(filter(lambda s: len(s) > 0, lines)) +\
                super(LocalCCodePrinter, self)._print_Assignment(expr)

        else:

# comment = '/* Assignment {0}={1} */'.format(self._print(expr.lhs),
# expr.rhs)
            comment = ''
#            for sa in self._preconditions(expr.lhs, expr.rhs):
#                comment += '{0}\n'.format(sa)

            expr_rhs = piecewise_fold(expr.rhs)

            if expr_rhs.is_Piecewise:

                expressions, conditions = zip(*expr_rhs.args)
                #                    '\n/*@ assert \\valid(&{0}); */\n\n/*@\nrequires \\valid(&{0});\nassigns {1};\nensures \is_finite({0});*/\n'.format(self._print(expr.lhs), self._needed_symbols(expr_rhs)) +\
#                    '\n/*@ assert \\valid(&{0}); */\n\n/*@\nrequires \\valid(&{0});\nassigns {1};\nensures \is_finite({0});*/\n'.format(self._print(expr.lhs), self._needed_symbols(expr.lhs)) +\

                return comment + \
                    self._print_declarations(conditions) +\
                    self._print_affectations(conditions) +\
                    self._print_declarations(expressions) + '\n' +\
                    super(LocalCCodePrinter, self)._print_Assignment(expr) + '\n' +\
                    asserts(self._postconditions(expr.lhs, expr.rhs))
            else:
                expressions = [expr]

                assigns = ''
                if self._assertions:
                    assigns = '/*@ assigns {0}; */\n'.format(
                        self._print(expr.lhs))

                if self._contracts:
                    copen = '/*@\n'
                    cclose = '*/\n'
                    assigns = 'assigns {0};\n'.format(self._print(expr.lhs))
                    pc = self._postconditions(expr.lhs, expr.rhs)
                    if len(pc) > 0:
                        ensures = '\n'.join(
                            ['ensures {0};'.format(e) for e in pc])
                    else:
                        ensures = ''

                    prec = self._postconditions(
                        expr.lhs, expr.rhs, on_lhs=False)
                    if len(prec) > 0:
                        requires = '\n'.join(
                            ['requires {0};'.format(p) for p in prec])
                        pre_asserts = '\n'.join(
                            ['/*@ assert {0};*/'.format(a) for a in prec])
                    else:
                        requires = ''
                        pre_asserts = ''
                else:
                    copen = ''
                    cclose = ''
                    requires = ''
                    ensures = ''
                    pre_asserts = ''
                return comment + \
                    self._print_declarations(expressions) +\
                    self._print_affectations(expressions) +\
                    pre_asserts +\
                    copen +\
                    requires +\
                    assigns +\
                    ensures +\
                    cclose +\
                    super(LocalCCodePrinter, self)._print_Assignment(expr) + '\n' +\
                    asserts(
                        self._postconditions(expr.lhs, expr.rhs))

    def _print_Function(self, expr):
        return '{0}({1})'.format(type(expr), ','.join((self._print(a) for a in expr.args)))

    def _print_Rational(self, expr):
        '''we only change the way rational are output: do not force long
        double temporary hack ...
        '''
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

    def _print_fsqrt(self, expr):

        return self._print(sqrt(expr.args[0]))

    def _print_fmul(self, expr):

        return self._print(Mul(*expr.args))

    def _print_Pow(self, expr):

        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)

        if expr.exp == -1:
            if expr.base.is_Atom:
                return '1.0/{0}'.format(expr.base)
            else:
                return '1.0/({0})'.format(expr.base)

        elif expr.exp < 0:
            return '1.0/({0})'.format(self._print(Pow(expr.base, -expr.exp)))

        elif expr.exp == 0.5:
            return 'sqrt({0})'.format(self._print(expr.base))

        elif Is(expr.exp).integer:

            if expr.base.is_Atom:
                plh = '{0}'
            else:
                plh = '({0})'
            return '*'.join([plh.format(self._print(expr.base))] * expr.exp)

        elif Is(expr.exp).Rational:

            return self._print(Pow(Pow(expr.base, Rational(1, expr.exp.q)),
                                   expr.exp.p, evaluate=False))

        else:

            return 'pow({0}, {1})'.format(self._print(expr.base),
                                          self._print(expr.exp))

    def _print_Float(self, expr):
        # cf
        # http://stackoverflow.com/questions/25222681/scientific-exponential-notation-with-sympy-in-an-ipython-notebook
        return StrPrinter({'full_prec': False}).doprint(Float(expr, 52))

    def _print_Mul(self, expr):

        prec = precedence(expr)

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec) for x in a]
        b_str = [self.parenthesize(x, prec) for x in b]

        if len(b) == 0:
            return sign + '*'.join(a_str)

        elif len(b) == 1:
            assert_p = ''.format(b_str[0])

            return assert_p + sign + '*'.join(a_str) + "/" + b_str[0] + ''
        else:
            assert_p = ''.format('*'.join(b_str))
            return assert_p + sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str) + ''

    def _print_Piecewise(self, expr):

        lines = []

        if expr.has(Assignment):

            try:
                p_expr = piecewise_fold(expr)
            except:
                p_expr = expr

            expressions, conditions = zip(*p_expr.args)

            lconds = len(conditions)

            if self._assertions:
                lines.append('/*@ assert {0}; */'.format(
                    ' || '.join([self._print(c) for c in conditions])))

            for num_cond, cond in enumerate(conditions):

                if num_cond == 0:
                    if_st = 'if ({0})'

                else:
                    if_st = 'else if ({0})'

                lines.append(if_st.format(self._print(cond)))
                lines.append('{')
# lines.append('DEBUG_PRINT("Case ({0}) is True.\\n");'.format(cond))
                lines.append(
                    self._print_affectations([expressions[num_cond]], cond))
                self._current_condition = cond
                lines.append(self._print(expressions[num_cond]))
                self._current_condition = None

                lines.append('}')

            return '\n'.join(filter(lambda s: len(s) > 0, lines))

        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = ["((%s) ? (%s)" % (self._print(c), self._print(e))
                       for e, c in expr.args[:-1]]
            last_line = ": (%s)" % self._print(expr.args[-1].expr)

            return ": ".join(ecpairs) + last_line + " ".join([")" * len(ecpairs)])

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [line.lstrip(' \t') for line in code]

        increase = [int(any(map(line.endswith, inc_token))) for line in code]
        decrease = [int(any(map(line.startswith, dec_token)))
                    for line in code]

        pretty = []
        level = self._level
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (self._tab * (level + 1), line))
            level += increase[n]
        return pretty


def localccode(expr, assign_to=None, tab='    ', level=0,
               array_format='C', epsilon_nan=0., epsilon_inf=0.,
               epsilon_power=1.,
               assertions=False,
               contracts=False, postcheck_hooks=False, do_cse=True,
               settings={}):

    return (LocalCCodePrinter(
        settings=settings,
        tab=tab, level=level, array_format=array_format,
        epsilon_nan=epsilon_nan, epsilon_inf=epsilon_inf,
        epsilon_power=epsilon_power,
        assertions=assertions, contracts=contracts,
        postcheck_hooks=postcheck_hooks, do_cse=do_cse).doprint(
        expr, assign_to))


def files_generation(name, variables, data, result='result',
                     float_type='double'):

    with open('{0}.h'.format(name), 'w') as header:
        header.write(
            """
#ifndef {0}_H
#define {0}_H
void {1}({2}, {3});
#endif
""".format(name.upper(), name, ', '.join(
           ['{0} {1}'.format(float_type, v) for v in variables]),
           '{0}* {1}'.format(float_type, result)))

    with open('{0}.c'.format(name), 'w') as implem:
        implem.write('#include "{0}.h"\n'.format(name))
        implem.write(data)


def funcodegen(name, expr, variables=None, intervals=None,
               maxval=1e6, result='result', float_type='double',
               tab='    ',
               level=0, array_format='C', epsilon_nan=0., epsilon_inf=0.,
               epsilon_power=1.,
               with_files_generation=False,
               assertions=False, contracts=False,
               postcheck_hooks=False, main_check=False, do_cse=True,
               settings={}):

    if variables is None:
        variables_str = dict()
        variables = list()
        for s in expr.free_symbols:
            if not s.__str__() in variables_str:
                variables.append(s)
                variables_str[s.__str__()] = True
        variables = sorted(variables, key=lambda s: s.__str__())

    if intervals is None:
        intervals = {v: Interval(-maxval, maxval) for v in variables}
    else:
        for v in variables:
            if not v in intervals:
                intervals[v] = Interval(-maxval, maxval)

    if not Is(expr).Matrix:
        return funcodegen(name, Matrix([expr]),
                          variables=variables,
                          intervals=intervals,
                          maxval=maxval,
                          result=result,
                          float_type=float_type,
                          tab=tab,
                          level=level,
                          array_format=array_format,
                          epsilon_nan=epsilon_nan,
                          epsilon_inf=epsilon_inf,
                          epsilon_power=epsilon_power,
                          with_files_generation=with_files_generation,
                          assertions=assertions,
                          contracts=contracts,
                          postcheck_hooks=postcheck_hooks,
                          main_check=main_check,
                          do_cse=do_cse,
                          settings=settings)

    if main_check:
        include_check = \
            '#include "funcodegen.h"\n'
        requires_check = \
            '/*@\n' +\
            '\n'.join(['requires ({0} <= {1} <= {2});'.format(print_float(intervals[v].inf), v.__str__(), print_float(intervals[v].sup)) for v in variables]) + '\n' +\
            'assigns {0}[0..{1}];\n'.format(result, (expr.shape[0] * expr.shape[1]) - 1) +\
            '\n'.join(['ensures \is_finite((double) {0}[{1}]);'.format(result, ind) for ind in range(expr.shape[0] * expr.shape[1])]) +\
            '*/\n'
        str_check = \
            '#ifdef __FRAMAC__\n' +\
            'int main()\n{\n' +\
            '\n'.join(['{0}{1} {2} =  Frama_C_double_interval({3}, {4});'.format(
                       tab, float_type, v.__str__(), print_float(intervals[v].inf), print_float(intervals[v].sup)) for v in variables]) +\
            '\n' +\
            '{0}{1} {2}[{3}];'.format(tab, float_type, result, expr.shape[0] * expr.shape[1]) +\
            '\n' +\
            '{0}{1}({2}, {3});'.format(tab, name, ', '.join((v.__str__() for v in variables)), result) +\
            '\n}' +\
            '\n#endif'
    else:
        include_check = ''
        requires_check = ''
        str_check = ''
    code_result = \
        include_check +\
        requires_check +\
        'void ' + name + '(\n' +\
        '\n'.join(['{0}{1} {2},'.format(
            tab, float_type, v.__str__()) for v in variables]) +\
        '\n' +\
        '{0}{1} *{2}'.format(tab, float_type, result) +\
        ')' +\
        '\n{\n' +\
        localccode(expr=expr, assign_to=result, tab=tab, level=level,
                   array_format=array_format, epsilon_nan=epsilon_nan,
                   epsilon_inf=epsilon_inf,
                   epsilon_power=epsilon_power,
                   assertions=assertions,
                   contracts=contracts,
                   postcheck_hooks=postcheck_hooks, do_cse=do_cse,
                   settings=settings) +\
        '\n}\n' + str_check
    if with_files_generation:
        return files_generation(name, variables, code_result, result=result,
                                float_type=float_type)
    else:
        return code_result
