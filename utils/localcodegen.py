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

def append(d,k,v):
    if k in d:
        if not v in d[k]:
            d[k].append(v)
    else:
        d[k] = list()
        d[k].append(v)

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
            yield Symbol('{0}{1}'.format(self._name, self._counter))



class LocalCCodePrinter(CCodePrinter):


    def __init__(self, settings={}, tab='    ', level=0, array_format='C'):
        CCodePrinter.__init__(self, settings)
        self._some_vars = SomeVars()
        self._value_type = 'double'
        self._relational_type = 'int'
        self._assignments_values = dict()
        self._decls = dict()
        self._varid = dict()
        self._declared = dict()
        self._affcts = dict()
        self._tab = tab
        self._level = level
        self._do_cse = True
        self._array_format = array_format
        self._inside_condition = 0


    def _print_MatrixElement(self, expr):
        if self._array_format == 'C':
            return "{0}[{1}]".format(expr.parent, expr.j +
                                     expr.i*expr.parent.shape[1])

        elif self._array_format == 'Fortran':
             return "{0}[{1}]".format(expr.parent, expr.i +
                                      expr.j*expr.parent.shape[0])

    def _traverse_matrix_indices(self, mat):

        rows, cols = mat.shape
        if self._array_format == 'C':
            return ((i, j) for i in range(rows) for j in range(cols))
        elif self._array_format == 'Fortran':
            return ((i, j) for j in range(cols) for i in range(rows))

    def _needed_symbols(self, expr):

        l = set()

        p_expr = expr

        all_conditions = []
        for i, subexpr in enumerate(postorder_traversal(p_expr)):
            if Is(subexpr).Piecewise:
                all_conditions += [c for _, c in subexpr.args]

        if len(all_conditions)>0:
            symbols = set([item for conditions in all_conditions
                           for item in conditions.free_symbols])
        else:
            symbols = expr.free_symbols

        for symb in symbols:
            if symb in self._decls:
                l.add(symb)
                l = l.union(self._needed_symbols(self._decls[symb]))
        return l

    def _print_declarations(self, lexpr):

        l = set()
        decls = []
        for expr in lexpr:
            l = l.union(self._needed_symbols(expr))

        for var, _ in sorted([(s, self._varid[s]) for s in l], key=lambda t: t[1]):
            if var in self._decls and not var in self._declared:
                if self._decls[var].is_Relational or isinstance(self._decls[var], (Logic, And, Or, Not)):
                    decls.append('{0} {1};'.format(self._relational_type, var))
                else:
                    decls.append('{0} {1};'.format(self._value_type, var))

                if self._inside_condition == 0:
                    self._declared[var] = True

        if len(decls)>0:
            return '\n'.join(decls) + '\n'
        else:
            return ''

    def _print_affectations(self, lexpr):

        l = set()
        affcts = []
        for expr in lexpr:
            l = l.union(self._needed_symbols(expr))


        for var, _ in sorted([(s, self._varid[s]) for s in l], key=lambda t: t[1]):
            if var in self._decls and not var in self._affcts:
                affcts.append('{0}={1};'.format(var, super(LocalCCodePrinter, self)._print(self._decls[var])))

            if self._inside_condition == 0:
                self._affcts[var] = True

        if len(affcts)>0:
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

        if len(exprs)>0:
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

        # avoid bad (cond) ? (...) : (...) sequences
        if Is(expr).Matrix:
            expr = Matrix(expr.shape[0], expr.shape[1],
                          lambda i, j: piecewise_fold(expr[i, j]))


        # cse is meaningful only if assignment is specified
        if self._do_cse and assign_to is not None:

            (assignments, substitued_expr) = cse(expr, self._some_vars)

            subs = []

            for (i, (var, val)) in enumerate(assignments):

                self._varid[var] = i
                self._decls[var] = val

                if val in self._assignments_values:
                    subs += (var, self._assignments_values[val])
                else:
                    self._assignments_values[val] = var

            for var in self._decls:
                if not self._decls[var].is_Symbol:

                    try:
                        self._decls[var] = self._decls[var].subs(subs)
                    except:
                        pass

            # idem...
            try:
                reworked_expr = substitued_expr[0].subs(subs)
            except:
                reworked_expr = substitued_expr[0]

            decls = ''
            affcts = ''

        else:
            reworked_expr = expr
            decls = ''
            affcts = ''

        return self.indent_code(decls + affcts + super(LocalCCodePrinter, self).doprint(
            piecewise_fold(reworked_expr), assign_to))


    def _print_Assignment(self, expr):

        comment = '// Assignment {0}={1}\n'.format(expr.lhs, expr.rhs)

        if not expr.rhs.is_Piecewise:

            if expr.rhs.is_Matrix:
                return super(LocalCCodePrinter, self)._print_Assignment(expr)

            else:
                return comment + self._print_declarations([expr.rhs]) +\
                    self._print_affectations([expr.rhs]) +\
                    super(LocalCCodePrinter, self)._print_Assignment(expr)

        else:

         conditions = [c for _, c in expr.rhs.args]

         return comment + self._print_declarations(conditions) +\
                self._print_affectations(conditions) +\
                super(LocalCCodePrinter, self)._print_Assignment(expr)



    def _print_Rational(self, expr):
        '''we only change the way rational are output: do not force long
        double temporary hack ...
        '''
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

    def _print_Pow(self, expr):

        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)

        if expr.exp == -1:
            return '1.0/(assert(NOT_ZERO({0})), {0})'.format(self.parenthesize(expr.base, PREC))

        elif expr.exp == 0.5:
            return '(assert(POSITIVE({0})), sqrt({0}))'.format(self._print(expr.base))

        elif expr.exp == 2:

            return '{0}*{0}'.format(self._print(expr.base))

        elif expr.exp < 1 and expr.exp > -1:

            return '(assert(POSITIVE({0}), pow({0}, {1})'.format(self._print(expr.base),
                                          self._print(expr.exp))

        else:
            return 'pow(%s, %s)' % (self._print(expr.base),
                                 self._print(expr.exp))


    def _print_Float(self, expr):
        # cf http://stackoverflow.com/questions/25222681/scientific-exponential-notation-with-sympy-in-an-ipython-notebook
        return StrPrinter({'full_prec': False}).doprint(Float(expr, 20))

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
            assert_p = '(assert(NOT_ZERO({0})), '.format(b_str[0])

            return assert_p + sign + '*'.join(a_str) + "/" + b_str[0] + ")"
        else:
            assert_p = '(assert(NOT_ZERO({0})),'.format('*'.join(b_str))
            return assert_p + sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str) + ")"

    def _print_Piecewise(self, expr):

        lines = []

        self._inside_condition += 1

        if expr.has(Assignment):

            try:
                p_expr = piecewise_fold(expr)
            except:
                p_expr = expr

            expressions, conditions = zip(*p_expr.args)

            for num_cond, cond in enumerate(conditions):

                if num_cond == 0:
                    if_st = 'if'
                elif num_cond == 1 and cond is True:
                    if_st = 'else'
                else:
                    if_st = 'else if'

                lines.append(if_st + ' ({0})'.format(self._print(cond)))
                lines.append('{')
                lines.append('DEBUG_PRINT("Case ({0}) is True.\\n");'.format(cond))

                lines.append(self._print(expressions[num_cond]))

                lines.append('}')

            self._inside_condition -= 1
            return "\n".join(lines)

        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = ["((%s) ? (%s)" % (self._print(c), self._print(e))
                    for e, c in expr.args[:-1]]
            last_line = ": (%s)" % self._print(expr.args[-1].expr)

            self._inside_condition -= 1
            return ": ".join(ecpairs) + last_line + " ".join([")"*len(ecpairs)])


    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [ line.lstrip(' \t') for line in code ]

        increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
        decrease = [ int(any(map(line.startswith, dec_token)))
                     for line in code ]

        pretty = []
        level = self._level
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (self._tab*level, line))
            level += increase[n]
        return pretty



def localccode(expr, assign_to=None, tab='    ', level=0,
               array_format='C', **settings):

    return (LocalCCodePrinter(tab=tab, level=level, array_format=array_format, settings=settings).doprint(expr, assign_to))
