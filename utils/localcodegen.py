# c code generation tools :

# - calls deeper common subexpression elimination
# - code from Piecewise function with conditions factorization
# - some pow(x,y) optimization

from sympy import *
from sympy.utilities.iterables import postorder_traversal

#
# some dict operations
#
def get(d,k):
    try:
        return d[k]
    except KeyError:
        return None

def append(d,k,v):
    if k in d:
        if not v in d[k]:
            d[k].append(v)
    else:
        d[k] = list()
        d[k].append(v)

def append_list(d,k,l):
    if k in d:
        for v in l:
            if not v in d[k]:
                d[k].append(v)
    else:
        d[k] = l

#
# add dependent conditions
#
def add_depend_conditions(var, expr, conds_var):
    for subtree in [ str(s) for s in postorder_traversal(expr) ]:
        if subtree in conds_var:
            append_list(conds_var, var, conds_var[subtree])

#
# assemble conditions (with logical 'and' operator) so they can be
# dumped in groups
#
def conditions_fold(var, cond, and_conds, conds):
    # cond is (expr = cond[0], condition = cond[1]). If expr is
    # Piecewise assemble conditions with a recursion
    if cond[0].is_Piecewise:
        for subcond in cond[0].__getnewargs__():
            conditions_fold(var,subcond,[cond[1]] + and_conds, conds)
    else:
        # else append to key = terminal conditions (tand_conds) (var,
        # expr=cond[0])
        tand_conds = and_conds+[cond[1]]
        append(conds, tuple(tand_conds), (var, cond[0]))
        append(conds, "all", (tuple(tand_conds)))

def conditions_dependencies(list_var, conds_var, conds):
    from_conds = dict()
    conds_ordered = list()
    for var in [ str(v) for v in list_var ]:
        if var in conds_var:
            for lc in conds_var[var]:
                print conds_var[var],lc
    return from_conds, conds_ordered

def dump_conditions(done, conds, from_conds, list_kconds, remove_conds, tab, nbops):

    for icond in list_kconds:

        if icond is True:
            continue

        if icond not in conds:
            continue

        if icond in done:
            continue

        done[icond]=True

        lcond = list(icond)

        for c in remove_conds:
            lcond.remove(c)
        cond = tuple(lcond)

        if len(cond) == 1:
            print tab + "if ({0})".format(local_ccode(nbops, cond[0]))
            print tab + "{"

        else:
            if len(cond)>1:
                s = tab + "if (({0}) ".format(local_ccode(nbops, cond[0]))
                for c in cond[1:]:
                    s = s + "&& ({0})".format(local_ccode(nbops, c))
                    s = s + ")"
                print s
                print tab + "{"

        for iexpr in conds[icond]:
            assert not iexpr[1].is_Piecewise
            print tab + tab + "{0}={1};".format(iexpr[0],local_ccode(nbops, iexpr[1]))

        for _cond in cond:
            if _cond in from_conds:
                dump_conditions(done, conds, from_conds, from_conds[_cond], [_cond], tab + "  ", nbops)

        if len(cond)>=1:
            print tab + "};"

#
# something = value in dump code
#
def set_var(var, expr, conds, conds_var, declared=None):

    if declared is not None:
        declared[var] = True

    iexpr = piecewise_fold(expr)

    add_depend_conditions(var, iexpr, conds_var)

    if var in conds_var:
        iexpr = Piecewise(*[ (iexpr, c) for c in conds_var[var] ])

    if iexpr.is_Piecewise:
        for cond in iexpr.__getnewargs__():
            append(conds_var, var, cond[1])
            conditions_fold(var, cond, [], conds)
            if declared is not None:
                conds[True].append((var,None))
    else:
        assert iexpr is not None
        conds[True].append((var,iexpr))

#
# Pow(x,n) -> x*x*x ... *x  n times
#
def intpow2mult(v_i,n_i):
    v = sympify(v_i)
    n = sympify(n_i)
    assert (n.is_Integer)
    if n==0:
        return "1."
    N = abs(n)
    if n<0:
        s = "1./("+ccode(v)
    else:
        s = ccode(v)
    for i in range(0,N-1):
        s += "*" + ccode(v)
    if n<0:
        s += ")"

    return s


def local_count_op(c, expr):

    # .count_ops() does not work on relationals (bug?)

    count = 0
    symbs = numbered_symbols()
    nexpr = expr
    args = list(expr.args)
    for i,e in enumerate(args):
        if not e.is_Number:
            sym = symbs.next()
            nexpr = nexpr.subs(e,sym)
            for j in range(i+1, len(args)):
                args[j] = args[j].subs(e,sym)

    for subexpr in [nexpr] + list(expr.args):

        if subexpr.args == ():
            continue

        if subexpr.is_Function:
            count += Symbol(str(type(subexpr)).upper())
            continue

        if subexpr.is_Mul:
            if subexpr.args[0] == -1:
                if len(subexpr.args)>0:
                    count += Symbol("MUL")
            else:
                count += Symbol("MUL")
            continue

        if hasattr(subexpr, "get_relational_class"):
            count += reduce(local_count_op,subexpr.args,c)
            continue
        else:
            count += subexpr.count_ops()



    return c+count

#
# apply pow transformation
#
def local_ccode(count, expr):



    if expr.is_Pow:
        (v,exp) = expr.args
        if exp.is_Integer and v.is_Symbol:
            count[0] += Symbol('MUL') * (exp-1)
            return intpow2mult(v,exp)
        else:
            if not exp.is_Integer and (exp*Rational(2)).is_Integer:
                count[0] += Symbol('SQRT')
                return "sqrt("+local_ccode(count, expr**2) +")"
            else:
                count[0] = local_count_op(count[0], expr)
                return ccode(expr)
    else:
        count[0] = local_count_op(count[0], expr)
        return ccode(expr)



#
# 2D output array layout can be c, fortran or i,j
#
def output_format(nrow,ncol,format):
    if format=='fortran':
        return lambda i,j: "{0}".format(i+j*nrow)
    if format=='c':
        return lambda i,j: "{0}".format(i*ncol+j)
    if format=='ij':
        return lambda i,j: "{0},{1}".format(i,j)
    assert(False)


#
# dump c code output is a 2D array
#
def dump_ccode(expr, array_format='ij', result_open='result[', result_close=']'):

    # all conditions : key = assembled conditions, value = list of value expressions
    conds = dict()

    # always true init
    conds[True]=list()

    # variables conditions : key = variable, value = list of conditions
    conds_var = dict()

    # all declared variables
    declared = dict()
    list_var = list()
    nbops = [0]

    # operate on matrix
    if not hasattr(expr, 'shape'):
        return dump_ccode(Matrix([[expr]]),
                          array_format,
                          result_open,
                          result_close)

    (w,g) = cse(expr)

    for sub in w:
        set_var(str(sub[0]), piecewise_fold(sub[1]), conds, conds_var, declared)

    assert hasattr(expr, 'shape')

    (nrow, ncol) = expr.shape
    oformat = output_format(nrow, ncol, array_format)

    for i in range(nrow):
        for j in range(ncol):
            current_expr = piecewise_fold(g[0][ncol*i+j])
            set_var("{0}{1}{2}".format(result_open,oformat(i,j),result_close), current_expr, conds, conds_var)


    none_dumped = dict()
    is_relational = dict()

    for var,val in conds[True]:

        for cond in conds:
            for x in conds[cond]:
                if var == x[0]:
                    if x[1] is not None:
                        is_relational[var] = x[1].is_Relational


        if val is None:
            if not var in none_dumped:
                if var in is_relational and is_relational[var]:
                    print "  int {0} = 0;".format(var)
                else:
                    print "  double {0};".format(var)
                none_dumped[var] = True

        else:
            assert (not val.is_Piecewise)
            if not hasattr(val,"cond") and not hasattr(val, "lhs"):  # not Relational and not ExprCondPair
                if var in declared:
                    print "  double {0}={1};".format(var,local_ccode(nbops, val))
            else:
                print     "  int    {0}={1};".format(var, local_ccode(nbops, val))


    from_conds, conds_ordered = conditions_dependencies(list_var, conds_var, conds)

    done = dict()



    dump_conditions(done, conds, from_conds, conds["all"], (), "  ", nbops)

#    dump_conditions(done, conds, from_conds, list( set(list(conds))-set(list(from_conds))), (), "  ", nbops)

    for var,val in conds[True]:
        if val is not None:
            if var not in declared:
                print "  {0}={1};".format(var,local_ccode(nbops, val))

#    print "  /* operations : {0} */".format(nbops[0])
