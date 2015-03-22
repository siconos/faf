# latex generation tools :

# - dump Piecewise function with conditions factorization

from sympy import *
from localcodegenold import set_var, piecewise_fold

def print_latex_by_conditions(expr, result):
    conds = dict()
    conds[True]=[]
    conds_var=dict()
    if hasattr(expr,'shape'):
        (nrow,ncol) = expr.shape
        for i in range(nrow):
            for j in range(ncol):
                set_var(result(i,j),piecewise_fold(expr[ncol*i+j]),conds, conds_var)

    for var,val in conds[True]:
        print "${0}={1}$".format(var,latex(val))

    for cond in list(conds):

        if hasattr(cond,'index') and not isinstance(cond,str):
            assert (not cond[0].is_Piecewise)

            if len(cond) ==1:
                print r"\subsection{$ \text{for}\: %s$}" % latex(cond[0])
            else:
                s = r"\subsection{$ \text{for}\: %s "% latex(cond[0])
                for c in cond[1:]:
                    s = s + r"\:\text{and}\:%s" % latex(c)
                    s = s + r"$}"
                print s
            for iexpr in conds[cond]:
                assert not iexpr[1].is_Piecewise
                print r"\paragraph{}$%s=%s$" % (iexpr[0],latex(iexpr[1]))
