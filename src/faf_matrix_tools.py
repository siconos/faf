from faf_timeout import *
from faf_tools import *

import numpy as np
import sys
from scipy.sparse import csr_matrix, linalg, diags

from numpy.linalg import matrix_rank,svd
try:
    from scipy.linalg.interpolative import estimate_rank
except:
    def estimate_rank(*args):
        print("warning : estimate_rank undefined")

try:
    from scipy.sparse.linalg import lsmr
except:
    def lsmr(*args):
        print("warning : lsmr undefined")

try:
    from scipy.sparse.linalg import svds
except:
    def svds(*args):
        print("warning : svds undefined")



@timeout(100)
def dense_matrix_rank(M):
    return matrix_rank(M)

@timeout(1000)
def sparse_matrix_svd(A,k):
    return svds(A,k)

@timeout(1000)
def dense_matrix_rank_estimate(A,tol):
    return estimate_rank(A,tol)



# estimate of condition number and norm from lsmr
# http://www.stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf
@timeout(200)
def _norm_cond(A):
    print('in norm cond')
    #problem = read_fclib_format(problem_filename)[1]
    #A = csr_matrix(N.SBM_to_sparse(problem.M)[1])
    #A = csr_matrix(N.SBM_to_sparse(M)[1])
    
    print("A=", A)
    print("A.shape", A.shape)
    print("Computev lsmr ...")
    r = lsmr(A, np.ones([A.shape[0], 1]))  # solve Ax = 1
    norm_lsmr=r[5]
    cond_lsmr=r[6]
    print("norm_lsr=", norm_lsmr)
    print("cond_lsr=", cond_lsmr)
    #print "r=", r
    try:
        print ("Computev svds(A,1) ...")
        _svd = svds(A,1)[1]
        eps = sys.float_info.epsilon
        tol = _svd.max() * max(A.shape) * eps
    except Exception as e :
        print ("-->   svds failed to compute the maximum singular value")
        print ("-->" ,  e)
        _svd = [1.0]
        eps = sys.float_info.epsilon
        tol = max(A.shape) * eps
    #print tol
    #print "============"
    # from scipy.sparse.linalg import LinearOperator
    # def mv(v):
    #     return A.dot(v)
    # A_LO = LinearOperator( A.shape, matvec=mv )
    # print isinstance(A_LO, LinearOperator)
    rank_estimate = np.nan
    # try:
    #     print("Compute rank estimate ...")
    #     print("dense ...")
    #     B= A.todense()
    #     print("dense ...")
    #     rank_estimate=dense_matrix_rank_estimate(B,tol)
    #     #rank_estimate=estimate_rank(A.todense(), tol)
    # except Exception as e :
    #     print ("--> rank_estimate", e)
    # print("rank_estimate", rank_estimate)

    # #print "svd dense method", svd(A.todense())[1]

    rank_dense = np.nan
    # try:
    #     print("Compute rank dense ...")
    #     rank_dense = dense_matrix_rank(A.todense())
    # except Exception as e :
    #     print ("--> dense_matrix_rank", e)
    # print ("rank_dense", rank_dense)

    if not np.isnan(rank_estimate):
        k=min(rank_estimate,A.shape[0]-1)
    else:
        k = A.shape[0]-1
    try:
        print("Compute svds(A,k) for  ",k," singular values  ...")
        _svd = sparse_matrix_svd(A,k)[1]
        #print "_svd",_svd
    except Exception as e :
        print("--> sparse_matrix_svd  failed to compute ",k," singular values")
        print("--> sparse_matrix_svd " ,  e)

        _svd =[]
    # compute rank with http://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.linalg.matrix_rank.html
    rank_svd=np.nan
    nonzero_sv=[]
    import math
    for sv in _svd:
        #print sv,tol
        if (sv >=tol and (not math.isnan(sv))):
            nonzero_sv.append(sv)
    nonzero_sv.sort(reverse=True)
    #print(nonzero_sv)
    if (len(nonzero_sv) >0):
        rank_svd = len(nonzero_sv)
    else:
        rank_svd = np.nan
    #print "rank_svd", rank_svd

    if not math.isnan(rank_dense):
        rank=rank_dense
    else:
        if not math.isnan(rank_svd):
            rank=rank_svd
        else :
            rank=rank_estimate

    if not math.isnan(rank):
        nonzero_sv = nonzero_sv[0:rank]

    if (len(nonzero_sv) >0):
        max_nonzero_sv= max(nonzero_sv)
        min_nonzero_sv= min(nonzero_sv)
        ratio_max_min_nonzero_sv =  max(nonzero_sv)/min(nonzero_sv)
    else:
        max_nonzero_sv= np.nan
        min_nonzero_sv= np.nan
        ratio_max_min_nonzero_sv =  np.nan
    # http://personales.unican.es/beltranc/archivos/CNmatricesF.pdf
    return norm_lsmr, cond_lsmr, max_nonzero_sv, min_nonzero_sv,  ratio_max_min_nonzero_sv, rank, rank_dense, rank_svd, rank_estimate

norm_cond = Memoize(_norm_cond)


def _cond_problem(filename):
    problem = read_fclib_format(filename)[1]
    return float(problem.numberOfContacts * 3) / float(numberOfDegreeofFreedom(filename))

cond_problem = Memoize(_cond_problem)

def _cond(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['cond']
        except:
            r = np.nan
    print ("cond=",r)
    return r

cond = Memoize(_cond)
def _cond_lsmr(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['cond_lsmr']
        except:
            r = np.nan
    print ("cond_lsmr=",r)
    return r
cond_lsmr = Memoize(_cond_lsmr)

def _rank_dense(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['rank_dense']
        except:
            r = np.nan
    print ("rank_dense=",r)
    return r

rank_dense = Memoize(_rank_dense)

def _rank_estimate(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['rank_estimate']
        except:
            r = np.nan
    print ("rank_estimate=",r)
    return r

rank_estimate = Memoize(_rank_estimate)

def dd(X):
    Sum_values_in_given_row = np.ravel((np.sum(abs(X),axis=1)))

    X_diag = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        X_diag[i] = abs(X[i,i])
    #print('X_diag', X_diag)

    test = 2*X_diag - Sum_values_in_given_row
    print('test',test) 
    
    if np.all( test >= 0.0):
        print('matrix is diagonally dominant')
        return True, min(test)
    else:
        print('NOT diagonally dominant')
        return False, min(test)


def _test_symmetry_W(problem_filename):
    problem = read_fclib_format(problem_filename)[1]
    A = csr_matrix(N.SBM_to_sparse(problem.M)[1])

    #print("A=", A)
    #print("A=", np.transpose(A))
    #print("A-A^T",A-np.transpose(A))
    print("norm A-A^T", linalg.norm(A-np.transpose(A)) )
    is_symmetric =False
    symmetry_test = linalg.norm(A-np.transpose(A),np.inf)
    print("symmetry_test === ", symmetry_test)
    
    is_symmetric = (symmetry_test <= 1e-12)
    is_dd, dd_test =  dd(A)
    print("is_dd", is_dd)
    print("is_symmetric", is_symmetric)
    return is_symmetric, symmetry_test, is_dd, dd_test
    
test_symmetry_W = Memoize(_test_symmetry_W)
