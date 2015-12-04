#!/usr/bin/env python

# parallel usage :
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--files={}'

#
# comp.py --max-problems=10 --no-compute --no-collect # output problems.txt
# cat problems.txt | parallel comp.py --timeout=100 --no-collect '--files={}'
#
#



from glob import glob
from itertools import product
import numpy as np
import random
import siconos.numerics as N
N.setNumericsVerbose(0)
import siconos.fclib as FCL

from numpy.linalg import matrix_rank,svd

try:
    from scipy.linalg.interpolative import estimate_rank
except:
    pass

from scipy.sparse import csr_matrix

try:
    from scipy.sparse.linalg import lsmr
except:
    def lsmr(*args):
        print "lsmr undefined"

try:
    from scipy.sparse.linalg import svds
except:
    def svds(*args):
        print "svds undefined"

from subprocess import check_call

import os

import multiprocessing
import multiprocessing.pool
import time


import h5py
import getopt
import sys
import hashlib
import shlex

#print os.path.join(os.path.dirname(sys.argv[0]), 'external/build')

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), 'external/build'))

try:
    import BogusInterface
except:
    pass
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.INFO)


class Memoize():
    def __init__(self, fun):
        self._fun = fun
        self._done = dict()

    def __call__(self, *args):
        if args in self._done:
            return self._done[args]
        else:
            try:
                r = self._fun(*args)
                self._done[args] = r
                return r
            except Exception as e:
                self._done[args] = e
                return e

class WithCriterium():

    def __init__(self, condmin, condmax):
        self._condmin = condmin
        self._condmax = condmax

    def __call__(self, filename):
        r = cond_problem(filename)
        return r > self._condmin and r < self._condmax

def list_from_file(filename):
    with open(filename, 'r') as f:
        return f.read().lstrip().rstrip().split('\n')

def subsample_problems(filenames, proba, maxp, cond, overwrite=False):


    def addext(f):
        if f.endswith('.hdf5'):
            return f
        else:
            return '{0}.hdf5'.format(f)
        
    _filenames = [addext(f) for f in filenames]

    if proba is not None:
        __r = random.sample(_filenames, int(len(_filenames) * proba))
    else:
        __r = _filenames

    if maxp is not None:
        _r = random.sample(__r, min(maxp, len(__r)))
        print _r

    else:
        _r = __r

    if cond is not None:
        r = filter(WithCriterium(cond[0], cond[1]), _r)
    else:
        r = _r

    # overwrite
    if overwrite:
        with open('problems.txt','w') as problems_txt:
            problems_txt.write('{0}\n'.format('\n'.join(r)))

    return r

def extern_guess(problem_filename, solver_name, iteration, h5file):
    data = h5file['data']
    comp_data = data['comp']

    reaction = comp_data[solver_name][problem_filename]['reactions'][iteration]
    velocity = comp_data[solver_name][problem_filename]['velocities'][iteration]
    return reaction, velocity


def split(s, sep, maxsplit=-1):

    try:
        return [shlex.split(kw)[0]
                 for kw in s.split(sep, maxsplit)]
    except Exception:
        sys.stderr.write('don\'t know how to split {0} with {1}\n'
                         .format(s, sep))
        return None


measure = 'flop'
clean = False
compute = True
display = False
display_convergence = False
display_distrib = False
display_distrib_var = False
no_matplot=False
gnuplot_profile = False
logscale=False
gnuplot_distrib = False
gnuplot_with_color = True
output_dat=False
user_filenames = []
user_solvers = []
user_solvers_exact = []
utimeout = 10
keep_files = False
output_errors = False
output_velocities = False
output_reactions = False
measure_name = 'flpops'
ask_compute = True
ask_collect = True
maxiter = 1000000
maxiterls = 10
precision = 1e-8
domain = np.arange(1, 100, .1)
ref_solver_name = 'NonsmoothGaussSeidel'
random_sample_proba = None
max_problems = None
cond_nc = None
with_guess = True
with_mumps = 0
file_filter=None
gnuplot_separate_keys = False
list_contents=False
compute_hardness = False
compute_cond_rank = False
adhoc= False
def usage():
  print "\n \n"
  print 'Usage: '+sys.argv[0]+'[option]'
  print "Options : "
  print " --help "   
  print "   display this message"
  print " --verbose "   
  print "   enable verbose mode equal to 1 for Siconos Numerics"
  print " --no-collect "   
  print "   leave the result into separate file that are named according the solver and the name of the problem"
  print " --just-collect"
  print "   collect all the result into comp.hdf5"
  print " --timeout=n"  
  print "   set the maximum time of computation for each problem to n seconds (default",utimeout,"s)"
  print " --maxiter=n"  
  print "   set the maximum number of iterations for each problem to n (default",maxiter,")"
  print " --maxiterls=n"  
  print "   set the maximum number of iterations for each problem to n (default",maxiterls,")"
  print " --domain='a:d:b'"
  print "   restrict the domain of the performance profile to the interval [a,b] with a step of d (default",domain[0],":",domain[1]-domain[0],":",domain[-1]+domain[1]-domain[0],")"
  print "   or a perfomance profile a should be greater or equal 1"
  print " --measure=value"
  print "   select the value  as the measure for the perfomance profile. Possible values are time, iter, flpops"
  print " --display"
  print "   perform the computation of performance profile and display it in matplotlib"
  print " --display-distrib='from-files' or "
  print "   perform the computation of distribution and display it in matplotlib"
  print " --new"
  print "   remove comp.hdf5 file"
  print " --solvers=string"
  print "   use keyworks in s separated by comma for filtering solvers"
  print " --solvers-exact=string"
  print "   use exact names of solvers in s separated by comma for filtering solvers"
  print " --with-mumps"
  print "   use mumps as linear system solver"
  print " --max-problems=<max>"
  print "   Randomly select <max> problems in current directory." 
  print "   The problems list is written in problems.txt file"
  print " --gnuplot-profile"
  print "   output gnuplot command file profile.gp for plotting profiles woth gnuplot" 
  print " --gnuplot-distrib"
  print "   output gnuplot command file distrib.gp for plotting distribution woth gnuplot" 
  print " --gnuplot-separate-keys"
  print "   output keys anf legend for gnuplot in a separate file."
  print " --display-distrib='from-files' "
  print " --list-contents"
  print "   list contents of comp.hdf5 file"
  print " --compute-cond-rank"
  print "   compute the rank (vairous numerical methods) and condition number of W and store it in the problem file"
  print " --compute-hardness"
  print "   compute the average performance of the best solver on a set of problem divided by the average number of contact"
  
  

  print " Other options have to be documented" 
  print " "
  print " Usage examples:"

  toto = """
  comp.py   --display --time --domain='1:0.1:10'  comp.hdf5
  
  comp.py --display --measure=time --solvers=Gauss,Tresca,SOCLCP,ACLM --domain=1:0.1:100
  """
  print toto
  

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'flop', 'iter', 'time', 'verbose=','no-guess',
                                    'clean', 'display', 'display-convergence','no-matplot',
                                    'files=', 'solvers-exact=', 'solvers=',
                                    'random-sample=', 'max-problems=',
                                    'timeout=', 'maxiter=', 'maxiterls=', 'precision=',
                                    'keep-files', 'new', 'errors',
                                    'velocities', 'reactions', 'measure=',
                                    'just-collect', 'cond-nc=', 'display-distrib=',
                                    'no-collect', 'no-compute', 'domain=',
                                    'replace-solvers-exact=','replace-solvers=',
                                    'gnuplot-profile','gnuplot-distrib', 'logscale', 'gnuplot-separate-keys',
                                    'output-dat', 'with-mumps', 'file-filter=',
                                    'list-contents',
                                    'add-precision-in-comp-file','add-timeout-in-comp-file',
                                    'compute-cond-rank','compute-hardness','adhoc'])


except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)
for o, a in opts:
    if o == '--verbose':
        N.setNumericsVerbose(int(a))
    if o == '--help':
        usage()
        exit(2)
    elif o == '--timeout':
        utimeout = float(a)
    elif o == '--maxiter':
        maxiter = int(a)
    elif o == '--maxiterls':
        maxiterls = int(a)
    elif o == '--precision':
        precision = float(a)
    elif o == '--clean':
        clean = True
    elif o == '--display':
        display = True
        compute = False
    elif o == '--list-contents':
        list_contents = True
        compute = False
    elif o == '--display-convergence':
        display_convergence = True
        compute = False
    elif o == '--measure':
        measure_name = a
    elif o == '--random-sample':
        if os.path.exists('problems.txt'):
            os.remove('problems.txt')
        random_sample_proba = float(a)
    elif o == '--max-problems':
        if os.path.exists('problems.txt'):
            os.remove('problems.txt')
        max_problems = int(a)
    elif o == '--keep-files':
        keep_files = True
    elif o == '--errors':
        output_errors = True
    elif o == '--velocities':
        output_velocities = True
    elif o == '--reactions':
        output_reactions = True
    elif o == '--just-collect':
        ask_compute = False
    elif o == '--no-collect':
        ask_collect = False
    elif o == '--no-compute':
        ask_compute = False
    elif o == '--cond-nc':
        cond_nc = [float (x) for x in split(a,':')]
    elif o == '--display-distrib':
        display_distrib = True
        compute = False
        display_distrib_var = a
    elif o == '--domain':
        urange = [float (x) for x in split(a,':')]
        domain = np.arange(urange[0], urange[2], urange[1])
    elif o == '--no-matplot':
        no_matplot=True
    elif o == '--solvers':
        user_solvers = split(a, ',')
    elif o == '--solvers-exact':
        user_solvers_exact = split(a, ',')
    elif o == '--replace-solvers-exact':
        replace_solvers = split(a, ',')
        try:
            with h5py.File('comp.hdf5','r+') as comp_file:
                solver_in_compfile =  list(comp_file['data']['comp'])
                #print "list(comp_file['data']['comp'])",  solver_in_compfile
                replace_solver_in_compfile = filter(lambda s: any(us == s for us in replace_solvers), solver_in_compfile)
                print "replace solver in comp file", replace_solver_in_compfile
                for s in replace_solver_in_compfile:
                    del comp_file['data']['comp'][s]
        except Exception as e:
            print e
    elif o == '--replace-solvers':
        replace_solvers = split(a, ',')
        #print "replace_solvers",  replace_solvers
        try:
            with h5py.File('comp.hdf5','r+') as comp_file:
                solver_in_compfile =  list(comp_file['data']['comp'])
                #print "list(comp_file['data']['comp'])",  solver_in_compfile
                replace_solver_in_compfile = filter(lambda s: any(us in s for us in replace_solvers), solver_in_compfile)
                print "replace solver in comp file", replace_solver_in_compfile
                for s in replace_solver_in_compfile:
                    del comp_file['data']['comp'][s]
        except Exception as e:
            print e
    elif o == '--gnuplot-profile':
        gnuplot_profile=True
    elif o == '--logscale':
        logscale=True
    elif o == '--gnuplot-distrib':
        gnuplot_distrib=True
    elif o == '--gnuplot-separate-keys':
        gnuplot_separate_keys = True
    elif o == '--output-dat':
        output_dat=True
    elif o == '--with-mumps':
        with_mumps=1
    elif o == '--new':
        try:
            os.remove('comp.hdf5')
        except:
            pass

    elif o == '--files':

        files = split(a, ',')

        for f in files:

            if os.path.exists(f):
                user_filenames += [f]
            else:
                if os.path.exists('{0}.hdf5'.format(f)):
                    user_filenames += ['{0}.hdf5'.format(f)]
    elif o == '--file-filter':
        file_filter=split(a, ',')
        
    elif o == '--no-guess':
        with_guess = False
    elif o == '--add-precision-in-comp-file':
        with h5py.File('comp.hdf5','r+') as comp_file:
            create_attrs_precision_in_comp_file(comp_file,float(a))
    elif o == '--add-timeout-in-comp-file':
        with h5py.File('comp.hdf5','r+') as comp_file:
            create_attrs_timeout_in_comp_file(comp_file,float(a))
    elif o == '--compute-hardness':
        compute_hardness = True
        compute = False
    elif o == '--compute-cond-rank':
        compute_cond_rank = True
        compute = False
    elif o == '--adhoc':
        adhoc = True
        compute = False
    

from ctypes import cdll, c_float, c_longlong, byref
try:
    with_papi = True
    papi=cdll.LoadLibrary('/usr/local/lib/libpapi.so')
except:
    try:
        with_papi = True
        papi=cdll.LoadLibrary('/usr/lib/x86_64-linux-gnu/libpapi.so.5.3.0.0')
    except:
        try:
            with_papi = True
            papi=cdll.LoadLibrary('/home/bremond/faf/install/lib/libpapi.so')
        except:
            with_papi = False


def init_flop():
    if with_papi:
        ireal_time = c_float()
        iproc_time = c_float()
        iflpops = c_longlong()
        imflops = c_float()
        papi.PAPI_flops(byref(ireal_time), byref(iproc_time), byref(iflpops),
                        byref(imflops))


def get_flop(real_time, proc_time, flpops, mflops):
    if with_papi:
        r = papi.PAPI_flops(byref(real_time), byref(proc_time), byref(flpops),
                            byref(mflops))
    else:
        real_time.value = np.nan
        proc_time.value = np.nan
        flpops.value = -1
        mflops.value = np.nan


class TimeoutException(Exception):
    pass


class RunableProcessing(multiprocessing.Process):
    def __init__(self, func, *args, **kwargs):
        self.queue = multiprocessing.Queue(maxsize=1)
        args = (func,) + args
        multiprocessing.Process.__init__(self, target=self.run_func, args=args, kwargs=kwargs)

    def run_func(self, func, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
            self.queue.put((True, result))
        except Exception as e:
            self.queue.put((False, e))

    def done(self):
        return self.queue.full()

    def result(self):
        return self.queue.get()


def timeout(seconds, force_kill=True):
    if seconds==0:
        def wrapper(function):
            return function
        return wrapper
    else:
        def wrapper(function):
            def inner(*args, **kwargs):
                now = time.time()
                proc = RunableProcessing(function, *args, **kwargs)
                proc.start()
                proc.join(seconds)
                if proc.is_alive():
                    if force_kill:
                        proc.terminate()
                    runtime = int(time.time() - now)
                    raise TimeoutException('timed out after {0} seconds'.format(runtime))
                if not proc.done():
                    proc.terminate()
                assert proc.done()
                success, result = proc.result()
                if success:
                    #return time.time() - now, result
                    return result
                else:
                    #raise time.time() - now, result
                    raise result
            return inner
        return wrapper

@timeout(5)
def dense_matrix_rank(M):
    return matrix_rank(M)



# estimate of condition number and norm from lsmr
# http://www.stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf
#@timeout(20)
def _norm_cond(problem_filename):
    problem = read_fclib_format(problem_filename)[1]
    A = csr_matrix(N.SBMtoSparse(problem.M)[1])
    #print "A=", A
    print "A.shape", A.shape
    print "computev lsmr ..."
    r = lsmr(A, np.ones([A.shape[0], 1]))  # solve Ax = 1
    norm_lsmr=r[5]
    cond_lsmr=r[6]
    print "norm_lsr=", norm_lsmr
    print "cond_lsr=", cond_lsmr
    #print "r=", r
    try:
        print "computev svds(A,1) ..."
        _svd = svds(A,1)[1]
        eps = sys.float_info.epsilon
        tol = _svd.max() * max(A.shape) * eps
    except Exception as e :
        print "-->   svds failed to compute the maximum singular value"
        print "-->" ,  e
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
    try:
        print "Compute rank estimate ..."
        rank_estimate=estimate_rank(A.todense(), tol)
        print "rank_estimate", rank_estimate
    except Exception as e :
        print e
    
    #print "svd dense method", svd(A.todense())[1]

    rank_dense = np.nan
    try:
        print "Compute rank dense ..."
        rank_dense = dense_matrix_rank(A.todense())
    except Exception as e :
        print "--> dense_matrix_rank", e
    print "rank_dense", rank_dense
    
    k=min(rank_estimate,A.shape[0]-1)
    try:
        print "Compute svds(A,k)  ..."
        _svd = svds(A,k)[1]
        #print "_svd",_svd
    # except Warning as w:
    #     print "-->   svds warnig in computing ",k," singular values"
    #     print "-->" ,  w
    except Exception as e :
        print "-->   svds failed to compute ",k," singular values"
        print "-->" ,  e

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
    rank_svd = len(nonzero_sv)
    print "rank_svd", len(nonzero_sv)

    if not math.isnan(rank_dense):
        rank=rank_dense
    else:
        if not math.isnan(rank_svd):
            rank=rank_svd
        else :
            rank=rank_estimate
            
    nonzero_sv = nonzero_sv[0:rank]
    # http://personales.unican.es/beltranc/archivos/CNmatricesF.pdf
    return norm_lsmr, cond_lsmr, max(nonzero_sv), min(nonzero_sv), max(nonzero_sv)/min(nonzero_sv), rank, rank_dense, rank_svd, rank_estimate

norm_cond = Memoize(_norm_cond)

def _read_fclib_format(filename):
    #fc_problem = FCL.fclib_read_local(f)
    #solution = FCL.fclib_read_solution(f)
    #print FCL.fclib_merit_local(fc_problem, FCL.MERIT_1, solution)
    #print fc_problem.W.m
    #print solution.u
    #solution.u = np.zeros(fc_problem.W.m * fc_problem.spacedim)
    #solution.r = np.zeros(fc_problem.W.m * fc_problem.spacedim)
    #print FCL.fclib_merit_local(fc_problem, FCL.MERIT_1, solution)
    #problem = N.frictionContact_fclib_read(f)
    #    print '{0}: M.m={1}, numberOfContacts*3/dim = {2}'.format(f, problem.W.m, problem.numberOfContacts*3/problem.W.m)
    fclib_problem = FCL.fclib_read_local(filename)
    numerics_problem =  N.from_fclib_local(fclib_problem)
    return fclib_problem, numerics_problem

def _numberOfDegreeofFreedom(f):
    with h5py.File(f, 'r') as fclib_file:
        
        try:
            r = 6*fclib_file['fclib_local']['info'].attrs['numberOfInvolvedDS']
        except:
            try:
                r = fclib_file['fclib_local']['info'].attrs['numberOfDegreeOfFreedom'][0]
            except:
                r = np.nan
            #print "r=",r
    return r


# def _numberOfInvolvedDS(f):
#     with h5py.File(f, 'r') as fclib_file:
#         try:
#             r = fclib_file['fclib_local']['info'].attrs['numberOfInvolvedDS']
#         except:
#             r = np.nan
#     return r

def _dimension(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['spacedim'][0]
        except:
            r = np.nan
    return r

def _numberOfDegreeofFreedomContacts(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W']['m'][0]
        except:
            r = np.nan
    #print "r=",r
    return r

def _cond_problem(filename):
    problem = read_fclib_format(filename)[1]
    return float(problem.numberOfContacts * 3) / float(numberOfDegreeofFreedom(filename))

def _cond(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['cond']
        except:
            r = np.nan
    #print "r=",r
    return r

def _rank_dense(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W'].attrs['rank_dense']
        except:
            r = np.nan
    #print "r=",r
    return r

read_fclib_format = Memoize(_read_fclib_format)

numberOfDegreeofFreedom = Memoize(_numberOfDegreeofFreedom)

numberOfDegreeofFreedomContacts = Memoize(_numberOfDegreeofFreedomContacts)

dimension = Memoize(_dimension)

cond_problem = Memoize(_cond_problem)

cond = Memoize(_cond)

rank_dense = Memoize(_rank_dense)



class SolverCallback:
    def __init__(self, h5file, data):
        self._offset = 0
        self._data = data
        self._file = h5file

    def get_step(self, reaction, velocity, error):

        self._reactions = self._data['reactions']
        self._velocities = self._data['velocities']
        self._errors = self._data['errors']
        self._offset += 1
        if output_reactions:
                self._reactions.resize(self._offset, 0)
                self._reactions[self._offset - 1, :] = reaction

        if output_velocities:
                self._velocities.resize(self._offset, 0)
                self._velocities[self._offset - 1, :] = velocity

        if output_errors:
                self._errors.resize(self._offset, 0)
                self._errors[self._offset - 1, :] = error
        print "in get_step"

class Caller():

    def __init__(self):
        pass

    def _solver_call(self, solver, *args):
        return solver(*args)

    def __call__(self, tpl):

        solver, filename = tpl

        if hasattr(solver, 'read_fclib_format'):
            sproblem = solver.read_fclib_format(filename)
            problem = read_fclib_format(filename)[1]
        else:
            problem = read_fclib_format(filename)[1]
            sproblem = problem

        if (output_dat) :
            # write a dat file to create a test for Siconos/Numerics
            #N.frictionContact_display(problem)
            datfilename = filename + '.dat'
            N.frictionContact_printInFilename(problem,datfilename )



        pfilename = os.path.basename(os.path.splitext(filename)[0])

        output_filename = '{0}-{1}.hdf5'.format(solver.name(),
                                                pfilename)
        # considered as tmp file
        try:
                os.remove(output_filename)
        except:
                pass

        try:
            self._internal_call(solver, sproblem, filename, pfilename, output_filename)

        except Exception as e:

            print 'Exception in internal call', e

            try:
               os.remove(output_filename)
            except:
                pass

            with h5py.File(output_filename, 'w') as output:

                digest = hashlib.sha256(open(filename, 'rb').read()).digest()

                create_attrs_in_comp_file(output,precision,utimeout,measure_name)
                
                comp_data=output['data']['comp']
                solver_data = comp_data.create_group(solver.name())
                solver_problem_data = solver_data.create_group(pfilename)
                attrs = solver_problem_data.attrs
                info = 1
                time_s = np.nan
                iter = np.nan
                err = np.nan
                real_time = np.nan
                proc_time = np.nan
                flpops = np.nan
                mflops = np.nan

                attrs.create('filename', filename)
                attrs.create('nc',  numberOfDegreeofFreedomContacts(filename))
                attrs.create('nds', numberOfDegreeofFreedom(filename))
                attrs.create('cond_nc', cond_problem(filename))
                attrs.create('digest', digest)
                attrs.create('info', info)
                attrs.create('iter', iter)
                attrs.create('err', err)
                attrs.create('time', time_s)
                attrs.create('real_time', real_time)
                attrs.create('proc_time', proc_time)
                attrs.create('flpops', flpops)
                attrs.create('mflops', mflops)
                attrs.create('precision', precision)
                attrs.create('timeout', utimeout)

                print(filename, numberOfDegreeofFreedomContacts(filename), numberOfDegreeofFreedom(filename), cond_problem(filename), solver.name(), info, iter, err,
                      time_s, real_time, proc_time,
                      flpops, mflops,
                      precision, utimeout)

                with open('report.txt', "a") as report_file:
                    print   >> report_file , (filename, solver.name(), info, iter, err,
                                              time_s, real_time, proc_time,
                                              flpops, mflops,
                                              precision, utimeout)



    @timeout(utimeout)
    def _internal_call(self, solver, problem, filename, pfilename, output_filename):


        with h5py.File(output_filename, 'w') as output:

            create_attrs_in_comp_file(output,precision,utimeout,measure_name)
                
            comp_data=output['data']['comp']
            
            solver_data = comp_data.create_group(solver.name())

            solver_problem_data = solver_data.create_group(pfilename)
            attrs = solver_problem_data.attrs

            psize = numberOfDegreeofFreedomContacts(filename)

            info = None
            iter = None
            err = None
            time_s = None
            real_time = None
            proc_time = None
            flpops = None
            mflops = None

            digest = hashlib.sha256(open(filename, 'rb').read()).digest()

            if psize is not None:
                solver_problem_data.create_dataset('reactions',
                                                   (0, psize),
                                                   maxshape=(None, psize))

                solver_problem_data.create_dataset('velocities',
                                                   (0, psize),
                                                   maxshape=(None, psize))

                solver_problem_data.create_dataset('errors',
                                                   (0, 1),
                                                   maxshape=(None, 1))

            solver_problem_callback = \
              SolverCallback(output, solver_problem_data)

            # need a function, not an instance method for PyObjectCall...
            def pffff(r, v, e):
                solver_problem_callback.get_step(r, v, e)

            if output_errors or output_velocities or output_reactions:
                solver.SolverOptions().callback = pffff

            # get first guess or set guess to zero
            reactions, velocities = solver.guess(filename)

            _, guess_err = N.fc3d_compute_error(read_fclib_format(filename)[1],
                                                             reactions, velocities, precision, solver.SolverOptions())

#            print "guess error:", guess_err

            try:
                again = True
                info = 0
                iter = 0
                err = np.inf
                real_time = 0.
                proc_time = 0.
                flpops = 0.
                mflops = 0.

                # several call to solver if the precision is not reached
                while again:

                    t0 = time.clock()
                    #t0 = time.process_time()
                    result = solver(problem, reactions, velocities)

                    time_s = time.clock() - t0 # on unix, t is CPU seconds elapsed (floating point)
                    #time_s  = time.process_time() -t0
                    fclib_sol = FCL.fclib_solution()

                    fclib_sol.v = None
                    fclib_sol.r = reactions
                    fclib_sol.u = velocities
                    fclib_sol.l = None

                    #nerr = FCL.fclib_merit_local(read_fclib_format(filename)[0],
                    #                             FCL.MERIT_1, fclib_sol)

#                    _, xerr = N.FrictionContact3D_compute_error(read_fclib_format(filename)[1],
#                                                                reactions, velocities, precision, solver.SolverOptions())


                    i_info, i_iter, i_err, i_real_time, i_proc_time, i_flpops, i_mflops = result

                    info = i_info
                    iter += i_iter
                    err = i_err
                    real_time += i_real_time
                    proc_time += i_proc_time
                    flpops += i_flpops
                    mflops = (mflops + i_mflops)/2.

#                   if info == 0 and xerr >= precision:
#                        solver.SolverOptions().iparam[0]=1
#                        solver.SolverOptions().dparam[0]=solver.SolverOptions().dparam[1]/10
#                        print 'precision not reached : ', xerr, '>', precision
#                        again = False
#                    else:
#                        again = False

                    again = False

                    if info != 0:
                        time_s = np.nan
                        iter = np.nan
                        err = np.nan
                        real_time = np.nan
                        proc_time = np.nan
                        flpops = np.nan
                        mflops = np.nan

            except Exception as exception:
                print exception
                info = 1
                time_s = np.nan
                iter = np.nan
                err = np.nan
                real_time = np.nan
                proc_time = np.nan
                flpops = np.nan
                mflops = np.nan

            attrs.create('filename', filename)
            attrs.create('nc', numberOfDegreeofFreedomContacts(filename))
            attrs.create('nds', numberOfDegreeofFreedom(filename))
            attrs.create('cond_nc', cond_problem(filename))
            attrs.create('digest', digest)
            attrs.create('info', info)
            attrs.create('iter', iter)
            attrs.create('err', err)
            attrs.create('time', time_s)
            attrs.create('real_time', real_time)
            attrs.create('proc_time', proc_time)
            attrs.create('flpops', flpops)
            attrs.create('mflops', mflops)
            attrs.create('precision', precision)
            attrs.create('timeout', utimeout)
            # filename, solver name, revision svn, parameters, nb iter, err
            print(filename, numberOfDegreeofFreedomContacts(filename), numberOfDegreeofFreedom(filename), cond_problem(filename), solver.name(), info, iter, err,
                  time_s, real_time, proc_time,
                  flpops, mflops, precision, utimeout)

                    # if info == 1:
                    #     measure_v = np.inf
                    # else:
                    #     if measure == 'flop':
                    #         measure_v = flpops
                    #     elif measure == 'iter':
                    #         measure_v = iter
                    #     elif measure == 'time':
                    #         measure_v = time_s

                    # measure[solver][ip] = measure_v

                    # min_measure[fileproblem] = min(measure_v,
                    #                               min_measure[fileproblem])
                    # ip += 1

                    # comp_file.flush()


class SiconosSolver():
    _name = None
    _gnuplot_name = None
    _API = None
    _TAG = None
    _iparam_iter = None
    _dparam_err = None
    _SO = None

    def _get(self, tab, index):
        if index is not None:
            return tab[index]
        else:
            return None

    def __init__(self, name=None, gnuplot_name=None, API=None, TAG=None, iparam_iter=None,
                 dparam_err=None, maxiter=maxiter, precision=precision):
        self._name = name
        if (gnuplot_name==None):
            self._gnuplot_name = self._name
        else:
            self._gnuplot_name = gnuplot_name
        self._API = API
        self._TAG = TAG
        self._iparam_iter = iparam_iter
        self._dparam_err = dparam_err
        self._SO = N.SolverOptions(TAG)  # set default solver options
        self._SO.iparam[0] = maxiter
        self._SO.dparam[0] = precision

    def SolverOptions(self):
        return self._SO

    def __call__(self, problem, reactions, velocities):
        real_time = c_float()
        proc_time = c_float()
        flpops = c_longlong()
        mflops = c_float()
        init_flop()
        info = self._API(problem, reactions, velocities, self._SO)

        get_flop(real_time, proc_time, flpops, mflops)
        return (info, self._get(self._SO.iparam, self._iparam_iter),
                self._get(self._SO.dparam, self._dparam_err),
                real_time.value, proc_time.value,
                flpops.value, mflops.value)

    def guess(self, filename):
        problem = read_fclib_format(filename)[1]

        with h5py.File(filename, 'r') as f:
            psize = numberOfDegreeofFreedomContacts(filename)
            if with_guess and 'guesses' in f:
                number_of_guesses = f['guesses']['number_of_guesses'][0]
                velocities = f['guesses']['1']['u'][:]
                reactions = f['guesses']['1']['r'][:]
            else:
                # guess is missing
                reactions = np.zeros(psize)
                velocities = np.zeros(psize)

        return reactions, velocities

    def name(self):
        return self._name
    def gnuplot_name(self):
        return self._gnuplot_name

class BogusSolver(SiconosSolver):
    def read_fclib_format(self, filename):
        return BogusInterface.FCLib.fclib_read_local(filename)

def wrap_bogus_solve(problem, reactions, velocities, SO):
    res = BogusInterface.solve_fclib(problem, reactions, velocities, SO)
    return res

bogusPureNewton = BogusSolver(name="BogusPureNewton", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusPureNewton.SolverOptions().iparam[4]=0


bogusPureEnumerative = BogusSolver(name="BogusPureEnumerative", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusPureEnumerative.SolverOptions().iparam[4]=1


bogusHybrid = BogusSolver(name="BogusHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusHybrid.SolverOptions().iparam[4]=2


bogusRevHybrid = BogusSolver(name="BogusRevHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusRevHybrid.SolverOptions().iparam[4]=3



class SiconosHybridSolver(SiconosSolver):

    def guess(self, filename):
        pfilename = os.path.splitext(filename)[0]
        with h5py.File('comp.hdf5', 'r') as comp_file:
            return extern_guess(pfilename, 'NSGS', 1, comp_file)


class SiconosWrappedSolver(SiconosSolver):
    def __call__(self, problem, reactions, velocities):
        return self._API(problem, reactions, velocities, self._SO)


#
# Some solvers
#
nsn_acSTD = SiconosSolver(name="NSN-AlartCurnier",
                           gnuplot_name="NSN-AC",
                           API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                           TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                           iparam_iter=1,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)

nsn_acSTD.SolverOptions().iparam[10] = 0;
nsn_acSTD.SolverOptions().iparam[11] = 0;
nsn_acSTD.SolverOptions().iparam[12] = maxiterls;
nsn_acSTD.SolverOptions().iparam[13] = with_mumps;
nsn_acSTD.SolverOptions().iparam[3] = 10000000

nsn_acSTD_nls = SiconosSolver(name="NSN-AlartCurnier-NLS",
                               gnuplot_name="NSN-AC-NLS",
                               API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                               TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                               iparam_iter=1,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)

nsn_acSTD_nls.SolverOptions().iparam[10] = 0;
nsn_acSTD_nls.SolverOptions().iparam[11] = -1;
nsn_acSTD_nls.SolverOptions().iparam[12] = 0;
nsn_acSTD_nls.SolverOptions().iparam[13] = with_mumps;
nsn_acSTD_nls.SolverOptions().iparam[3] = 10000000

nsn_acJeanMoreau = SiconosSolver(name="NSN-JeanMoreau",
                                  gnuplot_name="NSN-JM",
                                  API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                  TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                  iparam_iter=1,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

nsn_acJeanMoreau.SolverOptions().iparam[10] = 1;
nsn_acJeanMoreau.SolverOptions().iparam[11] = 0;
nsn_acJeanMoreau.SolverOptions().iparam[12] = maxiterls;
nsn_acJeanMoreau.SolverOptions().iparam[13] = with_mumps;
nsn_acJeanMoreau.SolverOptions().iparam[3] = 10000000

nsn_acJeanMoreau_nls = SiconosSolver(name="NSN-JeanMoreau-NLS",
                                      gnuplot_name="NSN-JM-NLS",
                                      API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                      TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                      iparam_iter=1,
                                      dparam_err=1,
                                      maxiter=maxiter, precision=precision)

nsn_acJeanMoreau_nls.SolverOptions().iparam[10] = 1;
nsn_acJeanMoreau_nls.SolverOptions().iparam[11] = -1;
nsn_acJeanMoreau_nls.SolverOptions().iparam[12] = 0;
nsn_acJeanMoreau_nls.SolverOptions().iparam[13] = with_mumps;
nsn_acJeanMoreau_nls.SolverOptions().iparam[3] = 10000000

nsn_acSTDGenerated = SiconosSolver(name="NSN-AlartCurnier-Generated",
                                    gnuplot_name="NSN-AC-Generated",
                                    API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                    TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                    iparam_iter=1,
                                    dparam_err=1,
                                    maxiter=maxiter, precision=precision)

nsn_acSTDGenerated.SolverOptions().iparam[10] = 2;
nsn_acSTDGenerated.SolverOptions().iparam[11] = 0;
nsn_acSTDGenerated.SolverOptions().iparam[12] = maxiterls;
nsn_acSTDGenerated.SolverOptions().iparam[13] = with_mumps;
nsn_acSTDGenerated.SolverOptions().iparam[3] = 10000000

nsn_acSTDGenerated_nls = SiconosSolver(name="NSN-AlartCurnier-Generated-NLS",
                                        gnuplot_name="NSN-AC-Generated-NLS",
                                        API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                        TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                        iparam_iter=1,
                                        dparam_err=1,
                                        maxiter=maxiter, precision=precision)

nsn_acSTDGenerated_nls.SolverOptions().iparam[10] = 2;
nsn_acSTDGenerated_nls.SolverOptions().iparam[11] = -1;
nsn_acSTDGenerated_nls.SolverOptions().iparam[12] = 0;
nsn_acSTDGenerated_nls.SolverOptions().iparam[13] = with_mumps;
nsn_acSTDGenerated_nls.SolverOptions().iparam[3] = 10000000

nsn_acJeanMoreauGenerated = SiconosSolver(name="NSN-JeanMoreau-Generated",
                                           gnuplot_name="NSN-JM-Generated",
                                           API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                           TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                           iparam_iter=1,
                                           dparam_err=1,
                                           maxiter=maxiter, precision=precision)

nsn_acJeanMoreauGenerated.SolverOptions().iparam[10] = 3;
nsn_acJeanMoreauGenerated.SolverOptions().iparam[11] = 0;
nsn_acJeanMoreauGenerated.SolverOptions().iparam[12] = maxiterls;
nsn_acJeanMoreauGenerated.SolverOptions().iparam[13] = with_mumps;
nsn_acJeanMoreauGenerated.SolverOptions().iparam[3] = 10000000


nsn_acJeanMoreauGenerated_lusol = None
if with_mumps:
    # reference
    nsn_acJeanMoreauGenerated_lusol = SiconosSolver(name="NSN-JeanMoreau-Generated-lusol",
                                                     gnuplot_name="NSN-JM-Generated-LUSOL",
                                                     API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                     TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                     iparam_iter=1,
                                                     dparam_err=1,
                                                     maxiter=maxiter, precision=precision)


    nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[10] = 3;
    nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[11] = 0;
    nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[12] = maxiterls;
    nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[13] = 0;
    nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[3] = 10000000


nsn_acJeanMoreauGenerated_nls = SiconosSolver(name="NSN-JeanMoreau-Generated-NLS",
                                               gnuplot_name="NSN-JM-Generated-NLS",
                                               API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                               TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                               iparam_iter=1,
                                               dparam_err=1,
                                               maxiter=maxiter, precision=precision)

nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[10] = 3;
nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[11] = -1;
nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[12] = maxiterls;
nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[13] = with_mumps;
nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[3] = 10000000

nsn_acJeanMoreauGenerated_nls_lusol=None
if with_mumps:
    nsn_acJeanMoreauGenerated_nls_lusol = SiconosSolver(name="NSN-JeanMoreau-Generated-NLS-lusol",
                                                         gnuplot_name="NSN-JM-Generated-NLS-LUSOL",
                                                         API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                         TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                         iparam_iter=1,
                                                         dparam_err=1,
                                                         maxiter=maxiter, precision=precision)

    nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[10] = 3;
    nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[11] = -1;
    nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[12] = maxiterls;
    nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[13] = 0;
    nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[3] = 10000000

nsn_fb_gp = SiconosSolver(name="NSN-FischerBurmeister-GP",
                           gnuplot_name="NSN-FB-GP",
                           API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                           TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                           iparam_iter=1,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)

nsn_fb_gp.SolverOptions().iparam[3] = 1000000
nsn_fb_gp.SolverOptions().iparam[11] = 0
nsn_fb_gp.SolverOptions().iparam[12] = maxiterls
nsn_fb_gp.SolverOptions().iparam[13] = with_mumps

nsn_fb_gp_lusol = None
if with_mumps:
    nsn_fb_gp_lusol = SiconosSolver(name="NSN-FischerBurmeister-GP-lusol",
                                     gnuplot_name="NSN-FB-GP-LUSOL",
                                     API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                     TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                     iparam_iter=1,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)

    nsn_fb_gp_lusol.SolverOptions().iparam[3] = 1000000
    nsn_fb_gp_lusol.SolverOptions().iparam[11] = 0
    nsn_fb_gp_lusol.SolverOptions().iparam[12] = maxiterls
    nsn_fb_gp_lusol.SolverOptions().iparam[13] = 0

nsn_fb_fblsa = SiconosSolver(name="NSN-FischerBurmeister-FBLSA",
                              gnuplot_name="NSN-FB-FBLSA",
                              API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                              TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                              iparam_iter=1,
                              dparam_err=1,
                              maxiter=maxiter, precision=precision)

nsn_fb_fblsa.SolverOptions().iparam[3] = 1000000
nsn_fb_fblsa.SolverOptions().iparam[11] = 1
nsn_fb_fblsa.SolverOptions().iparam[12] = maxiterls
nsn_fb_fblsa.SolverOptions().iparam[13] = with_mumps

nsn_fb_nls = SiconosSolver(name="NSN-FischerBurmeister-NLS",
                            gnuplot_name="NSN-FB-NLS",
                            API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                            TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                            iparam_iter=1,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)

nsn_fb_nls.SolverOptions().iparam[3] = 1000000
nsn_fb_nls.SolverOptions().iparam[11] = -1
nsn_fb_nls.SolverOptions().iparam[12] = 0
nsn_fb_nls.SolverOptions().iparam[13] = with_mumps

nsn_fb_nls_lusol = None
if with_mumps:
    nsn_fb_nls_lusol = SiconosSolver(name="NSN-FischerBurmeister-NLS-lusol",
                                      gnuplot_name="NSN-FB-NLS-LUSOL",
                                      API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                      TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                      iparam_iter=1,
                                      dparam_err=1,
                                      maxiter=maxiter, precision=precision)

    nsn_fb_nls_lusol.SolverOptions().iparam[3] = 1000000
    nsn_fb_nls_lusol.SolverOptions().iparam[11] = -1
    nsn_fb_nls_lusol.SolverOptions().iparam[12] = 0
    nsn_fb_nls_lusol.SolverOptions().iparam[13] = 0


nsn_nm_gp = SiconosSolver(name="NSN-NaturalMap-GP",
                           gnuplot_name="NSN-NM-GP",
                           API=N.fc3d_nonsmooth_Newton_NaturalMap,
                           TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                           iparam_iter=1,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)

nsn_nm_gp.SolverOptions().iparam[3] = 1000000
nsn_nm_gp.SolverOptions().iparam[11] = 0
nsn_nm_gp.SolverOptions().iparam[12] = maxiterls
nsn_nm_gp.SolverOptions().iparam[13] = with_mumps

nsn_nm_gp_lusol = None
if with_mumps:
    nsn_nm_gp_lusol = SiconosSolver(name="NSN-NaturalMap-GP-lusol",
                                     gnuplot_name="NSN-NM-GP-LUSOL",
                                     API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                     TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                     iparam_iter=1,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)

    nsn_nm_gp_lusol.SolverOptions().iparam[3] = 1000000
    nsn_nm_gp_lusol.SolverOptions().iparam[11] = 0
    nsn_nm_gp_lusol.SolverOptions().iparam[12] = maxiterls
    nsn_nm_gp_lusol.SolverOptions().iparam[13] = 0

nsn_nm_fblsa = SiconosSolver(name="NSN-NaturalMap-FBLSA",
                              gnuplot_name="NSN-NM-FBLSA",
                              API=N.fc3d_nonsmooth_Newton_NaturalMap,
                              TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                              iparam_iter=1,
                              dparam_err=1,
                              maxiter=maxiter, precision=precision)

nsn_nm_fblsa.SolverOptions().iparam[3] = 1000000
nsn_nm_fblsa.SolverOptions().iparam[11] = 1
nsn_nm_fblsa.SolverOptions().iparam[12] = maxiterls
nsn_nm_fblsa.SolverOptions().iparam[13] = with_mumps

nsn_nm_nls = SiconosSolver(name="NSN-NaturalMap-NLS",
                            gnuplot_name="NSN-NM-NLS",
                            API=N.fc3d_nonsmooth_Newton_NaturalMap,
                            TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                            iparam_iter=1,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)

nsn_nm_nls.SolverOptions().iparam[3] = 1000000
nsn_nm_nls.SolverOptions().iparam[11] = -1
nsn_nm_nls.SolverOptions().iparam[12] = 0
nsn_nm_nls.SolverOptions().iparam[13] = with_mumps

nsn_nm_nls_lusol = None
if with_mumps:
    nsn_nm_nls_lusol = SiconosSolver(name="NSN-NaturalMap-NLS-lusol",
                                      gnuplot_name="NSN-NM-NLS-LUSOL",
                                      API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                      TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                      iparam_iter=1,
                                      dparam_err=1,
                                      maxiter=maxiter, precision=precision)

    nsn_nm_nls_lusol.SolverOptions().iparam[3] = 1000000
    nsn_nm_nls_lusol.SolverOptions().iparam[11] = -1
    nsn_nm_nls_lusol.SolverOptions().iparam[12] = 0
    nsn_nm_nls_lusol.SolverOptions().iparam[13] = 0


hnsn_ac = SiconosHybridSolver(name = "HLocalAlartCurnier",
                               API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                               TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                               iparam_iter=1,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)
hnsn_ac.SolverOptions().iparam[3] = 10000000




nsgs = SiconosSolver(name="NSGS-AC",
                     API=N.fc3d_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsgs.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_AC
nsgs.SolverOptions().internalSolvers.iparam[10]=0


nsgs_ac_gp = SiconosSolver(name="NSGS-AC-GP",
                     API=N.fc3d_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsgs_ac_gp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP
nsgs_ac_gp.SolverOptions().internalSolvers.iparam[10]=0

nsgs_jm = SiconosSolver(name="NSGS-JM",
                     API=N.fc3d_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsgs_jm.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_AC
nsgs_jm.SolverOptions().internalSolvers.iparam[10]=1


nsgs_jm_gp = SiconosSolver(name="NSGS-JM-GP",
                     API=N.fc3d_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsgs_jm_gp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP
nsgs_jm_gp.SolverOptions().internalSolvers.iparam[10]=1


snsgs = SiconosSolver(name="NSGS-AC-GP-Shuffled",
                      gnuplot_name="NSGS-AC-GP Shuffled",
                      API=N.fc3d_nsgs,
                      TAG=N.SICONOS_FRICTION_3D_NSGS,
                      iparam_iter=7,
                      dparam_err=1,
                      maxiter=maxiter, precision=precision)
snsgs.SolverOptions().iparam[5] = 1
#print snsgs.SolverOptions().iparam[6]

nsgs_sfull = SiconosSolver(name="NSGS-AC-GP-Shuffled-full",
                           gnuplot_name="NSGS-AC-GP Fully shuffled",
                           API=N.fc3d_nsgs,
                           TAG=N.SICONOS_FRICTION_3D_NSGS,
                           iparam_iter=7,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)
nsgs_sfull.SolverOptions().iparam[5] = 2
#nsgs_sfull.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
#N.printSolverOptions(nsgs_sfull.SolverOptions())


nsgs_pli = SiconosSolver(name="NSGS-PLI",
                         gnuplot_name="NSGS-FP-VI-UPK iter=100",
                         API=N.fc3d_nsgs,
                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                         iparam_iter=7,
                         dparam_err=1,
                         maxiter=maxiter, precision=precision)
nsgs_pli.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
nsgs_pli.SolverOptions().internalSolvers.iparam[0] = 100

nsgs_pli_10 = SiconosSolver(name="NSGS-PLI-10",
                            gnuplot_name="NSGS-FP-VI-UPK iter=10",
                            API=N.fc3d_nsgs,
                            TAG=N.SICONOS_FRICTION_3D_NSGS,
                            iparam_iter=7,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)
nsgs_pli_10.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
nsgs_pli_10.SolverOptions().internalSolvers.iparam[0] = 10


nsgs_p = SiconosSolver(name="NSGS-P",
                       gnuplot_name="NSGS-FP-DS-One",
                       API=N.fc3d_nsgs,
                       TAG=N.SICONOS_FRICTION_3D_NSGS,
                       iparam_iter=7,
                       dparam_err=1,
                       maxiter=maxiter, precision=precision)
nsgs_p.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone

nsgs_pd = SiconosSolver(name="NSGS-PD",
                        gnuplot_name="NSGS-FP-DS-One  D",
                        API=N.fc3d_nsgs,
                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                        iparam_iter=7,
                        dparam_err=1,
                        maxiter=maxiter, precision=precision)
nsgs_pd.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization

nsgs_pr = SiconosSolver(name="NSGS-PR",
                        gnuplot_name="NSGS-FP-DS-One  R",
                        API=N.fc3d_nsgs,
                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                        iparam_iter=7,
                        dparam_err=1,
                        maxiter=maxiter, precision=precision)
nsgs_pr.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization




local_tol_values = [1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16]
#local_tol_values = [1e-2,1e-6,1e-10,1e-16]
nsgs_series=[]
for local_tol in local_tol_values:
    str1 = "{0:1.0e}".format(local_tol).replace("1e","10\^{")+"}"
    nsgs_solver = SiconosSolver(name="NSGS-AC-GP-"+str(local_tol),
                                gnuplot_name="NSGS-AC-GP \$tol\_{local}="+str1+"\$",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=7,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)
    nsgs_solver.SolverOptions().internalSolvers.dparam[0] = local_tol
    nsgs_series.append(nsgs_solver)

for local_tol in local_tol_values:
    str1 = "{0:1.0e}".format(local_tol).replace("1e","10\^{")+"}"
    nsgs_solver = SiconosSolver(name="NSGS-PLI-"+str(local_tol),
                                gnuplot_name="NSGS-FP-VI-UPK \$tol\_{local}="+str1+"\$",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=7,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)
    nsgs_pli.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
    nsgs_solver.SolverOptions().internalSolvers.dparam[0] = local_tol
    nsgs_solver.SolverOptions().internalSolvers.iparam[0] = 100
    nsgs_series.append(nsgs_solver)


snsgs_series=[]
for i in range(10):
    snsgs_solver = SiconosSolver(name="NSGS-AC-Shuffled-"+str(i),
                                 gnuplot_name="NSGS-AC-GP Shuffled "+str(i),
                                 API=N.fc3d_nsgs,
                                 TAG=N.SICONOS_FRICTION_3D_NSGS,
                                 iparam_iter=7,
                                 dparam_err=1,
                                 maxiter=maxiter, precision=precision)
    snsgs_solver.SolverOptions().iparam[5] = 1
    snsgs_solver.SolverOptions().iparam[6] = (1237*i)*(1237*i)
    #print snsgs_solver.SolverOptions().iparam[6]
    snsgs_series.append(snsgs_solver)



 
# only dense
nsgsv = SiconosSolver(name="NSGS-Velocity",
                      API=N.fc3d_nsgs_velocity,
                      TAG=N.SICONOS_FRICTION_3D_NSGSV,
                      iparam_iter=7,
                      dparam_err=1,
                      maxiter=maxiter, precision=precision)


omega=1.5
psor = SiconosSolver(name="PSOR-AC",
                     gnuplot_name="PSOR-AC",
                     API=N.fc3d_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
psor.SolverOptions().iparam[4] = 1
psor.SolverOptions().dparam[8] = omega

omega_values = [0.5, 0.8, 1.0, 1.1, 1.3, 1.5, 1.8]
psor_series=[]
for omega in omega_values:
    psor_solver = SiconosSolver(name="PSOR-AC-"+str(omega),
                                gnuplot_name="PSOR-AC \$\\\omega="+str(omega)+"\$",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=7,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)
    psor_solver.SolverOptions().iparam[4] = 1
    psor_solver.SolverOptions().dparam[8] = omega
    psor_series.append(psor_solver)

TrescaFixedPoint = SiconosSolver(name="TrescaFixedPoint-NSGS-PLI",
                                 API=N.fc3d_TrescaFixedPoint,
                                 TAG=N.SICONOS_FRICTION_3D_TFP,
                                 iparam_iter=7,
                                 dparam_err=1,
                                 maxiter=maxiter, precision=precision)

ACLMFixedPoint = SiconosSolver(name="ACLMFixedPoint-SOCLCP-NSGS-PLI",
                               API=N.fc3d_ACLMFixedPoint,
                               TAG=N.SICONOS_FRICTION_3D_ACLMFP,
                               iparam_iter=7,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)

SOCLCP = SiconosSolver(name="SOCLCP-NSGS-PLI",
                       API=N.fc3d_SOCLCP,
                       TAG=N.SICONOS_FRICTION_3D_SOCLCP,
                       iparam_iter=7,
                       dparam_err=1,
                       maxiter=maxiter, precision=precision)

DeSaxceFixedPoint = SiconosSolver(name="FixedPoint-DeSaxce",
                                  API=N.fc3d_DeSaxceFixedPoint,
                                  TAG=N.SICONOS_FRICTION_3D_DSFP,
                                  iparam_iter=7,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

ExtraGrad = SiconosSolver(name="ExtraGradient",
                          API=N.fc3d_ExtraGradient,
                          TAG=N.SICONOS_FRICTION_3D_EG,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

FixedPointProjection = SiconosSolver(name="FixedPoint-Projection",
                                     API=N.fc3d_fixedPointProjection,
                                     TAG=N.SICONOS_FRICTION_3D_FPP,
                                     iparam_iter=7,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)

VIExtraGrad = SiconosSolver(name="ExtraGradient-VI",
                            API=N.fc3d_VI_ExtraGradient,
                            TAG=N.SICONOS_FRICTION_3D_VI_EG,
                            iparam_iter=7,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)
# VIExtraGrad.SolverOptions().dparam[4]=0.6
# VIExtraGrad.SolverOptions().dparam[5]=1/0.7
# VIExtraGrad.SolverOptions().dparam[6]=0.9
# VIExtraGrad.SolverOptions().dparam[7]=0.3

iparam1_values = [0,1]

iparam2_values = [0,1]

iparam3_values = [0]

VIExtraGrad_series=[]
for i1 in iparam1_values:
    for i2 in iparam2_values:
        for i3 in iparam3_values:
            if i1 == 0 :
                g_name="EG-VI-UPK"
            elif i1 == 1 :
                g_name="EG-VI-UPTS"
            elif i1 == 2:
                g_name="EG-VI-UPHS"

            if i2 == 0:
                g_name = g_name + " False"
            elif i2 == 1:
                g_name = g_name + " True"
                
            if i3 == 0:
                g_name = g_name 
            elif i3 == 1:
                g_name = g_name + " min" 

            VIExtraGrad_solver= SiconosSolver(name="ExtraGrad-VI-"+str(i1)+str(i2)+str(i3),
                                                         gnuplot_name=g_name,
                                                         API=N.fc3d_VI_ExtraGradient,
                                                         TAG=N.SICONOS_FRICTION_3D_VI_EG,
                                                         iparam_iter=7,
                                                         dparam_err=1,
                                                         maxiter=maxiter, precision=precision)
            VIExtraGrad_solver.SolverOptions().iparam[1] = i1
            VIExtraGrad_solver.SolverOptions().iparam[2] = i2
            VIExtraGrad_solver.SolverOptions().iparam[3] = i3
            VIExtraGrad_series.append(VIExtraGrad_solver)

            
VIFixedPointProjection = SiconosSolver(name="FixedPoint-VI",
                                       API=N.fc3d_VI_FixedPointProjection,
                                       TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                       iparam_iter=7,
                                       dparam_err=1,
                                       maxiter=maxiter, precision=precision)

iparam1_values = [0,1,2]
iparam1_values = [0,1]

iparam2_values = [0]
iparam2_values = [0,1 ]


iparam3_values = [0,1]
iparam3_values = [0]
VIFixedPointProjection_series=[]
for i1 in iparam1_values:
    for i2 in iparam2_values:
        for i3 in iparam3_values:
            if i1 == 0 :
                g_name="FP-VI-UPK"
            elif i1 == 1 :
                g_name="FP-VI-UPTS"
            elif i1 == 2:
                g_name="FP-VI-UPHS"

            if i2 == 0:
                g_name = g_name + " False"
            elif i2 == 1:
                g_name = g_name + " True"
            elif i2 == 2:
                g_name = g_name + " semi False"
                
            if i3 == 0:
                g_name = g_name 
            elif i3 == 1:
                g_name = g_name + " min" 

            VIFixedPointProjection_solver= SiconosSolver(name="FixedPoint-VI-"+str(i1)+str(i2)+str(i3),
                                                         gnuplot_name=g_name,
                                                         API=N.fc3d_VI_FixedPointProjection,
                                                         TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                                         iparam_iter=7,
                                                         dparam_err=1,
                                                         maxiter=maxiter, precision=precision)
            VIFixedPointProjection_solver.SolverOptions().iparam[1] = i1
            VIFixedPointProjection_solver.SolverOptions().iparam[2] = i2
            VIFixedPointProjection_solver.SolverOptions().iparam[3] = i3
            VIFixedPointProjection_series.append(VIFixedPointProjection_solver)

Prox = SiconosSolver(name="PROX-NSN-AC",
                     gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Prox.SolverOptions().internalSolvers.iparam[3] = 100


Proxfixed = SiconosSolver(name="PROX-NSN-AC-fixed",
                     gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu=1, \\\sigma=5.0\$ fixed",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Proxfixed.SolverOptions().dparam[3] = -1e4
Proxfixed.SolverOptions().internalSolvers.iparam[3] = 100


regul_value=1e03
str1 = "{0:1.0e}".format(regul_value).replace("1e","10\^{")+"}"
str2 = "{0:1.0e}".format(regul_value)
Regul_variable = SiconosSolver(name="PROX-NSN-AC-regulVar-"+str2,
                     gnuplot_name="PPA-NSN-AC-GP   \$ \\\mu=1, \\\sigma=5.0\$ regulVar "+str1,
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Regul_variable.SolverOptions().dparam[3] = regul_value
Regul_variable.SolverOptions().iparam[9] = 1
Regul_variable.SolverOptions().internalSolvers.iparam[3] = 100


regul_values= [ 1e4, 1e6, 1e8, 1e10]
regul_series =[]
for rr in regul_values:
    str1 = "{0:1.0e}".format(rr).replace("1e","10\^{")+"}"
    str2 = "{0:1.0e}".format(rr)
    regul_solver  = SiconosSolver(name="PROX-NSN-AC-regul-"+str2,
                                  gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu="+str1+"\$",
                                  API=N.fc3d_proximal,
                                  TAG=N.SICONOS_FRICTION_3D_PROX,
                                  iparam_iter=7,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)
    regul_solver.SolverOptions().dparam[3] = - rr
    regul_solver.SolverOptions().iparam[9] = 1
    regul_solver.SolverOptions().internalSolvers.iparam[3] = 100
    regul_series.append(regul_solver)

Prox_nls = SiconosSolver(name="PROX-NSN-AC-NLS",
                     gnuplot_name="PPA-NSN-AC  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Prox_nls.SolverOptions().internalSolvers.iparam[11] = -1

ProxFB = SiconosSolver(name="PROX-NSN-FB-GP",
                     gnuplot_name="PPA-NSN-FB-GP  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
ProxFB.SolverOptions().internalSolvers.iparam[3] = 100
ProxFB.SolverOptions().dparam[4]=5.0 # sigma
ProxFB.SolverOptions().dparam[5]=1.0 # nu
nsn_fb_gp_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
nsn_fb_gp_inprox.iparam[3] = 1000000
nsn_fb_gp_inprox.iparam[11] = 0
nsn_fb_gp_inprox.iparam[12] = 6 #maxiterls
nsn_fb_gp_inprox.iparam[13] = with_mumps

ProxFB.SolverOptions().internalSolvers = nsn_fb_gp_inprox
ProxFB.SolverOptions().internalSolvers.iparam[3] = 100



ProxFB_fblsa = SiconosSolver(name="PROX-NSN-FB-FBLSA",
                     gnuplot_name="PPA-NSN-FB-FBLSA  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsn_fb_fblsa_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
nsn_fb_fblsa_inprox.iparam[3] = 1000000
nsn_fb_fblsa_inprox.iparam[11] = 1
nsn_fb_fblsa_inprox.iparam[12] = 6 #maxiterls
nsn_fb_fblsa_inprox.iparam[13] = with_mumps
ProxFB_fblsa.SolverOptions().internalSolvers = nsn_fb_fblsa_inprox
ProxFB_fblsa.SolverOptions().internalSolvers.iparam[3] = 1000000

ProxFB_nls = SiconosSolver(name="PROX-NSN-FB-NLS",
                     gnuplot_name="PPA-NSN-FB-NLS  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsn_fb_nls_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
nsn_fb_nls_inprox.iparam[3] = 1000000
nsn_fb_nls_inprox.iparam[11] = -1
nsn_fb_nls_inprox.iparam[13] = with_mumps
ProxFB_nls.SolverOptions().internalSolvers = nsn_fb_nls_inprox
ProxFB_nls.SolverOptions().internalSolvers.iparam[3] = 1000000

ProxNSGS = SiconosSolver(name="PROX-NSGS-NSN-AC",
                     gnuplot_name="PROX-NSGS-NSN-AC \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.fc3d_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
local_nsgs_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSGS)
ProxNSGS.SolverOptions().internalSolvers = local_nsgs_inprox
ProxNSGS.SolverOptions().internalSolvers.iparam[0] = 10000







sigmavalues= [0.5, 1.0, 4.0, 5.0, 50, 100.0, 1000.0 ]
muvalues= [0.5, 1.0, 2.0]
prox_series =[]
for mu in muvalues:
    for sigma in sigmavalues:
        prox_solver  = SiconosSolver(name="PROX-NSN-AC-nu"+str(mu)+"-sigma"+str(sigma),
                                     gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu="+str(mu)+", \\\sigma="+str(sigma)+"\$",
                                     API=N.fc3d_proximal,
                                     TAG=N.SICONOS_FRICTION_3D_PROX,
                                     iparam_iter=7,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)
        prox_solver.SolverOptions().internalSolvers.iparam[3] = 1000000
        prox_solver.SolverOptions().dparam[4]=sigma # sigma
        prox_solver.SolverOptions().dparam[5]=mu # nu
        prox_series.append(prox_solver)

nsn_ac_wrapped = SiconosSolver(name="NSN-AlartCurnier-Wrapped",
                                API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                iparam_iter=1,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)

def fc3d_nsn_ac_r(problem, reactions, velocities, _SO):
    SO = N.SolverOptions(N.SICONOS_FRICTION_3D_VI_FPP)
    SO.iparam[3] = 1000
    N.fc3d_VI_FixedPointProjection(problem, reactions, velocities, SO)
    #    print '->',SO.dparam[3]
    nsn_ac_wrapped.SolverOptions().dparam[3] = SO.dparam[3]
    return nsn_ac_wrapped(problem, reactions, velocities)

# flop measure only on nsn_ac
nsn_acr = SiconosWrappedSolver(name="NSN-AlartCurnier-R",
                                API=fc3d_nsn_ac_r,
                                TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                iparam_iter=1,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)

#
quartic = SiconosSolver(name="NSGS-Quartic",
                        API=N.fc3d_nsgs,
                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                        iparam_iter=7,
                        dparam_err=1, maxiter=maxiter, precision=precision)

quartic3x3 = N.SolverOptions(N.SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU)

quartic.SolverOptions().internalSolvers = quartic3x3

# 1 contact
#AlartCurnierNewton = SiconosSolver(name="AlartCurnierNewton",
#                                   API=fc3d_AlartCurnierNewton,
#                                   iparam_iter=1,
#                                   iparam_iter=1)

# rho estimation needed
Prox._SO.dparam[3] = 1000

HyperplaneProjection = SiconosSolver(name="HyperplaneProjection",
                                     API=N.fc3d_HyperplaneProjection,
                                     TAG=N.SICONOS_FRICTION_3D_HP,
                                     iparam_iter=7,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)



# check for correct flop
#def pipo(*args,**kwargs):
#    a=1.
#    b=3.
#    c = a+b
#    return 0
#
#Pipo = SiconosSolver(name="pipo",
#                     API=pipo,
#                     TAG=SICONOS_FRICTION_3D_HP,
#                     iparam_iter=1,
#                     dparam_err=1)
# http://code.google.com/p/mpi4py/
#from mpi4py import MPI
#rank = MPI.COMM_WORLD.rank
#fc3d_sparseGlobalAlartCurnierInit(nsn_ac.SolverOptions())

#all_solvers = [nsgs, snsgs, TrescaFixedPoint, nsn_ac, Prox, DeSaxceFixedPoint,
#               FixedPointProjection, VIFixedPointProjection, ExtraGrad, VIExtraGrad]
#all_solvers = [nsgs, snsgs, quartic, TrescaFixedPoint, ACLMFixedPoint, DeSaxceFixedPoint, VIFixedPointProjection, VIFixedPointProjection1, VIFixedPointProjection2, VIFixedPointProjection3, VIExtraGrad, SOCLCP, Prox, Prox2, Prox3, Prox4, Prox5, nsn_acSTD, nsn_acSTDGenerated,  nsn_acr, nsn_acJeanMoreau, nsn_acJeanMoreauGenerated, nsn_fb_gp, nsn_fb_fblsa]

nsgs_solvers = [nsgs, nsgs_ac_gp, nsgs_jm, nsgs_jm_gp, nsgs_sfull, snsgs, nsgs_pli, nsgs_pli_10, nsgs_p,nsgs_pd,nsgs_pr , quartic]
# remove very nasty solver
#nsgs_solvers.remove(nsgs_p)
nsgs_solvers.remove(nsgs_pd)
nsgs_solvers.remove(nsgs_pr)
nsgs_solvers.remove(quartic)

nsn_solvers =  [nsn_acSTD, nsn_acSTD_nls, nsn_acSTDGenerated, nsn_acSTDGenerated_nls,  nsn_acr, nsn_acJeanMoreau, nsn_acJeanMoreau_nls, nsn_acJeanMoreauGenerated, nsn_acJeanMoreauGenerated_lusol,
                nsn_acJeanMoreauGenerated_nls, nsn_acJeanMoreauGenerated_nls_lusol,
                nsn_fb_gp, nsn_fb_gp_lusol, nsn_fb_nls, nsn_fb_nls_lusol,
                nsn_nm_gp, nsn_nm_gp_lusol, nsn_nm_nls, nsn_nm_nls_lusol]

all_solvers = list(nsgs_solvers)
all_solvers.extend(nsn_solvers)
all_solvers.extend( [ psor,
                      TrescaFixedPoint, DeSaxceFixedPoint,
                      VIFixedPointProjection, VIExtraGrad,
                      SOCLCP,
                      Prox,  Prox_nls, ProxFB,  ProxFB_nls, ProxNSGS, Proxfixed, Regul_variable, regul_series[0],
                      ACLMFixedPoint])

all_solver_unstable = [ProxFB_fblsa]
all_solvers.extend(all_solver_unstable)

all_solvers = filter(lambda s : s is not None, all_solvers)

###
# specific studies of solvers.

all_solvers.extend(VIFixedPointProjection_series)
all_solvers.extend(VIExtraGrad_series)
all_solvers.extend(psor_series)
all_solvers.extend(prox_series)
all_solvers.extend(regul_series)
all_solvers.extend(nsgs_series)


if (os.path.isfile(os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'))):
    execfile(os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'))
    #print "execfile(",os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'), ")"
if (os.path.isfile('adhoc_solverlist.py')):
    execfile('adhoc_solverlist.py')
    #print "execfile(adhoc_solverlist.py)"

    

solvers=[]
if user_solvers != []:
    #print "user_solvers", user_solvers
    solvers.extend( filter(lambda s: any(us in s._name for us in user_solvers), all_solvers))

    if solvers == []:
        raise RuntimeError ("Cannot find any matching solver")

elif user_solvers_exact != []:
    #print "user_solvers_exact", user_solvers_exact
    solvers.extend(filter(lambda s: any(us ==  s._name  for us in user_solvers_exact), all_solvers))

    if solvers == []:
        raise RuntimeError("Cannot find any solvers in specified list")
    
else:
    solvers= all_solvers


#solvers = [ProxFB]
    
def is_fclib_file(filename):
    r = False
    try:
        with h5py.File(filename, 'r') as f:
            r = 'fclib_local' in f or 'fclib_global' in f
    except:
        pass

    return r

def read_numerics_format(f):
    return N.frictionContactProblemFromFile(f)




from numpy import eye, array
keeper = []

NC = 1

M = eye(3*NC)


q = array([-1., 1., 3.])

mu = array([0.1]);

z = array([0.,0.,0.])

reactions = array([0.,0.,0.])

velocities = array([0.,0.,0.])

measure = dict()
solver_r = dict()

if not os.path.exists('problems.txt'):
    with open('problems.txt', 'w') as problems_txt:
        for f in glob('*.hdf5'):
            problems_txt.write('{0}\n'.format(f))

if user_filenames == []:
    if file_filter == None:
        all_filenames = list_from_file('problems.txt')
    else:
        all_filenames = filter(lambda f: any(uf in f for uf in file_filter), list_from_file('problems.txt'))
   
else:
    
        all_filenames = user_filenames
   
#all_filenames=['BoxesStack1-i9841-33.hdf5']
#ask_collect = False

_problem_filenames = filter(is_fclib_file,
                            all_filenames)
    
__problem_filenames = subsample_problems(_problem_filenames,
                                         random_sample_proba,
                                         max_problems, None, overwrite = (not display and not ask_compute and not ask_collect))


problem_filenames = subsample_problems(__problem_filenames,
                                       None,
                                       None, cond_nc)

n_problems = len(problem_filenames)

#problems = [read_fclib_format(f) for f in problem_filenames]



min_measure = dict()

for fileproblem in problem_filenames:
    min_measure[fileproblem] = np.inf

if clean:
    h5mode = 'w'
else:
    h5mode = 'a'

caller = Caller()


#pool = MyPool(processes=8)


def create_attrs_precision_in_comp_file(comp_file,precision_val):
    data = comp_file.get('data')
    if data == None :
        data = comp_file.create_group('data')
    comp_data = data.get('comp')
    if comp_data == None: 
        comp_data = data.create_group('comp')
    comp_data.attrs.create('precision',precision_val)

def create_attrs_timeout_in_comp_file(comp_file,utimeout_val):
    data = comp_file.get('data')
    if data == None :
        data = comp_file.create_group('data')
    comp_data = data.get('comp')
    if comp_data == None: 
        comp_data = data.create_group('comp')
    comp_data.attrs.create('timeout',utimeout_val)



def create_attrs_in_comp_file(comp_file,precision_val,utimeout_val,measure_name_val):
    create_attrs_precision_in_comp_file(comp_file,precision_val)
    create_attrs_timeout_in_comp_file(comp_file,utimeout_val)

def collect(tpl):

    solver, filename = tpl

    pfilename = os.path.basename(os.path.splitext(filename)[0])
    results_filename = '{0}-{1}.hdf5'.format(solver.name(),pfilename)
    #print "file=", results_filename
    if os.path.exists('comp.hdf5'):
        with h5py.File('comp.hdf5', 'r') as comp_file:
            comp_precision=comp_file['data']['comp'].attrs.get('precision')
            comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            comp_measure_name=comp_file['data']['comp'].attrs.get('mesaure_name')
            #print "comp_precision",comp_precision
            if comp_precision == None :
                raise RuntimeError ("Warning. precision information is missing in existing comp.hdf5 file (old version)\n      you must add it with --add-precision-in-comp-file=<val> ")
            if comp_utimeout == None :
                raise RuntimeError ("Warning. timeout information is missing in existing comp.hdf5 file (old version)\n      you must add it with --add-timeout-in-comp-file=<val> ")
    else:
        with h5py.File('comp.hdf5', 'w') as comp_file:
            create_attrs_in_comp_file(comp_file,precision,utimeout,measure_name)
            comp_precision=comp_file['data']['comp'].attrs.get('precision')
            comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            comp_measure_name=comp_file['data']['comp'].attrs.get('mesaure_name')
            
    if os.path.exists(results_filename) and not os.stat(results_filename).st_size == 0:
        try:
            if os.path.exists('comp.hdf5'):                
                with h5py.File( results_filename, 'r+') as result_file:
                    result_precision=result_file['data']['comp'].attrs.get('precision')
                    result_utimeout=result_file['data']['comp'].attrs.get('timeout')
                    result_measure_name=result_file['data']['comp'].attrs.get('measure_name')                                    
                    if comp_precision != result_precision:
                        raise RuntimeError ("Precision of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_precision,result_precision))
                    if comp_utimeout != result_utimeout:
                        raise RuntimeError ("Timeout of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_utimeout,result_utimeout))
                           
            check_call(['h5copy','-p','-i', results_filename,
                        '-ocomp.hdf5','-s/data/comp/{0}/{1}'.format(solver.name(),pfilename),
                        '-d/data/comp/{0}/{1}'.format(solver.name(),pfilename)])
            if not keep_files:
                os.remove('{0}-{1}.hdf5'.format(solver.name(),pfilename))
        except Exception as e:
            print e


class Results():
    def __init__(self, result_file):
        self._result_file = result_file

    def __call__(self, tpl):
        solver = tpl[0]
        problem_filename = os.path.splitext(tpl[1])[0]
        try:
            r = self._result_file['data']['comp'][solver.name()][problem_filename]
            # if abs(r.attrs.get('precision') -  precision) >= 1e-16 :
            #     raise RuntimeError()
            # if abs(r.attrs.get('timeout') -  utimeout) >= 1e-16 :
            #    raise RuntimeError()
            print "Already in comp file : ", (r.attrs['filename'], cond_problem(r.attrs['filename']), solver.name(), r.attrs['info'],
                                              r.attrs['iter'], r.attrs['err'], r.attrs['time'], r.attrs['real_time'], r.attrs['proc_time'],
                                              r.attrs['flpops'], r.attrs['mflops'],r.attrs.get('precision'),r.attrs.get('timeout'))
            return False
        except:
            return True

if __name__ == '__main__':

    if compute:
        print "Tasks will be run for solvers :", [ s._name for s in solvers]
        print " on files ",problem_filenames
        

        
        all_tasks = [t for t in product(solvers, problem_filenames)]

        if os.path.exists('comp.hdf5'):
            with h5py.File('comp.hdf5', 'r') as comp_file:
                tasks = filter(Results(comp_file), all_tasks)

        else:
            tasks = all_tasks

        if ask_compute:
            r = map(caller, tasks)

        if ask_collect:
            map(collect, tasks)
            
    if list_contents:
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            for item in comp_data.attrs.keys():
                print "comp_data attrs: ", item + ":", comp_data.attrs[item]
            print "Solvers :"
            for solvername in comp_data:
                print "  ",solvername
                for filename in comp_data[solvername]:
                    list_keys= comp_data[solvername][filename].attrs.keys()
                    if u'digest' in list_keys:
                        list_keys.remove(u'digest')
                    print "  ",solvername,   [comp_data[solvername][filename].attrs[item] for item in list_keys]

                    
            
    if display:
        print "Tasks will be run for solvers :", [ s._name for s in solvers]
        filename=None
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            # 1 n_problems
            n_problems = 0

            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    if file_filter == None:
                        all_filenames = comp_data[solver_name]
                    else:
                        all_filenames = filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name])

                    filenames = subsample_problems(all_filenames,
                                                   random_sample_proba,
                                                   max_problems, cond_nc)

                    n_problems = max(n_problems, len(filenames))

            # 2 measures & min_measure

            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    if file_filter == None:
                        all_filenames = comp_data[solver_name]
                    else:
                        all_filenames = filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name])

                    filenames = subsample_problems(all_filenames,
                                                   random_sample_proba,
                                                   max_problems, cond_nc)

                    
                    assert len(filenames) <= n_problems

                    measure[solver_name] = np.inf * np.ones(n_problems)
                    solver_r[solver_name] = np.inf * np.ones(n_problems)

                    ip = 0

                    for filename in filenames:
                        if filename not in min_measure:
                            min_measure[filename] = np.inf
                        try:
                            pfilename = os.path.splitext(filename)[0]
                            if comp_data[solver_name][pfilename].attrs['info'] == 0:
                                measure[solver_name][ip] =  comp_data[solver_name][pfilename].attrs[measure_name]
                                min_measure[filename] = min(min_measure[filename], measure[solver_name][ip])
                            else:
                                measure[solver_name][ip] = np.inf
                        except:
                            measure[solver_name][ip] = np.nan
                        ip += 1

            # 3 solver_r
            #            for solver_name in comp_data:
            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    if file_filter == None:
                        all_filenames = comp_data[solver_name]
                    else:
                        all_filenames = filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name])

                    filenames = subsample_problems(all_filenames,
                                                   random_sample_proba,
                                                   max_problems, cond_nc)

                    ip = 0
                    for filename in filenames:
                        pfilename = os.path.splitext(filename)[0]
                        try:
                            if comp_data[solver_name][pfilename].attrs['info'] == 0:
                                solver_r[solver_name][ip] = measure[solver_name][ip] / \
                                  min_measure[filename]

                            else:
                                solver_r[solver_name][ip] = np.inf
                        except:
                            solver_r[solver_name][ip] = np.inf
                        ip += 1

            # 4 rhos
            rhos = dict()
            #for solver_name in comp_data:
            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    assert min(solver_r[solver_name]) >= 1
                    rhos[solver_name] = np.empty(len(domain))
                    for itau in range(0, len(domain)):
                        rhos[solver_name][itau] = float(len(np.where( solver_r[solver_name] <= domain[itau] )[0])) / float(n_problems)


            if (gnuplot_profile and (filename != None)) :
                def write_report(r, filename):
                    with open(filename, "w") as input_file:
                        for k, v in r.items():
                            line = '{}, {}'.format(k, v)
                            print >> input_file, line

                out_data=np.empty([len(domain),len(comp_data)+1])
                write_report(rhos,'rhos.txt')
                write_report(solver_r,'solver_r.txt')
                def long_substr(data):
                    substr = ''
                    if len(data) > 1 and len(data[0]) > 0:
                        for i in range(len(data[0])):
                            for j in range(len(data[0])-i+1):
                                if j > len(substr) and is_substr(data[0][i:i+j], data):
                                    substr = data[0][i:i+j]
                    return substr

                def is_substr(find, data):
                    if len(data) < 1 and len(find) < 1:
                        return False
                    for i in range(len(data)):
                        if find not in data[i]:
                            return False
                    return True


                with open('profile.gp','w') as gp:
                    #  all_rhos = [ domain ] + [ rhos[solver_name] for solver_name in comp_data ]
                    all_rhos = [ domain ] + [ rhos[solver.name()] for solver in filter(lambda s: s._name in comp_data, solvers) ]
                    np.savetxt('profile.dat', np.matrix(all_rhos).transpose())
                    gp.write('resultfile = "profile.dat"\n')
                    test_name = long_substr(filenames).partition('-')[0]
                    if test_name.endswith('_'):
                        test_name  = test_name[:-1]
                    gp.write('basename="profile-{0}"\n'.format(test_name))
                    #print filename.partition('-')[0]
                    print test_name
                    gp.write('\n')
                    gp.write('term_choice_tikz=1\n')
                    gp.write('if (term_choice_tikz == 1) \\\n')
                    gp.write('set term tikz standalone monochrome  size 5in,3in font \'\\scriptsize\\sf\';  \\\n')
                    gp.write('extension = \'.tex\'; \\\n')
                    gp.write('extension_legend = \'_legend.tex\'; \\\n')
                    gp.write('set output basename.extension; \\\n')
                    gp.write('print "output = ", basename.extension; \\\n')

                    gp.write('else \\\n')
                    gp.write('set term aqua;\\\n')
                    gp.write('\n')
                    gp.write('set xrange [{0}:{1}]\n'.format(domain[0]-0.01, domain[len(domain)-1]))
                    gp.write('set yrange [-0.01:1.01]\n')
                    gp.write('set ylabel \'$\\rho(\\tau)$ \' \n')
                    maxrows=len(solvers)/2+1
                    gp.write('set key below right vertical maxrows {0}\n'.format(maxrows))
                   
                    if logscale:
                        gp.write('set logscale x\n')
                        gp.write('set xlabel \'$\\tau$ ({0}) (logscale)\' \n'.format(measure_name))
                    else:
                        gp.write('set xlabel \'$\\tau$ ({0})\' \n'.format(measure_name))

                    #gp.write('set title \'{0}\'\n'.format(filename.partition('-')[0]));
                    gp.write('plot ')
                    if gnuplot_separate_keys:
                        if (gnuplot_with_color):
                            gp.write(
                                ','.join(['resultfile using 1:{0} notitle w l  dashtype {1} linecolor {2}'.format(index + 2,index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
                            
                            gp.write('\n set output basename.extension_legend; \n')
                            gp.write('print "output = ", basename.extension_legend; \n \n')
                            gp.write('unset border; \n \n')
                            gp.write('unset tics; \n \n')
                            gp.write('unset xlabel; \n \n')
                            gp.write('unset ylabel; \n \n')
                            gp.write('set term tikz standalone monochrome  size 5in,1.5in font \'\\scriptsize\\sf\';  \\\n')
 
                            gp.write('set key right inside vertical maxrows {0}\n'.format(maxrows))

                            gp.write('\n plot [0:1] [0:1]')
                            gp.write(
                                ','.join([' NaN t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))

                            
                            
                        else:
                            gp.write(
                                ','.join(['resultfile using 1:{0} t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
                    else:
                        if (gnuplot_with_color):
                            gp.write(
                                ','.join(['resultfile using 1:{0} t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
                        else:
                            gp.write(
                                ','.join(['resultfile using 1:{0} t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
                        
                # all_rhos = [ rhos[solver_name] for solver_name in comp_data ]
                # g.plot(*all_rhos)

            if (gnuplot_profile and (filename == None)) :
                print "Warning: no problem corresponding to the required solver"
                if (os.path.isfile('profile.gp')):
                    os.remove('profile.gp')
                
            if not no_matplot:
                # 5 plot
                from matplotlib.pyplot import subplot, title, plot, grid, show, legend, figure, xlim, ylim, xscale

                #for solver_name in comp_data:
                for solver in solvers:
                    solver_name=solver.name()
                    if logscale:
                        xscale('log')

                    if solver_name in comp_data :
                        plot(domain, rhos[solver_name], label=solver_name)
                        ylim(0, 1.0001)
                        xlim(domain[0], domain[-1])
                        legend(loc=4)
                    grid()

                
    if display_convergence:
        from matplotlib.pyplot import subplot, title, plot, grid, show, legend, figure
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            for solver_name in comp_data:

                if user_filenames == []:
                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
                else:
                    filenames = user_filenames

                for filename in filenames:

                    try:
                        pfilename = os.path.splitext(filename)[0]
                        solver_problem_data = comp_data[solver_name][pfilename]
                        figure()

                        plot(np.arange(len(solver_problem_data['errors'][:])),
                             solver_problem_data['errors'], label='{0} - {1}'.format(solver_name, filename))
                        legend()
                        grid()
                    except:
                        pass
    if compute_cond_rank:
        print "Tasks will be run for", problem_filenames
        for problem_filename in problem_filenames:
            print "compute for", problem_filename,"...."
            with h5py.File(problem_filename, 'r+') as fclib_file:
                no_rank_info=True
                if (fclib_file['fclib_local']['W'].attrs.get('rank') == None) :
                    print "Rank info already not  in", problem_filename
                else:
                    print "Rank info already in", problem_filename
                    no_rank_info=False
            if no_rank_info:
                try:
                    [norm_lsmr, cond_lsmr, max_nz_sv, min_nz_sv, cond, rank, rank_dense, rank_svd, rank_estimate] = norm_cond(problem_filename)
                    print ( problem_filename, norm_lsmr, cond_lsmr, max_nz_sv, min_nz_sv, cond,  rank_dense, rank_svd, rank_estimate)
                    with h5py.File(problem_filename, 'r+') as fclib_file:
                        fclib_file['fclib_local']['W'].attrs.create('rank', rank)
                        fclib_file['fclib_local']['W'].attrs.create('rank_dense', rank_dense)
                        fclib_file['fclib_local']['W'].attrs.create('rank_svd', rank_svd)
                        fclib_file['fclib_local']['W'].attrs.create('rank_estimate', rank_estimate)
                        fclib_file['fclib_local']['W'].attrs.create('cond', cond)
                        fclib_file['fclib_local']['W'].attrs.create('max_nz_sv', max_nz_sv)
                        fclib_file['fclib_local']['W'].attrs.create('min_nz_sv', min_nz_sv)
                        fclib_file['fclib_local']['W'].attrs.create('norm_lsmr', norm_lsmr)
                        fclib_file['fclib_local']['W'].attrs.create('cond_lsmr', cond_lsmr)
                except Exception as e :
                    print "-->", e
                        
    if adhoc:
        print "script adhoc (convenient moulinette)"
        # for problem_filename in problem_filenames:
        #     print "treatment", problem_filename
        #     with h5py.File(problem_filename, 'r+') as fclib_file:
        #         try:
        #             import math
        #             rank =fclib_file['fclib_local']['W'].attrs.get('rank')
        #             if True :      #if rank == None:
        #                 rank_dense=fclib_file['fclib_local']['W'].attrs.get('rank_dense')
        #                 if rank_dense != None and not math.isnan(rank_dense):
        #                     print "rank := rank_dense"
        #                     fclib_file['fclib_local']['W'].attrs.create('rank', rank_dense)
        #                 rank_svd=fclib_file['fclib_local']['W'].attrs.get('rank_svd')
        #                 if rank_svd != None and not math.isnan(rank_svd):
        #                     print "rank := rank_svd"
        #                     fclib_file['fclib_local']['W'].attrs.create('rank', rank_svd)
                        
        #             else:
        #                 print "rank already present"

        #             r1 =fclib_file['fclib_local']['W'].attrs.get('r1')
        #             r2 =fclib_file['fclib_local']['W'].attrs.get('r2')
        #             if r1 != None:
        #                 print "r1 --> norm_lsmr"
        #                 fclib_file['fclib_local']['W'].attrs.create('norm_lsmr',r1)
        #                 fclib_file['fclib_local']['W'].attrs.__delitem__('r1')
        #             if r2 != None:
        #                 print "r2 --> cond_lsmr"
        #                 fclib_file['fclib_local']['W'].attrs.create('cond_lsmr',r2)
        #                 fclib_file['fclib_local']['W'].attrs.__delitem__('r2')
                    
        #         except Exception as e :
        #             print e
        #             pass
        with h5py.File('comp.hdf5', 'r+') as comp_file:
            data = comp_file['data']
            comp_data = data['comp']
            for solver in comp_data:
                if ( 'NSGS-AC-' in solver):
                    print "solver", solver
                    if ('NSGS-AC-GP' not in solver):
                        print "solver to be renamed", solver
                        new_solver= solver.replace("AC-","AC-GP-")
                        print "rename", solver, "in ",new_solver
                        data['comp'].move(solver,new_solver)

    if compute_hardness:
        nc = []
        nds = []
        cond_nc = []
        for problem_filename in problem_filenames:

            try:
                nc.append(numberOfDegreeofFreedomContacts(problem_filename)/3)
            except:
                pass
            try:
                nds.append(numberOfDegreeofFreedom(problem_filename))
            except:
                pass
            try:
                cond_nc.append(cond_problem(problem_filename))
            except:
                pass
        # compute other quantities
        print nc
        nc_avg = sum(nc)/float(len(nc))
        print "nc_avg", nc_avg
        with h5py.File('comp.hdf5', 'r') as comp_file:
            data = comp_file['data']
            comp_data = data['comp']
            for solver in solvers:
                solver_name=solver.name()


                if solver_name in comp_data :
                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, None, overwrite=False)
                    assert len(filenames) <= n_problems
                    measure[solver_name] = np.inf * np.ones(n_problems)
                    solver_r[solver_name] = np.inf * np.ones(n_problems)

                    ip = 0

                    for filename in filenames:
                        if filename not in min_measure:
                            min_measure[filename] = np.inf
                        try:
                            pfilename = os.path.splitext(filename)[0]
                            if comp_data[solver_name][pfilename].attrs['info'] == 0:
                                measure[solver_name][ip] =  comp_data[solver_name][pfilename].attrs[measure_name]
                                min_measure[filename] = min(min_measure[filename], measure[solver_name][ip])
                            else:
                                measure[solver_name][ip] = np.inf
                        except:
                            measure[solver_name][ip] = np.nan
                        ip += 1

            print "min_measure", min_measure
            avg_min_measure=0.0
            for k,v in min_measure.items():
                avg_min_measure +=v

            avg_min_measure = avg_min_measure/float(len(min_measure))
            #print         "avg_min_measure",avg_min_measure
            print         "Average min resolution measure by contact = {0:12.8e}".format(avg_min_measure/nc_avg)

               
                    
    if display_distrib:
        from matplotlib.pyplot import title, subplot, grid, show, legend, figure, hist, xlim, ylim, xscale
        if display_distrib_var == 'from-files':

            nc = []
            nds = []
            cond_nc = []
            cond_W = []
            rank_dense_W=[]
            rank_ratio =[]
            for problem_filename in problem_filenames:

                try:
                    nc.append(numberOfDegreeofFreedomContacts(problem_filename)/3)
                except:
                    pass
                try:
                    nds.append(numberOfDegreeofFreedom(problem_filename))
                except:
                    pass
                try:
                    cond_nc.append(cond_problem(problem_filename))
                except:
                    pass
                try:
                    cond_W.append(cond(problem_filename))
                except:
                    pass
                try:
                    rank_dense_W.append(rank_dense(problem_filename))
                except:
                    pass
                try:
                    rank_ratio.append(numberOfDegreeofFreedomContacts(problem_filename)/float(rank_dense(problem_filename)))
                except:
                    pass
            print "nds", nds
            print "rank_dense_W",  rank_dense_W
            print "cond_nc", cond_nc
            print "cond_W", cond_W

                
            figure()
            subplot(311)
            hist(nc, 100, label='nc', histtype='stepfilled')
            grid()
            legend()
            subplot(312)
            import math
            if not math.isnan(min(nds)):
                hist(nds, 100, label='nds', histtype='stepfilled')
            grid()
            legend()
            subplot(313)
            
            if not math.isnan(min(cond_nc)):
                hist(cond_nc, 100, label='cond_nc', histtype='stepfilled')
            grid()
            legend()

            figure()
            subplot(311)
            if not math.isnan(min(cond_W)):
                hist(cond_W, 100, label='cond(W)', histtype='stepfilled')
            grid()
            legend()
            subplot(312)
            if not math.isnan(min(rank_dense_W)):
                hist(rank_dense_W, 100, label='rank(W)', histtype='stepfilled')
            grid()
            legend()
            subplot(313)
            
            if not math.isnan(min(rank_ratio)):
                hist(rank_ratio, 100, label='rank_ratio(W)', histtype='stepfilled')
            grid()
            legend()

            



            
            if gnuplot_distrib :

                with open('distrib.gp','w') as gp:
                    #  all_rhos = [ domain ] + [ rhos[solver_name] for solver_name in comp_data ]
                    all_distrib = [ nc] + [nds] + [cond_nc]
                    np.savetxt('distrib.dat', np.matrix(all_distrib).transpose())
                    gp.write('resultfile = "distrib.dat"\n')
                    gp.write('basename="distrib-{0}"\n'.format(problem_filename.partition('-')[0]))
                    gp.write('\n')
                    gp.write('ismin(x) = (x<min)?min=x:0\n')
                    gp.write('ismax(x) = (x>max)?max=x:0\n')
                    gp.write('\n')
                    gp.write('max=-1e38;min=1e38;\n')
                    gp.write('plot resultfile u 1:(ismin($1)*ismax($1))\n')
                    gp.write('min_nc = min; max_nc = max\n')
                    gp.write('max=-1e38;min=1e38;\n')
                    gp.write('plot resultfile u 1:(ismin($2)*ismax($2))\n')
                    gp.write('min_ndof = min; max_ndof = max\n')
                    gp.write('max=-1e38;min=1e38;\n')
                    gp.write('plot resultfile u 1:(ismin($3)*ismax($3))\n')
                    gp.write('min_ncond = min; max_ncond = max\n')
                    gp.write('\n')
                    gp.write('\n')

                    gp.write('term_choice_tikz=1\n')
                    gp.write('if (term_choice_tikz == 1) \\\n')
                    gp.write('set term tikz standalone monochrome  size 5in,3in font \'\\small\\sf\';  \\\n')
                    gp.write('extension = \'.tex\'; \\\n')
                    gp.write('set output basename.extension; \\\n')
                    gp.write('print "output = ", basename.extension; \\\n')
                    gp.write('else \\\n')
                    gp.write('set term aqua;\\\n')
                    gp.write(' \n')
                    gp.write('set xtics offset 0,0.5 \n')
                    gp.write('set key left top\n')

                    gp.write('basheight = 0.36; heightoff = 0.0; winratio = 1.0; winheight = basheight*winratio ; trans = 0.9\n')
                    gp.write('set multiplot \n')
                    gp.write('set size winratio,winheight \n')
                    gp.write('\n')
                    # gp.write('set xrange [{0}:{1}]\n'.format(domain[0]-0.01, domain[len(domain)-1]))
                    # gp.write('set yrange [-0.01:1.01]\n')
                    gp.write('set ylabel \'\\shortstack{Number of \\\ problems} \' \n')
                    gp.write('bin(x, width) = width*floor(x/width) + binwidth/2.0\n')
                    #gp.write('set title \'{0}\'\n'.format(problem_filename.partition('-')[0]));
                    gp.write('\n')
                    gp.write('set origin 0.0,winheight*2.0*trans+heightoff\n')
                    gp.write('numberofbox=50\n')
                    gp.write('print \'max_nc =\', max_nc,\'min_nc =\', min_nc \n')
                    gp.write('binwidth = (max_nc-min_nc)/numberofbox\n')
                    gp.write('set boxwidth binwidth\n')

                    gp.write('print \'binwidth =\', binwidth \n')

                    gp.write('set xlabel \'number of contacts\' offset 0,1.2 \n')
                    #gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes title \'number of contacts\'  \n')
                    gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes notitle  \n')
                    gp.write('\n')

                    gp.write('binwidth = (max_ndof-min_ndof)/numberofbox\n')
                    gp.write('print \'binwidth =\', binwidth \n')

                    gp.write('set boxwidth binwidth\n')
                    gp.write('set origin 0.0,winheight*1.0*trans+heightoff\n')

                    gp.write('set xlabel \'number of degrees of freedom \' offset 0,1.2 \n')

                    #gp.write('plot resultfile u (bin($2, binwidth)):(1.0) smooth freq w boxes title  \'number of degrees of freedom \' \n')
                    gp.write('plot resultfile u (bin($2, binwidth)):(1.0) smooth freq w boxes notitle \n')
                    gp.write('\n')

                    gp.write('set origin 0.0,winheight*0.0+heightoff\n')
                    gp.write('binwidth = (max_ncond-min_ncond)/numberofbox\n')
                    gp.write('print \'binwidth =\', binwidth \n')

                    gp.write('set boxwidth binwidth\n')
                    gp.write('set xlabel \'ratio number of contacts unknowns/number of degrees of freedom\' offset 0,1.2 \n')
                    #gp.write('plot resultfile u (bin($3, binwidth)):(1.0) smooth freq w boxes title \'ratio number of contacts unknowns/number of degrees of freedom\' \n')
                    gp.write('plot resultfile u (bin($3, binwidth)):(1.0) smooth freq w boxes notitle \n')



                    #gp.write(','.join(['resultfile using 1:{0} t "{1}" w l'.format(index + 2, solver.name())
                    #                    for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
                # all_rhos = [ rhos[solver_name] for solver_name in comp_data ]
                # g.plot(*all_rhos)

        else:
            with h5py.File('comp.hdf5', 'r') as comp_file:

                data = comp_file['data']
                comp_data = data['comp']
                for solver_name in comp_data:

                    if user_filenames == []:
                        filenames = subsample_problems(comp_data[solver_name],
                                                       random_sample_proba,
                                                       max_problems, cond_nc)
                    else:
                        filenames = user_filenames

                    x = dict()
                    for filename in filenames:
                        pfilename = os.path.splitext(filename)[0]
                        if filename not in x:
                            try:
                                x[filename + solver_name] = comp_data[solver_name][pfilename].attrs[display_distrib_var]
                            except:
                                if display_distrib_var == 'cond-nc':
                                    print filename
                                    x[filename + solver_name] = cond_problem(filename)
                                else:
                                    x[filename + solver_name] = np.nan
                figure()
                l = [x[k] for k in x]
                l.sort()
                values = array(l)
                hist(values, 100, range=(min(values), max(values)), histtype='stepfilled')
                grid()

    if display or display_convergence or display_distrib:
        if not no_matplot:
            show()
