#!/usr/bin/env python

# parallel usage :
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--file={}'

#
# comp.py --max-problems=10 --no-compute --no-collect # output problems.txt
# cat problems.txt | parallel comp.py --timeout=100 --no-collect '--file={}'
#
#



from glob import glob
from itertools import product
import numpy as np
import random
import Siconos.Numerics as N
N.setNumericsVerbose(0)
import Siconos.FCLib as FCL

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

def subsample_problems(filenames, proba, maxp, cond):

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

    # overwrite!
    with open('problems.txt','w') as problems_txt:
        problems_txt.write('{0}\n'.format('\n'.join(r)))

    return r

def extern_guess(problem_filename, solver_name, iteration, h5file):
    data = h5file['data']
    comp_data = data['comp']

    reaction = comp_data[solver_name][problem_filename]['reactions'][iteration]
    velocity = comp_data[solver_name][problem_filename]['velocities'][iteration]
    return reaction, velocity

# estimate of condition number and norm from lsmr
# http://www.stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf
def _norm_cond(problem_filename):
    problem = read_fclib_format(problem_filename)[1]
    A = csr_matrix(N.SBMtoSparse(problem.M)[1])
    r = lsmr(A, np.ones([A.shape[0], 1]))  # solve Ax = 1
    svd = svds(A)[1]
    return r[5], r[6], max(svd), min(svd), max(svd)/min(svd)

norm_cond = Memoize(_norm_cond)

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
  print "   set the maximum time of iteration for each problem to n iterations (default",maxiter,")"
  print " --domain='a:d:b'"
  print "   restrict the domain of the performance profile to the interval [a,b] with a step of d (default",domain[0],":",domain[1]-domain[0],":",domain[-1]+domain[1]-domain[0],")"
  print "   or a perfomance profile a should be greater or equal 1"
  print " --measure=value"
  print "   select the value  as the measure for the perfomance profile. Possible values are time, iter, flpops"
  print " --display"
  print "   perform the computation of performance profile and display it in matplotlib"
  print " --new"
  print "   remove comp.hdf5 file"
  print " --solvers=string"
  print "   use keyworks in s separated by comma for filtering solvers"
  print " --solvers-exact=string"
  print "   use exact names of solvers in s separated by comma for filtering solvers"
  print " --with_mumps= 0 or 1"
  print "   use mumps as linear system solver"
  print " --max-problems=<max>"
  print "   Randomly select <max> problems in current directory." 
  print "   The problems list is written in problems.txt file"

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
                                    'clean', 'display', 'display-convergence',
                                    'files=', 'solvers-exact=', 'solvers=',
                                    'random-sample=', 'max-problems=',
                                    'timeout=', 'maxiter=', 'precision=',
                                    'keep-files', 'new', 'errors',
                                    'velocities', 'reactions', 'measure=',
                                    'just-collect', 'cond-nc=', 'display-distrib=',
                                    'no-collect', 'no-compute', 'domain=', 'replace-solver=',
                                    'gnuplot-profile','gnuplot-distrib', 'logscale', 'gnuplot-separate-keys',
                                    'output-dat', 'with_mumps=', 'file-filter='])


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
    elif o == '--precision':
        precision = float(a)
    elif o == '--clean':
        clean = True
    elif o == '--display':
        display = True
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
    elif o == '--solvers':
        user_solvers = split(a, ',')
    elif o == '--solvers-exact':
        user_solvers_exact = split(a, ',')
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
    elif o == '--replace-solver':
        try:
            with h5py.File('comp.hdf5','r+') as comp_file:
                print list(comp_file['data'])
                del comp_file['data']['comp'][a]
        except Exception as e:
            print e
    elif o == '--gnuplot-profile':
        gnuplot_profile=True
    elif o == '--logscale':
        logscale=True
    elif o == '--gnuplot-distrib':
        gnuplot_distrib=True
    elif o == 'gnuplot-separate-keys':
        gnuplot_separate_keys = True
    elif o == '--output-dat':
        output_dat=True
    elif o == '--with_mumps':
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

from ctypes import cdll, c_float, c_longlong, byref
try:
    with_papi = True
    papi=cdll.LoadLibrary('/usr/local/lib/libpapi.so')
except:
    try:
        with_papi = True
        papi=cdll.LoadLibrary('/usr/lib/x86_64-linux-gnu/libpapi.so.5.3.0.0')
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


def _numberOfInvolvedDS(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['info'].attrs['numberOfInvolvedDS']
        except:
            r = np.nan
    return r

def _dimension(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['spacedim'][0]
        except:
            r = np.nan
    return r

def _numberOfContacts(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W']['m'][0]
        except:
            r = np.nan
    return r

def _cond_problem(filename):
    problem = read_fclib_format(filename)[1]
    return float(problem.numberOfContacts * 3) / float(numberOfInvolvedDS(filename) * 6)


read_fclib_format = Memoize(_read_fclib_format)

numberOfInvolvedDS = Memoize(_numberOfInvolvedDS)

numberOfContacts = Memoize(_numberOfContacts)

dimension = Memoize(_dimension)

cond_problem = Memoize(_cond_problem)

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

            print '--->',e

            try:
               os.remove(output_filename)
            except:
                pass

            with h5py.File(output_filename, 'w') as output:

                digest = hashlib.sha256(open(filename, 'rb').read()).digest()
                data = output.create_group('data')
                comp_data = data.create_group('comp')
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
                attrs.create('nc', numberOfContacts(filename))
                attrs.create('nds', numberOfInvolvedDS(filename))
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

                print(filename, numberOfContacts(filename), numberOfInvolvedDS(filename), cond_problem(filename), solver.name(), info, iter, err,
                      time_s, real_time, proc_time,
                      flpops, mflops)

                with open('report.txt', "a") as report_file:
                    print   >> report_file , (filename, solver.name(), info, iter, err,
                                              time_s, real_time, proc_time,
                                              flpops, mflops)



    @timeout(utimeout)
    def _internal_call(self, solver, problem, filename, pfilename, output_filename):


        with h5py.File(output_filename, 'w') as output:


            data = output.create_group('data')
            comp_data = data.create_group('comp')
            solver_data = comp_data.create_group(solver.name())

            solver_problem_data = solver_data.create_group(pfilename)
            attrs = solver_problem_data.attrs

            psize = dimension(filename) * numberOfContacts(filename)

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

            _, guess_err = N.FrictionContact3D_compute_error(read_fclib_format(filename)[1],
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
            attrs.create('nc', numberOfContacts(filename))
            attrs.create('nds', numberOfInvolvedDS(filename))
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

            # filename, solver name, revision svn, parameters, nb iter, err
            print(filename, numberOfContacts(filename), numberOfInvolvedDS(filename), cond_problem(filename), solver.name(), info, iter, err,
                  time_s, real_time, proc_time,
                  flpops, mflops)

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
            psize = problem.dimension * problem.numberOfContacts
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

bogusPureNewton = BogusSolver(name="BogusPureNewton", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_LOCALFB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusPureNewton.SolverOptions().iparam[4]=0


bogusPureEnumerative = BogusSolver(name="BogusPureEnumerative", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_LOCALFB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusPureEnumerative.SolverOptions().iparam[4]=1


bogusHybrid = BogusSolver(name="BogusHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_LOCALFB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
bogusHybrid.SolverOptions().iparam[4]=2


bogusRevHybrid = BogusSolver(name="BogusRevHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_LOCALFB, iparam_iter=1, dparam_err=1, maxiter=maxiter, precision=precision)
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
localACSTD = SiconosSolver(name="NSN-AlartCurnier",
                           gnuplot_name="NSN-AC",
                           API=N.frictionContact3D_localAlartCurnier,
                           TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                           iparam_iter=1,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)

localACSTD.SolverOptions().iparam[10] = 0;
localACSTD.SolverOptions().iparam[11] = 0;
localACSTD.SolverOptions().iparam[12] = 10;
localACSTD.SolverOptions().iparam[13] = with_mumps;
localACSTD.SolverOptions().iparam[3] = 10000000


localACJeanMoreau = SiconosSolver(name="NSN-JeanMoreau",
                                  gnuplot_name="NSN-JM",
                                  API=N.frictionContact3D_localAlartCurnier,
                                  TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                  iparam_iter=1,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

localACJeanMoreau.SolverOptions().iparam[10] = 1;
localACJeanMoreau.SolverOptions().iparam[11] = 0;
localACJeanMoreau.SolverOptions().iparam[12] = 10;
localACJeanMoreau.SolverOptions().iparam[13] = with_mumps;
localACJeanMoreau.SolverOptions().iparam[3] = 10000000

localACSTDGenerated = SiconosSolver(name="NSN-AlartCurnier-Generated",
                                    gnuplot_name="NSN-AC-Generated",
                                    API=N.frictionContact3D_localAlartCurnier,
                                    TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                    iparam_iter=1,
                                    dparam_err=1,
                                    maxiter=maxiter, precision=precision)

localACSTDGenerated.SolverOptions().iparam[10] = 2;
localACSTDGenerated.SolverOptions().iparam[11] = 0;
localACSTDGenerated.SolverOptions().iparam[12] = 10;
localACSTDGenerated.SolverOptions().iparam[13] = with_mumps;
localACSTDGenerated.SolverOptions().iparam[3] = 10000000

localACJeanMoreauGenerated = SiconosSolver(name="NSN-JeanMoreau-Generated",
                                           gnuplot_name="NSN-JM-Generated",
                                           API=N.frictionContact3D_localAlartCurnier,
                                           TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                           iparam_iter=1,
                                           dparam_err=1,
                                           maxiter=maxiter, precision=precision)

localACJeanMoreauGenerated.SolverOptions().iparam[10] = 3;
localACJeanMoreauGenerated.SolverOptions().iparam[11] = 0;
localACJeanMoreauGenerated.SolverOptions().iparam[12] = 10;
localACJeanMoreauGenerated.SolverOptions().iparam[13] = with_mumps;
localACJeanMoreauGenerated.SolverOptions().iparam[3] = 10000000



localfb_gp = SiconosSolver(name="NSN-FischerBurmeister-GP",
                           gnuplot_name="NSN-FB-GP",
                           API=N.frictionContact3D_localFischerBurmeister,
                           TAG=N.SICONOS_FRICTION_3D_LOCALFB,
                           iparam_iter=1,
                           dparam_err=1,
                           maxiter=maxiter, precision=precision)

localfb_gp.SolverOptions().iparam[3] = 1000000
localfb_gp.SolverOptions().iparam[11] = 0
localfb_gp.SolverOptions().iparam[12] = 500
localfb_gp.SolverOptions().iparam[13] = with_mumps

localfb_fblsa = SiconosSolver(name="NSN-FischerBurmeister-FBLSA",
                              gnuplot_name="NSN-FB-FBLSA",
                              API=N.frictionContact3D_localFischerBurmeister,
                              TAG=N.SICONOS_FRICTION_3D_LOCALFB,
                              iparam_iter=1,
                              dparam_err=1,
                              maxiter=maxiter, precision=precision)

localfb_fblsa.SolverOptions().iparam[3] = 1000000
localfb_fblsa.SolverOptions().iparam[11] = 1
localfb_fblsa.SolverOptions().iparam[12] = 10
localfb_fblsa.SolverOptions().iparam[13] = with_mumps

localfb_nls = SiconosSolver(name="LocalFischerBurmeisterNLS",
                              API=N.frictionContact3D_localFischerBurmeister,
                              TAG=N.SICONOS_FRICTION_3D_LOCALFB,
                              iparam_iter=1,
                              dparam_err=1,
                              maxiter=maxiter, precision=precision)

localfb_nls.SolverOptions().iparam[3] = 1000000
localfb_nls.SolverOptions().iparam[11] = -1
localfb_nls.SolverOptions().iparam[12] = 10
localfb_nls.SolverOptions().iparam[13] = with_mumps


hlocalac = SiconosHybridSolver(name = "HLocalAlartCurnier",
                               API=N.frictionContact3D_localAlartCurnier,
                               TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                               iparam_iter=1,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)
hlocalac.SolverOptions().iparam[3] = 10000000

nsgs = SiconosSolver(name="NSGS-AC",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

snsgs = SiconosSolver(name="NSGS-AC-Shuffled",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
snsgs.SolverOptions().iparam[9] = 1

nsgs_pli = SiconosSolver(name="NSGS-PLI",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
nsgs_pli.SolverOptions().solverId = N.SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration


 
# only dense
nsgsv = SiconosSolver(name="NSGS-Velocity",
                      API=N.frictionContact3D_nsgs_velocity,
                      TAG=N.SICONOS_FRICTION_3D_NSGSV,
                      iparam_iter=7,
                      dparam_err=1,
                      maxiter=maxiter, precision=precision)


omega=1.5
psor = SiconosSolver(name="PSOR-AC",
                     gnuplot_name="PSOR-AC",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
psor.SolverOptions().iparam[8] = 1
psor.SolverOptions().dparam[8] = omega

omega_values = [0.5, 0.8, 1.0, 1.1, 1.3, 1.5, 1.8]
psor_series=[]
for omega in omega_values:
    psor_solver = SiconosSolver(name="PSOR-AC-"+str(omega),
                                gnuplot_name="PSOR-AC \$\\\omega="+str(omega)+"\$",
                                API=N.frictionContact3D_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=7,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)
    psor_solver.SolverOptions().iparam[8] = 1
    psor_solver.SolverOptions().dparam[8] = omega
    psor_series.append(psor_solver)

TrescaFixedPoint = SiconosSolver(name="TrescaFixedPoint-NSGS-PLI",
                                 API=N.frictionContact3D_TrescaFixedPoint,
                                 TAG=N.SICONOS_FRICTION_3D_TFP,
                                 iparam_iter=7,
                                 dparam_err=1,
                                 maxiter=maxiter, precision=precision)

ACLMFixedPoint = SiconosSolver(name="ACLMFixedPoint-SOCLCP-NSGS-PLI",
                               API=N.frictionContact3D_ACLMFixedPoint,
                               TAG=N.SICONOS_FRICTION_3D_ACLMFP,
                               iparam_iter=7,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)

SOCLCP = SiconosSolver(name="SOCLCP-NSGS-PLI",
                       API=N.frictionContact3D_SOCLCP,
                       TAG=N.SICONOS_FRICTION_3D_SOCLCP,
                       iparam_iter=7,
                       dparam_err=1,
                       maxiter=maxiter, precision=precision)

DeSaxceFixedPoint = SiconosSolver(name="FixedPoint-DeSaxce",
                                  API=N.frictionContact3D_DeSaxceFixedPoint,
                                  TAG=N.SICONOS_FRICTION_3D_DSFP,
                                  iparam_iter=7,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

ExtraGrad = SiconosSolver(name="ExtraGradient",
                          API=N.frictionContact3D_ExtraGradient,
                          TAG=N.SICONOS_FRICTION_3D_EG,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

FixedPointProjection = SiconosSolver(name="FixedPoint-Projection",
                                     API=N.frictionContact3D_fixedPointProjection,
                                     TAG=N.SICONOS_FRICTION_3D_FPP,
                                     iparam_iter=7,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)

VIExtraGrad = SiconosSolver(name="ExtraGradient-VI",
                            API=N.frictionContact3D_VI_ExtraGradient,
                            TAG=N.SICONOS_FRICTION_3D_VI_EG,
                            iparam_iter=7,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)
VIExtraGrad1 = SiconosSolver(name="ExtraGradient-VI1s",
                            API=N.frictionContact3D_VI_ExtraGradient,
                            TAG=N.SICONOS_FRICTION_3D_VI_EG,
                            iparam_iter=7,
                            dparam_err=1,
                            maxiter=maxiter, precision=precision)
VIExtraGrad1.SolverOptions().dparam[4]=0.6
VIExtraGrad1.SolverOptions().dparam[5]=1/0.7
VIExtraGrad1.SolverOptions().dparam[6]=0.9
VIExtraGrad1.SolverOptions().dparam[7]=0.3



iparam1_values = [0,1,2]
iparam2_values = [0,1]
iparam3_values = [0,1]
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
                g_name = g_name + " True" 
            elif i2 == 1:
                g_name = g_name + " False"
                
            if i3 == 0:
                g_name = g_name 
            elif i3 == 1:
                g_name = g_name + " min" 

            VIExtraGrad_solver= SiconosSolver(name="ExtraGrad-VI-"+str(i1)+str(i2)+str(i3),
                                                         gnuplot_name=g_name,
                                                         API=N.frictionContact3D_VI_ExtraGradient,
                                                         TAG=N.SICONOS_FRICTION_3D_VI_EG,
                                                         iparam_iter=7,
                                                         dparam_err=1,
                                                         maxiter=maxiter, precision=precision)
            VIExtraGrad_solver.SolverOptions().iparam[1] = i1
            VIExtraGrad_solver.SolverOptions().iparam[2] = i2
            VIExtraGrad_solver.SolverOptions().iparam[3] = i3
            VIExtraGrad_series.append(VIExtraGrad_solver)

            
VIFixedPointProjection = SiconosSolver(name="FixedPoint-VI",
                                       API=N.frictionContact3D_VI_FixedPointProjection,
                                       TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                       iparam_iter=7,
                                       dparam_err=1,
                                       maxiter=maxiter, precision=precision)

iparam1_values = [0,1,2]
iparam2_values = [0,1]
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
                g_name = g_name + " True" 
            elif i2 == 1:
                g_name = g_name + " False"
                
            if i3 == 0:
                g_name = g_name 
            elif i3 == 1:
                g_name = g_name + " min" 

            VIFixedPointProjection_solver= SiconosSolver(name="FixedPoint-VI-"+str(i1)+str(i2)+str(i3),
                                                         gnuplot_name=g_name,
                                                         API=N.frictionContact3D_VI_FixedPointProjection,
                                                         TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                                         iparam_iter=7,
                                                         dparam_err=1,
                                                         maxiter=maxiter, precision=precision)
            VIFixedPointProjection_solver.SolverOptions().iparam[1] = i1
            VIFixedPointProjection_solver.SolverOptions().iparam[2] = i2
            VIFixedPointProjection_solver.SolverOptions().iparam[3] = i3
            VIFixedPointProjection_series.append(VIFixedPointProjection_solver)

Prox = SiconosSolver(name="PROX-NSN-AC",
                     gnuplot_name="PPA-NSN-AC  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Prox.SolverOptions().internalSolvers.iparam[3] = 1000000

ProxFB = SiconosSolver(name="PROX-NSN-FB-GP",
                     gnuplot_name="PPA-NSN-FB-GP  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
ProxFB.SolverOptions().internalSolvers.iparam[3] = 1000000
ProxFB.SolverOptions().dparam[4]=5.0 # sigma
ProxFB.SolverOptions().dparam[5]=1.0 # nu
localfb_gp_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_LOCALFB)
localfb_gp_inprox.iparam[3] = 1000000
localfb_gp_inprox.iparam[11] = 0

localfb_gp_inprox.iparam[12] = 6

ProxFB.SolverOptions().internalSolvers = localfb_gp_inprox
ProxFB.SolverOptions().internalSolvers.iparam[3] = 1000000

ProxFB_fblsa = SiconosSolver(name="PROX-NSN-FB-FBLSA",
                     gnuplot_name="PPA-NSN-FB-FBLSA  \$ \\\mu=1, \\\sigma=5.0\$",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
localfb_fblsa_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_LOCALFB)
localfb_fblsa_inprox.iparam[3] = 1000000
localfb_fblsa_inprox.iparam[11] = 1
localfb_fblsa_inprox.iparam[12] = 6
ProxFB_fblsa.SolverOptions().internalSolvers = localfb_fblsa_inprox
ProxFB_fblsa.SolverOptions().internalSolvers.iparam[3] = 1000000




sigmavalues= [0.5, 1.0, 4.0, 5.0, 50, 100.0, 1000.0 ]
muvalues= [0.5, 1.0, 2.0]
prox_series =[]
for mu in muvalues:
    for sigma in sigmavalues:
        prox_solver  = SiconosSolver(name="PROX-NSN-AC-nu"+str(mu)+"-sigma"+str(sigma),
                                     gnuplot_name="PPA-NSN-AC  \$ \\\mu="+str(mu)+", \\\sigma="+str(sigma)+"\$",
                                     API=N.frictionContact3D_proximal,
                                     TAG=N.SICONOS_FRICTION_3D_PROX,
                                     iparam_iter=7,
                                     dparam_err=1,
                                     maxiter=maxiter, precision=precision)
        prox_solver.SolverOptions().internalSolvers.iparam[3] = 1000000
        prox_solver.SolverOptions().dparam[4]=sigma # sigma
        prox_solver.SolverOptions().dparam[5]=mu # nu
        prox_series.append(prox_solver)

localac_wrapped = SiconosSolver(name="NSN-AlartCurnier-Wrapped",
                                API=N.frictionContact3D_localAlartCurnier,
                                TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                iparam_iter=1,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)

def fc3d_localac_r(problem, reactions, velocities, _SO):
    SO = N.SolverOptions(N.SICONOS_FRICTION_3D_VI_FPP)
    SO.iparam[3] = 1000
    N.frictionContact3D_VI_FixedPointProjection(problem, reactions, velocities, SO)
    #    print '->',SO.dparam[3]
    localac_wrapped.SolverOptions().dparam[3] = SO.dparam[3]
    return localac_wrapped(problem, reactions, velocities)

# flop measure only on localac
localacr = SiconosWrappedSolver(name="NSN-AlartCurnier-R",
                                API=fc3d_localac_r,
                                TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                iparam_iter=1,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)

#
quartic = SiconosSolver(name="NSGS-Quartic",
                        API=N.frictionContact3D_nsgs,
                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                        iparam_iter=7,
                        dparam_err=1, maxiter=maxiter, precision=precision)

quartic3x3 = N.SolverOptions(N.SICONOS_FRICTION_3D_QUARTIC_NU)

quartic.SolverOptions().internalSolvers = quartic3x3

# 1 contact
#AlartCurnierNewton = SiconosSolver(name="AlartCurnierNewton",
#                                   API=frictionContact3D_AlartCurnierNewton,
#                                   iparam_iter=1,
#                                   iparam_iter=1)

# rho estimation needed
Prox._SO.dparam[3] = 1000

HyperplaneProjection = SiconosSolver(name="HyperplaneProjection",
                                     API=N.frictionContact3D_HyperplaneProjection,
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
#frictionContact3D_sparseGlobalAlartCurnierInit(localac.SolverOptions())

#all_solvers = [nsgs, snsgs, TrescaFixedPoint, localac, Prox, DeSaxceFixedPoint,
#               FixedPointProjection, VIFixedPointProjection, ExtraGrad, VIExtraGrad]
#all_solvers = [nsgs, snsgs, quartic, TrescaFixedPoint, ACLMFixedPoint, DeSaxceFixedPoint, VIFixedPointProjection, VIFixedPointProjection1, VIFixedPointProjection2, VIFixedPointProjection3, VIExtraGrad, SOCLCP, Prox, Prox2, Prox3, Prox4, Prox5, localACSTD, localACSTDGenerated,  localacr, localACJeanMoreau, localACJeanMoreauGenerated, localfb_gp, localfb_fblsa]

all_solvers = [nsgs, nsgs_pli, snsgs, quartic, psor,
               TrescaFixedPoint, DeSaxceFixedPoint,
               VIFixedPointProjection, VIExtraGrad,VIExtraGrad1,
               SOCLCP,
               localACSTD,localACSTDGenerated,  localacr, localACJeanMoreau, localACJeanMoreauGenerated,
               Prox,  ProxFB,
               ACLMFixedPoint, localfb_gp]

all_solver_unstable = [localfb_fblsa, ProxFB_fblsa]

all_solvers.extend(all_solver_unstable)

# specific studies of solvers.
#all_solvers.extend(VIFixedPointProjection_series)
#all_solvers.extend(VIExtraGrad_series)
#all_solvers.extend(psor_series)
#all_solvers.extend(prox_series)
all_solvers.remove(quartic)



solvers=[]
if user_solvers != []:
    #print "user_solvers", user_solvers
    solvers.extend( filter(lambda s: any(us in s._name for us in user_solvers), all_solvers))

if user_solvers_exact != []:
    #print "user_solvers_exact", user_solvers_exact
    solvers.extend(filter(lambda s: any(us ==  s._name  for us in user_solvers_exact), all_solvers))

if solvers == []:
    solvers= all_solvers

#solvers = [ProxFB]
print "Operations will be run for solvers :", [ s._name for s in solvers]
    
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
    
__problem_filenames = subsample_problems(all_filenames,
                                         random_sample_proba,
                                         max_problems, None)

_problem_filenames = filter(is_fclib_file,
                           __problem_filenames)


problem_filenames = subsample_problems(_problem_filenames,
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

def collect(tpl):

    solver, filename = tpl

    pfilename = os.path.basename(os.path.splitext(filename)[0])
    results_filename = '{0}-{1}.hdf5'.format(solver.name(),pfilename)
    if os.path.exists(results_filename) and not os.stat(results_filename).st_size == 0:
        try:
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
            print (r.attrs['filename'], cond_problem(r.attrs['filename']), solver.name(), r.attrs['info'], r.attrs['iter'], r.attrs['err'], r.attrs['time'], r.attrs['real_time'], r.attrs['proc_time'], r.attrs['flpops'], r.attrs['mflops'])
            return False
        except:
            return True

if __name__ == '__main__':

    if compute:
        all_tasks = [t for t in product(solvers, problem_filenames)]

        if os.path.exists('comp.hdf5'):
            with h5py.File('comp.hdf5', 'r') as result_file:
                tasks = filter(Results(result_file), all_tasks)

        else:
            tasks = all_tasks

        if ask_compute:
            r = map(caller, tasks)

        if ask_collect:
            map(collect, tasks)

    if display:
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']

            # 1 n_problems
            n_problems = 0

            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
                    n_problems = max(n_problems, len(filenames))

            # 2 measures & min_measure

            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :

                    filenames = subsample_problems(comp_data[solver_name],
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
                    filenames = subsample_problems(comp_data[solver_name],
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


            if gnuplot_profile :
                def write_report(r, filename):
                    with open(filename, "w") as input_file:
                        for k, v in r.items():
                            line = '{}, {}'.format(k, v)
                            print >> input_file, line

                out_data=np.empty([len(domain),len(comp_data)+1])
                write_report(rhos,'rhos.txt')
                write_report(solver_r,'solver_r.txt')


                with open('profile.gp','w') as gp:
                    #  all_rhos = [ domain ] + [ rhos[solver_name] for solver_name in comp_data ]
                    all_rhos = [ domain ] + [ rhos[solver.name()] for solver in filter(lambda s: s._name in comp_data, solvers) ]
                    np.savetxt('profile.dat', np.matrix(all_rhos).transpose())
                    gp.write('resultfile = "profile.dat"\n')
                    gp.write('basename="profile-{0}"\n'.format(filename.partition('-')[0]))
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
                    print filename.partition('-')[0]
                    if logscale:
                        gp.write('set logscale x\n')
                        gp.write('set xlabel \'$\\tau$ ({0}) (logscale)\' \n'.format(measure_name))
                    else:
                        gp.write('set xlabel \'$\\tau$ ({0})\' \n'.format(measure_name))

                    #gp.write('set title \'{0}\'\n'.format(filename.partition('-')[0]));
                    gp.write('plot ')
                    gnuplot_separate_keys =True
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

    if display_distrib:
        from matplotlib.pyplot import title, subplot, grid, show, legend, figure, hist
        if display_distrib_var == 'from-files':

            nc = []
            nds = []
            cond_nc = []

            for problem_filename in problem_filenames:

                try:
                    nc.append(read_fclib_format(problem_filename)[1].numberOfContacts)
                except:
                    pass
                try:
                    nds.append(6*numberOfInvolvedDS(problem_filename))
                except:
                    pass
                try:
                    cond_nc.append(cond_problem(problem_filename))
                except:
                    pass

            figure()
            subplot(311)
            hist(nc, 100, label='nc', histtype='stepfilled')
            grid()
            legend()
            subplot(312)
            hist(nds, 100, label='nds', histtype='stepfilled')
            grid()
            legend()
            subplot(313)
            hist(cond_nc, 100, label='cond_nc', histtype='stepfilled')
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
        show()
