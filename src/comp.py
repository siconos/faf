#!/usr/bin/env python

# parallel usage :
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--file={}'

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
        _r = random.sample(__r, maxp)
    else:
        _r = __r

    if cond is not None:
        r = filter(WithCriterium(cond[0], cond[1]), _r)
    else:
        r = _r
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
output_dat=False
user_filenames = []
user_solvers = []
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
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'flop', 'iter', 'time',
                                    'clean', 'display', 'display-convergence',
                                    'files=', 'solvers=',
                                    'random-sample=', 'max-problems=',
                                    'timeout=', 'maxiter=', 'precision=',
                                    'keep-files', 'new', 'errors',
                                    'velocities', 'reactions', 'measure=',
                                    'just-collect', 'cond-nc=', 'display-distrib=',
                                    'no-collect', 'domain=', 'replace-solver=',
                                    'gnuplot-profile','gnuplot-distrib', 'logscale',
                                    'output-dat' ])


except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)
for o, a in opts:
    if o == '--flop':
        measure = 'flop'
    elif o == '--iter':
        measure = 'iter'
    elif o == '--time':
        measure = 'time'
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
        random_sample_proba = float(a)
    elif o == '--max-problems':
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
    elif o == '--just-collect':
        ask_compute = False
    elif o == '--no-collect':
        ask_collect = False
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
    elif o == '--output-dat':
        output_dat=True
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
                    result = solver(problem, reactions, velocities)

                    time_s = time.clock() - t0 # on unix, t is CPU seconds elapsed (floating point)

                    fclib_sol = FCL.fclib_solution()

                    fclib_sol.v = None
                    fclib_sol.r = reactions
                    fclib_sol.u = velocities
                    fclib_sol.l = None

                    nerr = FCL.fclib_merit_local(read_fclib_format(filename)[0],
                                                 FCL.MERIT_1, fclib_sol)

                    _, xerr = N.FrictionContact3D_compute_error(read_fclib_format(filename)[1],
                                                                reactions, velocities, precision, solver.SolverOptions())

                    print nerr, xerr

                    i_info, i_iter, i_err, i_real_time, i_proc_time, i_flpops, i_mflops = result

                    info = i_info
                    iter += i_iter
                    err = i_err
                    real_time += i_real_time
                    proc_time += i_proc_time
                    flpops += i_flpops
                    mflops = (mflops + i_mflops)/2.

                    if info == 0 and xerr >= precision:
#                        solver.SolverOptions().iparam[0]=1
                        solver.SolverOptions().dparam[0]=solver.SolverOptions().dparam[1]/10
                        again = False
                    else:
                        again = False


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

    def __init__(self, name=None, API=None, TAG=None, iparam_iter=None,
                 dparam_err=None, maxiter=maxiter, precision=precision):
        self._name = name
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
            if 'guesses' in f:
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
            return extern_guess(pfilename, 'NonsmoothGaussSeidel', 1, comp_file)


class SiconosWrappedSolver(SiconosSolver):
    def __call__(self, problem, reactions, velocities):
        return self._API(problem, reactions, velocities, self._SO)


#
# Some solvers
#
localACSTD = SiconosSolver(name="LocalAlartCurnierSTD",
                                  API=N.frictionContact3D_localAlartCurnier,
                                  TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                  iparam_iter=1,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

localACSTD.SolverOptions().iparam[10] = 0;


localACJeanMoreau = SiconosSolver(name="LocalAlartCurnierJeanMoreau",
                                  API=N.frictionContact3D_localAlartCurnier,
                                  TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                  iparam_iter=1,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

localACJeanMoreau.SolverOptions().iparam[10] = 1;

localACSTDGenerated = SiconosSolver(name="LocalAlartCurnierSTDGenerated",
                                    API=N.frictionContact3D_localAlartCurnier,
                                    TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                    iparam_iter=1,
                                    dparam_err=1,
                                    maxiter=maxiter, precision=precision)

localACSTDGenerated.SolverOptions().iparam[10] = 2;

localACJeanMoreauGenerated = SiconosSolver(name="LocalAlartCurnierJeanMoreauGenerated",
                                           API=N.frictionContact3D_localAlartCurnier,
                                           TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                           iparam_iter=1,
                                           dparam_err=1,
                                           maxiter=maxiter, precision=precision)

localACJeanMoreauGenerated.SolverOptions().iparam[10] = 3;




localfb = SiconosSolver(name="LocalFischerBurmeister",
                        API=N.frictionContact3D_localFischerBurmeister,
                        TAG=N.SICONOS_FRICTION_3D_LOCALFB,
                        iparam_iter=1,
                        dparam_err=1,
                        maxiter=maxiter, precision=precision)

#localac.SolverOptions().iparam[3] = 10000000


hlocalac = SiconosHybridSolver(name = "HLocalAlartCurnier",
                               API=N.frictionContact3D_localAlartCurnier,
                               TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                               iparam_iter=1,
                               dparam_err=1,
                               maxiter=maxiter, precision=precision)

hlocalac.SolverOptions().iparam[3] = 10000000

nsgs = SiconosSolver(name="NonsmoothGaussSeidel",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

snsgs = SiconosSolver(name="ShuffledNonsmoothGaussSeidel",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

snsgs.SolverOptions().iparam[9] = 1

# only dense
nsgsv = SiconosSolver(name="NonsmoothGaussSeidelVelocity",
                      API=N.frictionContact3D_nsgs_velocity,
                      TAG=N.SICONOS_FRICTION_3D_NSGSV,
                      iparam_iter=7,
                      dparam_err=1,
                      maxiter=maxiter, precision=precision)

TrescaFixedPoint = SiconosSolver(name="TrescaFixedPoint",
                                 API=N.frictionContact3D_TrescaFixedPoint,
                                 TAG=N.SICONOS_FRICTION_3D_TFP,
                                 iparam_iter=7,
                                 dparam_err=1,
                                 maxiter=maxiter, precision=precision)

DeSaxceFixedPoint = SiconosSolver(name="DeSaxceFixedPoint",
                                  API=N.frictionContact3D_DeSaxceFixedPoint,
                                  TAG=N.SICONOS_FRICTION_3D_DSFP,
                                  iparam_iter=7,
                                  dparam_err=1,
                                  maxiter=maxiter, precision=precision)

Prox = SiconosSolver(name="ProximalFixedPoint",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Prox.SolverOptions().internalSolvers.iparam[3] = 1000000

Prox2 = SiconosSolver(name="ProximalFixedPoint2",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

Prox2.SolverOptions().dparam[5]=2.0 # nu



Prox2.SolverOptions().internalSolvers.iparam[3] = 1000000

Prox3 = SiconosSolver(name="ProximalFixedPoint3",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

Prox3.SolverOptions().dparam[4]=1.1 # sigma
Prox3.SolverOptions().dparam[5]=2.0 # nu

Prox3.SolverOptions().internalSolvers.iparam[3] = 1000000


Prox4 = SiconosSolver(name="ProximalFixedPoint4",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)
Prox4.SolverOptions().dparam[4]=100 # sigma
Prox4.SolverOptions().dparam[5]=1.0 # nu

Prox4.SolverOptions().internalSolvers.iparam[3] = 1000000


Prox5 = SiconosSolver(name="ProximalFixedPoint5",
                     API=N.frictionContact3D_proximal,
                     TAG=N.SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)

Prox5.SolverOptions().dparam[4]=1000.0 # sigma
Prox5.SolverOptions().dparam[5]=1.0 # nu

Prox5.SolverOptions().internalSolvers.iparam[3] = 1000000


ExtraGrad = SiconosSolver(name="ExtraGradient",
                          API=N.frictionContact3D_ExtraGradient,
                          TAG=N.SICONOS_FRICTION_3D_EG,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

FixedPointProjection = SiconosSolver(name="FixedPointProjection",
                          API=N.frictionContact3D_fixedPointProjection,
                          TAG=N.SICONOS_FRICTION_3D_FPP,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

VIExtraGrad = SiconosSolver(name="VIExtraGradient",
                          API=N.frictionContact3D_VI_ExtraGradient,
                          TAG=N.SICONOS_FRICTION_3D_VI_EG,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

VIFixedPointProjection = SiconosSolver(name="VIFixedPointProjection",
                          API=N.frictionContact3D_VI_FixedPointProjection,
                          TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)



localac_wrapped = SiconosSolver(name="LocalAlartCurnierWrapped",
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
localacr = SiconosWrappedSolver(name="LocalacR",
                                API=fc3d_localac_r,
                                TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                                iparam_iter=1,
                                dparam_err=1,
                                maxiter=maxiter, precision=precision)

#
quartic = SiconosSolver(name="NonsmoothGaussSeidelQuartic",
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


all_solvers = [nsgs, snsgs, TrescaFixedPoint, Prox, Prox2, Prox3, Prox4, Prox5, localACSTD, localACSTDGenerated, localACJeanMoreau, localACJeanMoreauGenerated, localfb, localacr, DeSaxceFixedPoint, VIFixedPointProjection, VIExtraGrad, bogusPureEnumerative, bogusPureNewton, bogusHybrid, bogusRevHybrid, quartic]


if user_solvers == []:
    solvers = all_solvers
else:
    solvers = filter(lambda s: any(us in s._name for us in user_solvers), all_solvers)


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


if user_filenames == []:
    all_filenames = glob('*.hdf5')
else:
    all_filenames = user_filenames

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
                    gp.write('set output basename.extension; \\\n')
                    gp.write('print "output = ", basename.extension; \\\n')

                    gp.write('else \\\n')
                    gp.write('set term aqua;\\\n')
                    gp.write('\n')
                    gp.write('set xrange [{0}:{1}]\n'.format(domain[0]-0.01, domain[len(domain)-1]))
                    gp.write('set yrange [-0.01:1.01]\n')
                    gp.write('set ylabel \'$\\rho(\\tau)$ \' \n')
                    gp.write('set key below right vertical maxrows 4\n')
                    print filename.partition('-')[0]
                    if logscale:
                        gp.write('set logscale x\n')
                        gp.write('set xlabel \'$\\tau$ ({0}) (logscale)\' \n'.format(measure_name))
                    else:
                        gp.write('set xlabel \'$\\tau$ ({0})\' \n'.format(measure_name))

                    #gp.write('set title \'{0}\'\n'.format(filename.partition('-')[0]));
                    gp.write('plot ')
                    gp.write(','.join(['resultfile using 1:{0} t "{1}" w l'.format(index + 2, solver.name())
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
