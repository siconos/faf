#!/usr/bin/env python

# parallel usage : 
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--file={}'

from glob import glob
from itertools import product
import numpy as np
import Siconos.Numerics as N
#import Siconos.FCLib as FCL

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
display = False
display_convergence = False
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
maxiter = 100000
precision = 1e-8
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'flop', 'iter', 'time', 'clean','display','display-convergence','files=','solvers=',
                                'timeout=', 'maxiter=', 'precision=', 'keep-files', 'new', 'errors', 'velocities', 'reactions', 'measure=', 'just-collect', 'no-collect'])
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
    elif o == '--display-convergence':
        display_convergence = True
    elif o == '--measure':
        measure_name = a
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
    elif o == '--new':
        try:
            os.remove('comp.hdf5')
        except:
            pass

    elif o == '--files':

        files = split(a,',')

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


def read_fclib_format(f):
    #    fc_problem = FCL.fclib_read_local(f)
    #solution = FCL.fclib_read_solution(f)
    #print FCL.fclib_merit_local(fc_problem, FCL.MERIT_1, solution)
    #print fc_problem.W.m
    #print solution.u
    #solution.u = np.zeros(fc_problem.W.m * fc_problem.spacedim)
    #solution.r = np.zeros(fc_problem.W.m * fc_problem.spacedim)
    #print FCL.fclib_merit_local(fc_problem, FCL.MERIT_1, solution)

    return N.frictionContact_fclib_read(f)


pread_fclib_format = Memoize(read_fclib_format)

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
        problem = pread_fclib_format(filename)

        pfilename = os.path.splitext(filename)[0]

        output_filename = '{0}-{1}.hdf5'.format(solver.name(),
                                                pfilename)
        # considered as tmp file
        try:
                os.remove(output_filename)
        except:
                pass

        try:
            self._internal_call(solver, problem, filename, pfilename, output_filename)

        except Exception as e:

            print e

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
                attrs.create('digest', digest)
                attrs.create('info', info)
                attrs.create('iter', iter)
                attrs.create('err', err)
                attrs.create('time', time_s)
                attrs.create('real_time', real_time)
                attrs.create('proc_time', proc_time)
                attrs.create('flpops', flpops)
                attrs.create('mflops', mflops)

                print(filename, solver.name(), info, iter, err,
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

            psize = problem.dimension * problem.numberOfContacts

            info = None
            iter = None
            err = None
            time_s = None
            real_time = None
            proc_time = None
            flpops = None
            mflops = None

            digest = hashlib.sha256(open(filename, 'rb').read()).digest()

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
            def pffff(r,v,e):
                solver_problem_callback.get_step(r,v,e)

            if output_errors or output_velocities or output_reactions:
                solver.SolverOptions().callback = pffff
            # get first guess or set guess to zero
            f = h5py.File(filename, 'r')
            if 'guesses' in f:
                number_of_guesses = f['guesses']['number_of_guesses'][0]
                velocities = f['guesses']['1']['u'][:]
                reactions = f['guesses']['1']['r'][:]
            else:
                # guess is missing
                reactions = np.zeros(psize)
                velocities = np.zeros(psize)

            try:
                t0 = time.clock()
                result = self._solver_call(solver,
                                           *(problem, reactions, velocities))
                time_s = time.clock() - t0 # on unix, t is CPU seconds elapsed (floating point)
                
                info, iter, err, real_time, proc_time, flpops, mflops = result

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
            print(filename, solver.name(), info, iter, err,
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

    def name(self):
        return self._name



#
# Some solvers
#
localac = SiconosSolver(name="LocalAlartCurnier",
                        API=N.frictionContact3D_localAlartCurnier,
                        TAG=N.SICONOS_FRICTION_3D_LOCALAC,
                        iparam_iter=1,
                        dparam_err=1,
                        maxiter=maxiter, precision=precision)

localac.SolverOptions().iparam[3] = 10000000


nsgs = SiconosSolver(name="NonsmoothGaussSeidel",
                     API=N.frictionContact3D_nsgs,
                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1,
                     maxiter=maxiter, precision=precision)


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


ExtraGrad = SiconosSolver(name="ExtraGradient",
                          API=N.frictionContact3D_ExtraGradient,
                          TAG=N.SICONOS_FRICTION_3D_EG,
                          iparam_iter=7,
                          dparam_err=1,
                          maxiter=maxiter, precision=precision)

# 1 contact
#Quartic = SiconosSolver(name="Quartic",
#                        API=frictionContact3D_unitary_enumerative,
#                        iparam_iter=7,
#                        dparam_err=1)

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


all_solvers = [nsgs, TrescaFixedPoint, localac, Prox, DeSaxceFixedPoint, ExtraGrad]
if user_solvers == []:
    solvers = all_solvers
else:
    solvers = filter(lambda s : s._name in user_solvers, all_solvers)


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

problem_filenames = filter(is_fclib_file, all_filenames)

n_problems = len(problem_filenames)

problems = [read_fclib_format(f) for f in problem_filenames]



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
    pfilename = os.path.splitext(filename)[0]
    try:
        check_call(['h5copy','-p','-i{0}-{1}.hdf5'.format(solver.name(),pfilename),
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
            self._result_file['data']['comp'][solver.name()][problem_filename]
            return False
        except:
            return True

if __name__ == '__main__':

    if not display:
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
            for solver_name in comp_data:
                filenames = comp_data[solver_name]

                measure[solver_name] = np.empty(n_problems)
                solver_r[solver_name] = np.empty(n_problems)

                ip = 0
                for filename in filenames:
                    if filename not in min_measure:
                        min_measure[filename] = np.inf
                    try:
                        if comp_data[solver_name][filename].attrs['info'] == 0:
                            measure[solver_name][ip] =  comp_data[solver_name][filename].attrs[measure_name]
                            min_measure[filename] = min(min_measure[filename], measure[solver_name][ip])

                            solver_r[solver_name][ip] = measure[solver_name][ip] / \
                              min_measure[filename]
                        else:
                            solver_r[solver_name][ip] = np.inf
                    except:
                        solver_r[solver_name][ip] = np.nan
                    ip += 1

            domain = np.arange(1, 10, .1)
            rhos = dict()
            for solver_name in comp_data:
                rhos[solver_name] = np.empty(len(domain))
                for itau in range(0, len(domain)):
                    rhos[solver_name][itau] = float(len(np.where( solver_r[solver_name] < domain[itau] )[0])) / float(n_problems)

            from matplotlib.pyplot import subplot, title, plot, grid, show, legend, figure

            for solver_name in comp_data:
                plot(domain, rhos[solver_name], label=solver_name)
                legend()
            grid()
        show()

    if display_convergence:
        from matplotlib.pyplot import subplot, title, plot, grid, show, legend, figure
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            for solver_name in comp_data:

                if user_filenames == []:
                    filenames = comp_data[solver_name]
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
        show()
