#!/usr/bin/env python

from glob import glob
from itertools import product, chain, imap
import numpy as np
from Siconos.Numerics import *
import Siconos.FCLib as FCL

import os

import multiprocessing
import time
import logging

import h5py
import getopt
import sys
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.INFO)

measure = 'flop'

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'flop', 'iter'])
except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)
for o, a in opts:
    if o == '--flop':
        measure = 'flop'
    elif o == '--iter':
        measure = 'iter'


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
        flpops.value = np.nan
        mflops.value = np.nan

real_time = c_float()
proc_time = c_float()
flpops = c_longlong()
mflops = c_float()
#init_flop()
#a=2.
#b=3.
#c=a+b
#get_flop(real_time, proc_time, flpops, mflops)
#papi.PAPI_stop_counters()
#print real_time.value, proc_time.value, flpops.value, mflops.value
#papi.PAPI_start_counters()
#init_flop()
#a=1.
#b=2.
#c=a+b
#get_flop(real_time, proc_time, flpops, mflops)
#print real_time.value, proc_time.value, flpops.value, mflops.value



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


flops = dict()

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
                return result
            else:
                raise result
        return inner
    return wrapper


#pool = multiprocessing.Pool(processes=1)

class SiconosSolver():
    _name = None
    _API = None
    _TAG = None
    _iparam_iter = None
    _dparam_err = None
    _SO = None

    def _get(self,tab,index):
        if index is not None:
            return tab[index]
        else:
            return None

    def __init__(self, name = None, API = None, TAG = None, iparam_iter = None, dparam_err = None):
        self._name = name
        self._API = API
        self._TAG = TAG
        self._iparam_iter = iparam_iter
        self._dparam_err = dparam_err
        self._SO = SolverOptions(TAG) # set default solver options
        self._SO.iparam[0] = 100000
        self._SO.dparam[0] = 1e-8

    def SolverOptions(self):
        return self._SO

    @timeout(100)
    def __call__(self,problem,reactions,velocities):
        real_time = c_float()
        proc_time = c_float()
        flpops = c_longlong()
        mflops = c_float()
        init_flop()
        info = self._API(problem, reactions, velocities, self._SO)
        get_flop(real_time, proc_time, flpops, mflops)

        return (info, self._get(self._SO.iparam, self._iparam_iter), self._get(self._SO.dparam, self._dparam_err), real_time.value, proc_time.value, flpops.value, mflops.value)

    def name(self):
        return self._name



#
# Some solvers
#

localac = SiconosSolver(name="Local Alart Curnier",
                         API=frictionContact3D_localAlartCurnier,
                         TAG=SICONOS_FRICTION_3D_LOCALAC,
                         iparam_iter=1,
                         dparam_err=1)

localac.SolverOptions().iparam[3] = 10000000


nsgs = SiconosSolver(name="Nonsmooth Gauss Seidel",
                     API=frictionContact3D_nsgs,
                     TAG=SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1)


# only dense
nsgsv = SiconosSolver(name="Nonsmooth Gauss Seidel velocity",
                      API=frictionContact3D_nsgs_velocity,
                      TAG=SICONOS_FRICTION_3D_NSGSV,
                      iparam_iter=7,
                      dparam_err=1)

TrescaFixedPoint = SiconosSolver(name="Tresca Fixed Point",
                                 API=frictionContact3D_TrescaFixedPoint,
                                 TAG=SICONOS_FRICTION_3D_TFP,
                                 iparam_iter=7,
                                 dparam_err=1)

DeSaxceFixedPoint = SiconosSolver(name="DeSaxce Fixed Point",
                                  API=frictionContact3D_DeSaxceFixedPoint,
                                  TAG=SICONOS_FRICTION_3D_DSFP,
                                  iparam_iter=7,
                                  dparam_err=1)

Prox = SiconosSolver(name="Proximal fixed point",
                     API=frictionContact3D_proximal,
                     TAG=SICONOS_FRICTION_3D_PROX,
                     iparam_iter=7,
                     dparam_err=1)


ExtraGrad = SiconosSolver(name="Extra gradient",
                          API=frictionContact3D_ExtraGradient,
                          TAG=SICONOS_FRICTION_3D_EG,
                          iparam_iter=1,
                          dparam_err=1)

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

hyperplaneProjection = SiconosSolver(name="hyperplane projection",
                                     API=frictionContact3D_HyperplaneProjection,
                                     TAG=SICONOS_FRICTION_3D_HP,
                                     iparam_iter=7,
                                     dparam_err=1)


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
#frictionContact3D_sparseGlobalAlartCurnierInit(localac.SolverOptions())

setNumericsVerbose(0)

solvers = [nsgs, TrescaFixedPoint, localac, Prox, DeSaxceFixedPoint, ExtraGrad]


def is_fclib_file(filename):
    with h5py.File(filename, 'r') as f:
        return 'fclib_local' in f or 'fclib_global' in f

def read_numerics_format(f):
    return frictionContactProblemFromFile(f)

def read_fclib_format(f):
    return frictionContact_fclib_read(f)

def read_problem(f):
    try:
        ext = os.path.splitext(f)[1]
        if ext == ".dat":
            return f, read_numerics_format(f)
        else:
            if ext == ".hdf5":
                if is_fclib_file(f):
                    return f, read_fclib_format(f)
                else:
                    raise Exception
    except Exception as e:
        print e
        return None, None


fileproblems = imap(read_problem,glob("*.hdf5"))

rfileproblems = [ f for f in fileproblems ]

solver_flpops = dict()
solver_r = dict()
n_problems = len(rfileproblems)

for solver in solvers:
    solver_flpops[solver] = np.empty(n_problems)
    solver_r[solver] = np.empty(n_problems)

min_flpops = dict()

for fileproblem in rfileproblems:
    min_flpops[fileproblem] = np.inf         

for solver in solvers:
    ip = 0
    for fileproblem in rfileproblems:

        if fileproblem[0] is not None:

            problem = fileproblem[1]
            filename = fileproblem[0]

            f = h5py.File(filename, 'r')
            if 'guesses' in f:
                number_of_guesses = f['guesses']['number_of_guesses'][0]
                velocities = f['guesses']['1']['u']
                reactions = f['guesses']['1']['r']
            else:
                # guess is missing
                reactions = np.zeros(problem.dimension * problem.numberOfContacts)
                velocities = np.zeros(problem.dimension * problem.numberOfContacts)

            try:
                info, iter, err, real_time, proc_time, flpops, mflops = \
                    solver(problem, reactions, velocities)
            except Exception as exception:
                print exception
                info = 1
                iter = np.nan
                err = np.nan
                real_time = np.nan
                proc_time = np.nan
                flpops = np.nan
                mflops = np.nan

                # need output in a csv database

            # filename, solver name, revision svn, parameters, nb iter, err
            print(filename, solver.name(), info, iter, err, real_time, proc_time, 
                  flpops, mflops)

            if measure == 'flop':
                solver_flpops[solver][ip] = flpops
            elif measure == 'iter':
                solver_flpops[solver][ip] = iter
                
            min_flpops[fileproblem] = min(flpops, min_flpops[fileproblem])
            ip += 1
            

for solver in solvers:
    ip = 0
    for fileproblem in rfileproblems:
        solver_r[solver][ip] = solver_flpops[solver][ip] / min_flpops[fileproblem]
        ip += 1

domain = np.arange(1,10,.1)

rhos = dict()
for solver in solvers:
    rhos[solver] = np.empty(len(domain))
    for itau in range(0, len(domain)):
        rhos[solver][itau] = float(len(np.where( solver_r[solver] < domain[itau] )[0])) / float(n_problems)


from matplotlib.pyplot import subplot, title, plot, grid, show, legend

for solver in solvers:
    print solver.name(), rhos[solver]
    plot(domain, rhos[solver], label=solver.name())
    legend()
grid()
show()
