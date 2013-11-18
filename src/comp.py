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
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.INFO)

withpapi = False

if withpapi:
    from ctypes import cdll
    papi=cdll.LoadLibrary('/usr/local/lib/libpapi.so')
    from ctypes import *


    def init_flop():
        ireal_time = c_float()
        iproc_time = c_float()
        iflpops = c_longlong()
        imflops = c_float()
        papi.PAPI_flops(byref(ireal_time), byref(iproc_time), byref(iflpops), 
                    byref(imflops))
        
        def get_flop(real_time, proc_time, flpops, mflops):
            r = papi.PAPI_flops(byref(real_time), byref(proc_time), byref(flpops), 
                                byref(mflops))
            

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

    def __call__(self,problem,reactions,velocities):
        if withpapi:
            real_time = c_float()
            proc_time = c_float()
            flpops = c_longlong()
            mflops = c_float()
            init_flop()
            info = self._API(problem, reactions, velocities, self._SO)
            get_flop(real_time, proc_time, flpops, mflops)
            return (info, self._get(self._SO.iparam, self._iparam_iter), self._get(self._SO.dparam, self._dparam_err), real_time.value, proc_time.value, flpops.value, mflops.value)
        else :
            info = self._API(problem,reactions,velocities,self._SO)            
            return (info, self._get(self._SO.iparam, self._iparam_iter), self._get(self._SO.dparam, self._dparam_err), 0.0, 0.0, 0.0, 0.0)
            

    def name(self):
        return self._name



#
# Some solvers
#

localac = SiconosSolver(name="local Alart Curnier",
                         API=frictionContact3D_localAlartCurnier,
                         TAG=SICONOS_FRICTION_3D_LOCALAC,
                         iparam_iter=1,
                         dparam_err=1)

nsgs = SiconosSolver(name="nonsmooth Gauss Seidel",
                     API=frictionContact3D_nsgs,
                     TAG=SICONOS_FRICTION_3D_NSGS,
                     iparam_iter=7,
                     dparam_err=1)


# only dense
nsgsv = SiconosSolver(name="nonsmooth Gauss Seidel velocity",
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
                          iparam_iter=7,
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
# http://code.google.com/p/mpi4py/
#from mpi4py import MPI
#frictionContact3D_sparseGlobalAlartCurnierInit(localac.SolverOptions())

setNumericsVerbose(0)

solvers = [nsgs, localac, TrescaFixedPoint, Prox, DeSaxceFixedPoint, ExtraGrad]


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
                return f, read_fclib_format(f)
    except Exception as e:
        print e
        return None, None


fileproblems = imap(read_problem,glob("*.hdf5"))

for fileproblem,solver in product(fileproblems,solvers):

    if fileproblem[0] is not None:

        problem = fileproblem[1]
        filename = fileproblem[0]

        f = h5py.File(filename, 'r')
        if 'guesses' in f:
            number_of_guesses=f['guesses']['number_of_guesses'][0]
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
            iter = np.inf
            err = np.inf
            real_time = np.inf
            proc_time = np.inf
            flpops = np.inf
            mflops = np.inf
            
            # need output in a csv database

        # filename, solver name, revision svn, parameters, nb iter, err
        print(filename, solver.name(), info, iter, err, real_time, proc_time, flpops, mflops)
