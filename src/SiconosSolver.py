import siconos.numerics as N
import h5py
import numpy as np
from faf_tools import *
from faf_papi import *

class SiconosSolver():
    _name = None
    _gnuplot_name = None
    _API = None
    _TAG = None
    _iparam_iter = None
    _dparam_err = None
    _SO = None
    _mpi_comm = None
    _mumps_id = None

    def _get(self, tab, index):
        if index is not None:
            return tab[index]
        else:
            return None

    def __init__(self, name=None, gnuplot_name=None, API=None, TAG=None, iparam_iter=None,
                 dparam_err=None, maxiter=None, precision=None, with_guess=None, global_solver=False, mpi_comm=None, mumps_id=None):
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
        if (maxiter == None):
            raise RuntimeError("SiconosSolver() maxiter have to be specified.")
        else:
            self._SO.iparam[0] = maxiter
        if (precision==None):
            raise RuntimeError("SiconosSolver() precision have to be specified.")
        else:
            self._SO.dparam[0] = precision
        if (with_guess==None):
            raise RuntimeError("SiconosSolver() with_guess have to be specified.")
        else:
            self._with_guess = with_guess
        self._global_solver = global_solver
        self._mpi_comm = mpi_comm
        self._mumps_id = mumps_id
        
        
    def SolverOptions(self):
        return self._SO

    def __call__(self, problem, reactions, velocities, global_velocities= None):
#        N.frictionContact_display(problem)

        if self._mpi_comm is not None:
            N.NM_MPI_set_comm(problem.M, self._mpi_comm)
            print ('SET MPI:', problem, problem.M, self._mpi_comm)
            
        if self._mumps_id is not None:
            print('MUMPS_id:', self._mumps_id)
            N.NM_MUMPS_set_id(problem.M, self._mumps_id)
            print ('SET MUMPS id:', problem, problem.M, self._mumps_id)
            print ('MUMPS id icntl 14:', N.NM_MUMPS_icntl(problem.M, 14))

#        N.frictionContact_display(problem)
        real_time = c_float()
        proc_time = c_float()
        flpops = c_longlong()
        mflops = c_float()
        init_flop()

        if self._global_solver :
            info = self._API(problem, reactions, velocities, global_velocities, self._SO)
            pass
        else:
            info = self._API(problem, reactions, velocities, self._SO)

        get_flop(real_time, proc_time, flpops, mflops)
        return (info, self._get(self._SO.iparam, self._iparam_iter),
                self._get(self._SO.dparam, self._dparam_err),
                real_time.value, proc_time.value,
                flpops.value, mflops.value)

    def guess(self, filename):
        problem = read_fclib_format(filename)[1]

        with h5py.File(filename, 'r') as f:
            if 'fclib_global' in f:
                tag = 'fclib_global'
            else:
                tag = 'fclib_local'
            psize = numberOfDegreeofFreedomContacts(filename)
            nsize = numberOfDegreeofFreedom(filename)
            global_velocities = None
            if self._with_guess and 'guesses' in f:
                number_of_guesses = f['guesses']['number_of_guesses'][0]
                velocities = f['guesses']['1']['u'][:]
                reactions = f['guesses']['1']['r'][:]
                if tag == 'fclib_global':
                    global_velocities =f['guesses']['1']['v'][:]
                
            else:
                # guess is missing
                reactions = np.zeros(psize)
                velocities = np.zeros(psize)
                if tag == 'fclib_global':
                    global_velocities = np.zeros(nsize)

        return reactions, velocities, global_velocities

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




class SiconosHybridSolver(SiconosSolver):

    def guess(self, filename):
        pfilename = os.path.splitext(filename)[0]
        with h5py.File('comp.hdf5', 'r') as comp_file:
            return extern_guess(pfilename, 'NSGS', 1, comp_file)


class SiconosWrappedSolver(SiconosSolver):
    def __call__(self, problem, reactions, velocities):
        return self._API(problem, reactions, velocities, self._SO)
