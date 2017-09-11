import siconos.fclib as FCL
import siconos.numerics as N
import numpy as np
import h5py
from ctypes import cdll, c_float, c_longlong, byref

import shlex

def is_fclib_file(filename):
    r = False
    try:
        with h5py.File(filename, 'r') as f:
            r = 'fclib_local' in f or 'fclib_global' in f
    #except Exception as e:
    #    print(e)
    except :
        pass
    return r

def read_numerics_format(f):
    return N.frictionContactProblemFromFile(f)



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

    # filename must be of string type to be read correctly
    # if filename from an attribute in hdf5 file that have been encoded with np.string_
    # it has to be decoded with np.decode
    # print("_read_fclib_format, filename:", filename, type(filename))
    # print("_read_fclib_format, filename.decode():", filename.decode(), type(filename.decode()))
    # we do not prefer to force the conversion there
    
    fclib_problem = FCL.fclib_read_local(filename)

    numerics_problem =  N.from_fclib_local(fclib_problem)
    return fclib_problem, numerics_problem

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
                print('Exception in Memoize', e )
                self._done[args] = e
                return e


read_fclib_format = Memoize(_read_fclib_format)

def _numberOfDegreeofFreedom(f):
    with h5py.File(f, 'r') as fclib_file:

        try:
            r = 6*fclib_file['fclib_local']['info'].attrs['numberOfInvolvedDS']
        except:
            try:
                r = fclib_file['fclib_local']['info'].attrs['numberOfDegreeOfFreedom'][0]
            except Exception as e:
                print('Exception in _numberOfDegreeofFreedom', e)
                r = np.nan
                
            #print "r=",r
    return r

def _numberOfDegreeofFreedomContacts(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W']['m'][0]
        except Exception as e:
            print('Exception in _numberOfDegreeofFreedomContacts', e)
            r = np.nan
    #print "r=",r
    return r


numberOfDegreeofFreedom = Memoize(_numberOfDegreeofFreedom)

numberOfDegreeofFreedomContacts = Memoize(_numberOfDegreeofFreedomContacts)




class WithCriterium():

    def __init__(self, condmin, condmax):
        self._condmin = condmin
        self._condmax = condmax

    def __call__(self, filename):
        r = cond_problem(filename)
        return r > self._condmin and r < self._condmax

import random

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
        print(_r)

    else:
        _r = __r

    if cond is not None:
        r = list(filter(WithCriterium(cond[0], cond[1]), _r))
    else:
        r = _r

    # overwrite
    if overwrite:
        with open('problems.txt','w') as problems_txt:
            problems_txt.write('{0}\n'.format('\n'.join(r)))

    return r



def list_from_file(filename):
    with open(filename, 'r') as f:
        return f.read().lstrip().rstrip().split('\n')



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


from contextlib import contextmanager
import select

# to collect output from stdout
# note : callback is much cleaner (no parsing) but is broken...
@contextmanager
def catch_stdout(really=True):
    if really:
        sys.stdout.write(' \b')
        pipe_out, pipe_in = os.pipe()

        def read_pipe():
            def more():
                r, _, _ = select.select([pipe_out], [], [], 0)
                return bool(r)
            out = ''
            while more():
                out += os.read(pipe_out, 1024)
            return out

        stdout = os.dup(1)
        os.dup2(pipe_in, 1)

        yield read_pipe

        os.dup2(stdout, 1)
    else:
        def dummy():
            return ''
        yield dummy
        pass

def _dimension(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['spacedim'][0]
        except:
            r = np.nan
    return r

dimension = Memoize(_dimension)




import os
import platform

def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file)
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return stat.st_mtime



