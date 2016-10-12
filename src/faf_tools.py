import siconos.fclib as FCL
import siconos.numerics as N
import numpy as np
import h5py
from ctypes import cdll, c_float, c_longlong, byref
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
                print 'Exception in Memoize', e 
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
                print 'Exception in _numberOfDegreeofFreedom', e
                r = np.nan
                
            #print "r=",r
    return r

def _numberOfDegreeofFreedomContacts(f):
    with h5py.File(f, 'r') as fclib_file:
        try:
            r = fclib_file['fclib_local']['W']['m'][0]
        except Exception as e:
            print 'Exception in _numberOfDegreeofFreedomContacts', e
            r = np.nan
    #print "r=",r
    return r


numberOfDegreeofFreedom = Memoize(_numberOfDegreeofFreedom)

numberOfDegreeofFreedomContacts = Memoize(_numberOfDegreeofFreedomContacts)
