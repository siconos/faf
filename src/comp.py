#!/usr/bin/env python3
# parallel usage :
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--files={}'
#
# comp.py --max-problems=10 --no-compute --no-collect # output problems.txt
# cat problems.txt | parallel comp.py --timeout=100 --no-collect '--files={}'
#

from __future__ import print_function
import re
from glob import glob
from itertools import product
from subprocess import check_call
import os
import h5py
import getopt
import sys
import hashlib
from mpi4py import MPI
#from io import StringIO


import numpy as np

import siconos.numerics as N
numerics_verbose=0
N.numerics_set_verbose(numerics_verbose)
import siconos.fclib as FCL

#print os.path.join(os.path.dirname(sys.argv[0]), 'external/build')

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), 'external/build'))

try:
    import BogusInterface
except:
    pass
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.INFO)

from SiconosSolver import *
from faf_tools import *
from faf_timeout import *
from faf_matrix_tools import *
from faf_display_tools import *
from faf_default_values import *
from faf_display import *
from faf_preprocess import *
from faf_postprocess import *

#debugger
#import pdb
#pdb.set_trace()


def usage():
  print('Usage: '+sys.argv[0]+'[option]')

  options_doc = """
  Options 
   --help 
     display this message
   --verbose 
     enable verbose mode equal to 1 for Siconos Numerics
   --no-collect 
     leave the result into separate file that are named according the solver and the name of the problem
   --just-collect
     collect all the result into comp.hdf5
   --timeout=n
     set the maximum time of computation for each problem to n seconds (default,utimeout,s)
   --maxiter=n
     set the maximum number of iterations for each problem to n (default",maxiter,")
   --maxiterls=n
     set the maximum number of iterations for each problem to n (default",maxiterls,")
   --domain='a:d:b'
     restrict the domain of the performance profile to the interval [a,b] with a step of d (default",domain[0],":",domain[1]-domain[0],":",domain[-1]+domain[1]-domain[0],")
     or a perfomance profile a should be greater or equal 1
   --measure=value
     select the value  as the measure for the perfomance profile. Possible values are time, iter, flpops
   --display
     perform the computation of performance profile and display it in matplotlib
   --display-distrib or 
     perform the computation of distribution and display it in matplotlib
   --new
     remove comp.hdf5 file
   --solvers=string
     use keyworks in s separated by comma for filtering solvers
   --solvers-exact=string
     use exact names of solvers in s separated by comma for filtering solvers
   --with-mumps
     use mumps as linear system solver
   --max-problems=<max>
     Randomly select <max> problems in current directory.
     The problems list is written in problems.txt file
   --gnuplot
     output gnuplot command file profile.gp or distrib.gp for plotting profiles with gnuplot
   --gnuplot-separate-keys
     output keys anf legend for gnuplot in a separate file.
   --list-contents
     list contents of comp.hdf5 file
   --list-contents-solvers
     list solvers contents of comp.hdf5 file
   --compute-cond-rank
     compute the rank (vairous numerical methods) and condition number of W and store it in the problem file
   --compute-hardness
     compute the average performance of the best solver on a set of problem divided by the average number of contact
   --compute-cond-rank
     compute the conditioning number and the rank of the matrix in the problems

   Other options have to be documented
   
   Usage examples:

   1) running comparison

   comp.py --measure=time --precision=1e-4 --timeout=100 --solvers-exact='NSN-JeanMoreau-NLS','NSN-AlartCurnier-NLS','NSN-NaturalMap-NLS' --no-collect 

   2) collecting results

   comp.py --measure=time --precision=1e-4 --timeout=100 --just-collect

   3) displaying results

   comp.py --display --measure=time  --domain='1:0.1:10'  comp.hdf5

   comp.py --display --measure=time --solvers=Gauss,Tresca,SOCLCP,ACLM --domain=1:0.1:100

   4) For openmp comparison:

   comp.py --thread-list=1,2,3,4,5 --solvers=NSGS-AC-OPENMP

   comp.py --display-speedup --measure=time

  """
  print(options_doc)

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
        print("in get_step")

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

            self._internal_call(solver, sproblem, filename, pfilename,
                                output_filename)
            
        except Exception as e:

            import traceback
            
            print('Exception in internal call')

            traceback.print_exc(e)
            
            try:
               os.remove(output_filename)
            except:
                pass

            with h5py.File(output_filename, 'w') as output:

                digest = hashlib.sha256(open(filename, 'rb').read()).digest()
               
                create_attrs_in_comp_file(output,precision,utimeout,os.uname()[1],measure_name)

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

                attrs.create('filename', np.string_(filename))
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
                if numerics_has_openmp_solvers :
                    try:
                        attrs.create('n_threads', solver.SolverOptions().iparam[10] )
                    except :
                        attrs.create('n_threads',-1)

                list_keys= list(attrs.keys())
                if u'digest' in list_keys:
                    list_keys.remove(u'digest')                
                list_print=[solver.name()]
                list_print.extend([attrs[item] for item in list_keys])
                print(list_print)
            
                with open('report.txt', "a") as report_file:
                    print   (list_print, file=report_file)


#    @timeout(utimeout)
    def _internal_call(self, solver, problem, filename, pfilename, output_filename):

        #print("_internal_call")

        with h5py.File(output_filename, 'w') as output:

            create_attrs_in_comp_file(output,precision,utimeout,os.uname()[1],measure_name)
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
                solver_problem_data.create_dataset(np.string_('reactions'),
                                                   (0, psize),
                                                   maxshape=(None, psize))

                solver_problem_data.create_dataset(np.string_('velocities'),
                                                   (0, psize),
                                                   maxshape=(None, psize))

                solver_problem_data.create_dataset(np.string_('errors'),
                                                   (0, 1),
                                                   maxshape=(None, 1))

            solver_problem_callback = \
              SolverCallback(output, solver_problem_data)

            # need a function, not an instance method for PyObjectCall...
            def pffff(r, v, e):
                solver_problem_callback.get_step(r, v, e)

# callback is broken
#            try:
#                if output_errors or output_velocities or output_reactions:
#                    solver.SolverOptions().callback = pffff
#            except:
#                pass

            # get first guess or set guess to zero
            reactions, velocities, global_velocities = solver.guess(filename)

            normq = np.linalg.norm(problem.q)
            if global_problem:
                _, guess_err = N.gfc3d_compute_error(read_fclib_format(filename)[1],
                                                     reactions, velocities, global_velocities,
                                                     precision, solver.SolverOptions(), normq)

                pass
            else:
                _, guess_err = N.fc3d_compute_error(read_fclib_format(filename)[1],
                                                    reactions, velocities, precision, solver.SolverOptions(), normq)

#            print "guess error:", guess_err

            try:

                if numerics_verbose >1:
                    N.solver_options_print(solver.SolverOptions())

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

                    #t0 = time.process_time()

                    stdout_result = ''
                    t0 = time.time()
                    with catch_stdout(really=output_errors) as get_stdout:
                        if global_problem:
                            result = solver(problem, reactions, velocities, global_velocities)
                        else:
                            result = solver(problem, reactions, velocities)
                        time_s = time.time() - t0 # on unix, t is CPU seconds elapsed (floating point)

                        current_stdout = get_stdout()

                        # try:
                        #     cl = enumerate(filter(lambda s: '||F||' in s,
                        #                           current_stdout.split('\n')))

                        #     rdat = [re.split('=|,', l)[-3:] for i, l in cl]

                        #     dat = [float(r) for r, z, f in rdat]

                        #     for e in dat:
                        #         pffff(None, None, e)

                        # except Exception as e:
                        #     sys.stderr.write('||', type(e))
                            
                        stdout_result += current_stdout

                    
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
                print(exception)
                info = 1
                time_s = np.nan
                iter = np.nan
                err = np.nan
                real_time = np.nan
                proc_time = np.nan
                flpops = np.nan
                mflops = np.nan

            attrs.create('filename', np.string_(filename))
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

            
            if numerics_has_openmp_solvers :
                attrs.create('n_threads', solver.SolverOptions().iparam[10] )


            list_keys= list(attrs.keys())
            if u'digest' in list_keys:
                list_keys.remove(u'digest')                
            list_print=[solver.name()]
            list_print.extend([attrs[item] for item in list_keys])
            print(list_print)
            
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
            with open('report.txt', "a") as report_file:
                print   (list_print, file=report_file)





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
    
def create_attrs_hostname_in_comp_file(comp_file,hostname_val):
    data = comp_file.get('data')
    if data == None :
        data = comp_file.create_group('data')
    comp_data = data.get('comp')
    if comp_data == None:
        comp_data = data.create_group('comp')
    comp_data.attrs.create('hostname',np.string_(hostname_val))

def create_attrs_measure_name_in_comp_file(comp_file,measure_name_val):
    data = comp_file.get('data')
    if data == None :
        data = comp_file.create_group('data')
    comp_data = data.get('comp')
    if comp_data == None:
        comp_data = data.create_group('comp')
    comp_data.attrs.create('measure_name',np.string_(measure_name_val))




def create_attrs_in_comp_file(comp_file,precision_val,utimeout_val,hostname_val,measure_name_val):
    create_attrs_precision_in_comp_file(comp_file,precision_val)
    create_attrs_timeout_in_comp_file(comp_file,utimeout_val)
    create_attrs_hostname_in_comp_file(comp_file,hostname_val)
    create_attrs_measure_name_in_comp_file(comp_file,measure_name_val)

def collect(tpl):

    solver, filename = tpl

    pfilename = os.path.basename(os.path.splitext(filename)[0])
    results_filename = '{0}-{1}.hdf5'.format(solver.name(),pfilename)
    #print "file=", results_filename
    if os.path.exists('comp.hdf5'):
        with h5py.File('comp.hdf5', 'r') as comp_file:
            comp_precision=comp_file['data']['comp'].attrs.get('precision')
            comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            comp_measure_name=comp_file['data']['comp'].attrs.get('measure_name')
            comp_hostname=comp_file['data']['comp'].attrs.get('hostname')
            #print "comp_precision",comp_precision
            if comp_precision == None :
                raise RuntimeError ("Warning. precision information is missing in existing comp.hdf5 file (old version)\n      you must add it with --add-precision-in-comp-file=<val> ")
            if comp_utimeout == None :
                raise RuntimeError ("Warning. timeout information is missing in existing comp.hdf5 file (old version)\n      you must add it with --add-timeout-in-comp-file=<val> ")
            if comp_hostname == None :
                raise RuntimeError ("Warning. hostname information is missing in existing comp.hdf5 file (old version)\n      you must add it with --add-hostname-in-comp-file=<val> ")
    else:
        with h5py.File('comp.hdf5', 'w') as comp_file:
            # if  comp.hdf5 is not existing, we set the attributes precision, timeout, hostanme, measurename to the first file collected.
            with h5py.File( results_filename, 'r+') as result_file:
                    result_precision=result_file['data']['comp'].attrs.get('precision')
                    result_utimeout=result_file['data']['comp'].attrs.get('timeout')
                    result_hostname=result_file['data']['comp'].attrs.get('hostname')
                    result_measure_name=result_file['data']['comp'].attrs.get('measure_name')
            create_attrs_in_comp_file(comp_file,result_precision,result_utimeout,result_hostname,result_measure_name)
            comp_precision=comp_file['data']['comp'].attrs.get('precision')
            comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            comp_measure_name=comp_file['data']['comp'].attrs.get('measure_name')
            comp_hostname=comp_file['data']['comp'].attrs.get('hostname')

    if os.path.exists(results_filename) and not os.stat(results_filename).st_size == 0:
        try:
            if os.path.exists('comp.hdf5'):
                with h5py.File( results_filename, 'r+') as result_file:
                    result_precision=result_file['data']['comp'].attrs.get('precision')
                    result_utimeout=result_file['data']['comp'].attrs.get('timeout')
                    result_hostname=result_file['data']['comp'].attrs.get('hostname')
                    result_measure_name=result_file['data']['comp'].attrs.get('measure_name')
                    if comp_precision != result_precision:
                        raise RuntimeError ("Precision of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_precision,result_precision))
                    if comp_utimeout != result_utimeout:
                        raise RuntimeError ("Timeout of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_utimeout,result_utimeout))
                    if comp_hostname != result_hostname:
                        raise RuntimeError ("hostname of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_hostname,result_hostname))
                    if comp_measure_name != result_measure_name:
                        raise RuntimeError ("Measure of the result in comp.hdf5 ({0}) are not consistent result with the new computed result ({1}) \nWe dot not collect it\nCreate a new comp.hdf5 file".format(comp_measure_name,result_measure_name))

            check_call(['h5copy','-p','-i', results_filename,
                        '-ocomp.hdf5','-s/data/comp/{0}/{1}'.format(solver.name(),pfilename),
                        '-d/data/comp/{0}/{1}'.format(solver.name(),pfilename)])
            if not keep_files:
                os.remove('{0}-{1}.hdf5'.format(solver.name(),pfilename))
        except Exception as e:
            print(e)


class Results():
    def __init__(self, result_file):
        self._result_file = result_file

    def __call__(self, tpl):
        solver = tpl[0]
        problem_filename = os.path.splitext(tpl[1])[0]
        #print("Results: problem_filename:", problem_filename, type(problem_filename) )
        try:
            r = self._result_file['data']['comp'][solver.name()][problem_filename]
            # if abs(r.attrs.get('precision') -  precision) >= 1e-16 :
            #     raise RuntimeError()
            # if abs(r.attrs.get('timeout') -  utimeout) >= 1e-16 :
            #    raise RuntimeError()

            # list_print =[r.attrs['filename'], cond_problem(r.attrs['filename']), solver.name(), r.attrs['info'],
            #                                   r.attrs['iter'], r.attrs['err'], r.attrs['time'], r.attrs['real_time'], r.attrs['proc_time'],
            #                                   r.attrs['flpops'], r.attrs['mflops'],r.attrs.get('precision'),r.attrs.get('timeout')]
            # if numerics_has_openmp_solvers :
            #     list_print.append(r.attrs['n_threads'])
            # print("Already in comp file : ", list_print)


            list_keys= list(r.attrs.keys())
            if u'digest' in list_keys:
                list_keys.remove(u'digest')
            print("Already in comp file : ", [r.attrs[item] for item in list_keys])


            
            return False
        except:
            return True

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    A=N.NM_create(0,1,1)
    N.NM_MPI_set_comm(A, comm)
    N.NM_MUMPS_set_control_params(A)
    N.NM_MUMPS(A, -1)
    if (comm.Get_rank() > 0):
        print('MPI process stop: ',comm.Get_rank())
        exit(0)
    print('MPI process:', comm.Get_rank())

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                       ['help', 'verbose=','no-guess',
                                        'clean', 'display', 'display-convergence','no-matplot',
                                        'files=', 'solvers-exact=', 'solvers=',
                                        'global',
                                        'random-sample=', 'max-problems=',
                                        'timeout=', 'maxiter=', 'maxiterls=', 'precision=',
                                        'keep-files', 'new', 'errors',
                                        'velocities', 'reactions', 'measure=',
                                        'just-collect', 'cond-nc=', 'display-distrib',
                                        'no-collect', 'no-compute', 'domain=',
                                        'replace-solvers-exact=','replace-solvers=',
                                        'gnuplot-output','logscale', 'gnuplot-separate-keys',
                                        'output-dat', 'with-mumps', 'file-filter=', 'remove-files=',
                                        'list-contents','list-contents-solver',
                                        'add-precision-in-comp-file','add-timeout-in-comp-file',
                                        'compute-cond-rank','compute-hardness','test-symmetry','forced','adhoc',
                                        'display-speedup', 'thread-list=','estimate-optimal-timeout'])


    except getopt.GetoptError as err:
            sys.stderr.write('{0}\n'.format(str(err)))
            usage()
            exit(2)
    for o, a in opts:
        if o == '--verbose':
            numerics_verbose=int(a)
            N.numerics_set_verbose(numerics_verbose)
        if o == '--help':
            usage()
            exit(2)
        elif o == '--timeout':
            utimeout = float(a)
        elif o == '--maxiter':
            maxiter = int(a)
        elif o == '--maxiterls':
            maxiterls = int(a)
        elif o == '--global':
           global_problem=True
        elif o == '--precision':
            precision = float(a)
        elif o == '--clean':
            clean = True
        elif o == '--estimate-optimal-timeout':
            compute_rho=True
            estimate_optimal_timeout=True
            compute = False
        elif o == '--display':
            display = True
            compute_rho=True
            compute = False
        elif o == '--list-contents':
            list_contents = True
            compute = False
        elif o == '--list-contents-solver':
            list_contents_solver = True
            compute = False
        elif o == '--display-convergence':
            display_convergence = True
            compute = False
        elif o == '--display-speedup':
            display_speedup = True
            compute = False
        elif o == '--thread-list':
            print(a)
            thread_list =  [int (x) for x in split(a,',')]
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
                    replace_solver_in_compfile = list(filter(lambda s: any(us == s for us in replace_solvers), solver_in_compfile))
                    print("replace solver in comp file", replace_solver_in_compfile)
                    for s in replace_solver_in_compfile:
                        del comp_file['data']['comp'][s]
            except Exception as e:
                print(e)
        elif o == '--replace-solvers':
            replace_solvers = split(a, ',')
            #print "replace_solvers",  replace_solvers
            try:
                with h5py.File('comp.hdf5','r+') as comp_file:
                    solver_in_compfile =  list(comp_file['data']['comp'])
                    #print "list(comp_file['data']['comp'])",  solver_in_compfile
                    replace_solver_in_compfile = list(filter(lambda s: any(us in s for us in replace_solvers), solver_in_compfile))
                    print("replace solver in comp file", replace_solver_in_compfile)
                    for s in replace_solver_in_compfile:
                        del comp_file['data']['comp'][s]
            except Exception as e:
                print(e)
        elif o == '--gnuplot-output':
            gnuplot_output=True
        elif o == '--logscale':
            logscale=True
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

        elif o == '--remove-files':
            remove_file=split(a, ',')


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
        elif o == '--test-symmetry':
            test_symmetry = True
            compute = False
        elif o == '--adhoc':
            adhoc = True
            compute = False
        elif o == '--forced':
            forced=True

    numerics_has_openmp_solvers=False
    try:
        dir(N).index('fc3d_nsgs_openmp')
        numerics_has_openmp_solvers=True
    except ValueError:
        print("warning : fc3d_nsgs_openmp is not in siconos numerics")





    ##################################
    ## creation of solver list
    ##################################
    print("1 -- Creation of solver list")

    if global_problem :
        from faf_global_solvers import *
        fs = faf_global_solvers(maxiter, precision, maxiterls, with_guess, with_mumps, numerics_has_openmp_solvers, mpi_comm=N.NM_MPI_comm(A), mumps_id = N.NM_MUMPS_id(A))
        all_solvers = fs.create_solvers()
    else :
        from faf_solvers import *
        fs = faf_solvers(maxiter, precision, maxiterls, with_guess, with_mumps, numerics_has_openmp_solvers, mpi_comm=N.NM_MPI_comm(A), mumps_id = N.NM_MUMPS_id(A))
        all_solvers = fs.create_solvers()


    if (os.path.isfile(os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'))):
        execfile(os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'))
        #print "execfile(",os.path.join( os.path.dirname(__file__),'adhoc_solverlist.py'), ")"
    if (os.path.isfile('adhoc_solverlist.py')):
        execfile('adhoc_solverlist.py')
        #print "execfile(adhoc_solverlist.py)"

    solvers=[]
    if user_solvers != []:
        #print "user_solvers", user_solvers
        solvers.extend(list(filter(lambda s: any(us in s._name for us in user_solvers), all_solvers)))

        if solvers == []:
            raise RuntimeError ("Cannot find any matching solver")
    elif user_solvers_exact != []:
        #print "user_solvers_exact", user_solvers_exact
        solvers.extend(list(filter(lambda s: any(us ==  s._name  for us in user_solvers_exact), all_solvers)))

        if solvers == []:
            raise RuntimeError("Cannot find any solvers in specified list")
    else:
        solvers= all_solvers

    ##################################
    ## creation of problems list
    ##################################
    print("2 -- Creation of problem list")

    if not os.path.exists('problems.txt'):
         with open('problems.txt', 'w') as problems_txt:
             if global_problem :
                 for f in filter(is_fclib_file_global, glob('*.hdf5')):
                     problems_txt.write('{0}\n'.format(f))
             else:
                 for f in filter(is_fclib_file, glob('*.hdf5')):
                     problems_txt.write('{0}\n'.format(f))

    if user_filenames == []:
        if file_filter == None:
            all_filenames = list_from_file('problems.txt')
        else:
            all_filenames = list(filter(lambda f: any(uf in f for uf in file_filter), list_from_file('problems.txt')))

    else:
            all_filenames = user_filenames

    #all_filenames=['BoxesStack1-i9841-33.hdf5']
    #ask_collect = False
    #print("all_filenames",all_filenames)
    _problem_filenames = list(filter(is_fclib_file,
                                all_filenames))
    #print("_problems_filenames", _problem_filenames)
    __problem_filenames = subsample_problems(_problem_filenames,
                                             random_sample_proba,
                                             max_problems, None, overwrite = (not display and not ask_compute and not ask_collect))


    problem_filenames = subsample_problems(__problem_filenames,
                                           None,
                                           None, cond_nc)

    n_problems = len(problem_filenames)

    #problems = [read_fclib_format(f) for f in problem_filenames]




    # min_measure = dict()

    # for fileproblem in problem_filenames:
    #     min_measure[fileproblem] = np.inf

    if clean:
        h5mode = 'w'
    else:
        h5mode = 'a'


    caller = Caller()


    #pool = MyPool(processes=8)

    rhos=None
    filenames=None
    filename=None





    ### compute ####
    if compute:
        all_tasks = [t for t in product(solvers, problem_filenames)]
        if os.path.exists('comp.hdf5'):
            with h5py.File('comp.hdf5', 'r') as comp_file:
                tasks = list(filter(Results(comp_file), all_tasks))
        else:
            tasks = all_tasks
            
            
        print("3 -- Running computation and/or collecting tasks")
        print("     number of remaining tasks:", len(tasks))
        print("     for solvers :", [ s._name for s in solvers])
        print("     on files ",problem_filenames)
        #print(tasks)


        if ask_compute:
            print(" with precision=", precision, " timeout=", utimeout, "and maxiter = ", maxiter)
            outputs = list(map(caller, tasks))
            N.NM_MUMPS(A, 0)
        if ask_collect:
            list(map(collect, tasks))

    if list_contents or list_contents_solver:
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            for item in comp_data.attrs.keys():
                print("comp_data attrs: ", item + ":", comp_data.attrs[item])
            print("Solvers :")
            for solvername in comp_data:
                print("  ",solvername)
                if (list_contents):
                    for filename in comp_data[solvername]:
                        list_keys= list(comp_data[solvername][filename].attrs.keys())
                        if u'digest' in list_keys:
                            list_keys.remove(u'digest')
                        print("  ",solvername,   [comp_data[solvername][filename].attrs[item] for item in list_keys])

    ### pre-processing ####
    pre = Faf_preprocess('comp.hdf5', problem_filenames)
    if compute_cond_rank:
        pre.compute_cond_rank(forced)
    if test_symmetry:
        pre.test_symmetry()
                        
    ### post-processing ####
                        
    if compute_rho:
        print("3 -- Compute rho for solvers :", [ s._name for s in solvers])
        pp = Faf_postprocess('comp.hdf5', solvers, problem_filenames)
        rhos, solver_r, measure, min_measure, filenames, filename = pp.compute_rho(
            file_filter, remove_file,
            random_sample_proba,max_problems,
            cond_nc, measure_name, domain)
        
    if estimate_optimal_timeout:
        print("4 -- Estimate optimal timeout ")
        pp = Faf_postprocess('comp.hdf5', solvers, problem_filenames)
        pp.estimate_optimal_timeout()
                  
    if compute_hardness:
        # should be in preprocessing
        pp = Faf_postprocess('comp.hdf5', solvers, problem_filenames)
        pp.compute_hardness()

        


    #### display #####

    # print('#### display #####')
    # print(gnuplot_output)
    gnuplot_output=True
    d=Faf_display('comp.hdf5',
              solvers,
              time,
              domain,
              filenames,
              filename,
              rhos,
              gnuplot_output, gnuplot_with_color, gnuplot_separate_keys, no_matplot,logscale)

        
    if display:      
        d.default_display_task(solver_r)

    if display_convergence:
        d.display_convergence(user_filenames,random_sample_proba, max_problems, cond_nc)

 
    if display_distrib:
        print('problem_filenames', problem_filenames)
        d.display_distribution(problem_filenames,gnuplot_output)
        
       
    if display_speedup:
        # this one has to be updated
        d.display_speedup()
 
    display_bw=False
    from matplotlib.pyplot import show
    if display or display_convergence or display_distrib or display_speedup:
        if not no_matplot:
            if (display_bw):
                figs = list(map(figure, get_fignums()))
                for fig in figs:
                    setFigLinesBW(fig)

            show()

         
    if adhoc:
        print("script adhoc (convenient moulinette)")
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
                    print("solver", solver)
                    if ('NSGS-AC-GP' not in solver):
                        print("solver to be renamed", solver)
                        new_solver= solver.replace("AC-","AC-GP-")
                        print("rename", solver, "in ",new_solver)
                        data['comp'].move(solver,new_solver)
