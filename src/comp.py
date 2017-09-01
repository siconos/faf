#!/usr/bin/env python
# parallel usage :
# ls *.hdf5 | parallel comp.py --timeout=100 --no-collect '--files={}'
#
# comp.py --max-problems=10 --no-compute --no-collect # output problems.txt
# cat problems.txt | parallel comp.py --timeout=100 --no-collect '--files={}'
#


import re
from glob import glob
from itertools import product
from subprocess import check_call
import os
import h5py
import getopt
import sys
import hashlib

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
   --display-distrib='from-files' or 
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
   --gnuplot-profile
     output gnuplot command file profile.gp for plotting profiles woth gnuplot
   --gnuplot-distrib
     output gnuplot command file distrib.gp for plotting distribution woth gnuplot
   --gnuplot-separate-keys
     output keys anf legend for gnuplot in a separate file.
   --display-distrib='from-files' 
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


try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'verbose=','no-guess',
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
                                    'output-dat', 'with-mumps', 'file-filter=', 'remove-files=',
                                    'list-contents','list-contents-solver',
                                    'add-precision-in-comp-file','add-timeout-in-comp-file',
                                    'compute-cond-rank','compute-hardness','adhoc',
                                    'display-speedup', 'thread-list='])


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
    elif o == '--adhoc':
        adhoc = True
        compute = False

numerics_has_openmp_solvers=False
try:
    dir(N).index('fc3d_nsgs_openmp')
    numerics_has_openmp_solvers=True
except ValueError:
    print("warning : fc3d_nsgs_openmp is not in siconos numerics")



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

            print('Exception in internal call', e)

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

                list_print = [filename, numberOfDegreeofFreedomContacts(filename), numberOfDegreeofFreedom(filename), cond_problem(filename), solver.name(), info, iter, err,
                              time_s, real_time, proc_time,
                              flpops, mflops,
                              precision, utimeout]

                if numerics_has_openmp_solvers :
                    try:
                        attrs.create('n_threads', solver.SolverOptions().iparam[10] )
                        list_print.append(solver.SolverOptions().iparam[10])
                    except :
                        attrs.create('n_threads',-1)
                        list_print.append(-1)

                print(list_print)

                with open('report.txt', "a") as report_file:
                    print   ( (filename, solver.name(), info, iter, err,
                                   time_s, real_time, proc_time,
                                   flpops, mflops,
                                   precision, utimeout), file=report_file)




    @timeout(utimeout)
    def _internal_call(self, solver, problem, filename, pfilename, output_filename):


        print("_internal_call")

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
            reactions, velocities = solver.guess(filename)

            normq = np.linalg.norm(problem.q)
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

                    t0 = time.time()
                    #t0 = time.process_time()

                    stdout_result = ''

                    with catch_stdout(really=output_errors) as get_stdout:

                        result = solver(problem, reactions, velocities)

                        current_stdout = get_stdout()

                        try:

                            cl = enumerate(filter(lambda s: '||F||' in s,
                                                  current_stdout.split('\n')))

                            rdat = [re.split('=|,', l)[-3:] for i, l in cl]

                            dat = [float(r) for r, z, f in rdat]

                            for e in dat:
                                pffff(None, None, e)

                        except Exception as e:
                            sys.stderr.write('||', type(e))
                            
                        stdout_result += current_stdout

                    
                    time_s = time.time() - t0 # on unix, t is CPU seconds elapsed (floating point)
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
            # filename, solver name, revision svn, parameters, nb iter, err

            list_print = [filename, numberOfDegreeofFreedomContacts(filename), numberOfDegreeofFreedom(filename), cond_problem(filename), solver.name(), info, iter, err,
                          time_s, real_time, proc_time,
                          flpops, mflops,
                          precision, utimeout]

            if numerics_has_openmp_solvers :
                attrs.create('n_threads', solver.SolverOptions().iparam[10] )
                list_print.append(solver.SolverOptions().iparam[10])

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
                print   ( (filename, solver.name(), info, iter, err,
                               time_s, real_time, proc_time,
                               flpops, mflops,
                               precision, utimeout), file=report_file)





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
            comp_measure_name=comp_file['data']['comp'].attrs.get('measure_name')
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
            comp_measure_name=comp_file['data']['comp'].attrs.get('measure_name')

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


##################################
## creation of solver list
##################################
print("1 -- Creation of solver list")
from faf_solvers import *
fs = faf_solvers(maxiter, precision, maxiterls, with_guess, with_mumps, numerics_has_openmp_solvers)
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



measure = dict()
solver_r = dict()
min_measure = dict()

for fileproblem in problem_filenames:
    min_measure[fileproblem] = np.inf

if clean:
    h5mode = 'w'
else:
    h5mode = 'a'

caller = Caller()


#pool = MyPool(processes=8)




if __name__ == '__main__':

    if compute:
        all_tasks = [t for t in product(solvers, problem_filenames)]
        if os.path.exists('comp.hdf5'):
            with h5py.File('comp.hdf5', 'r') as comp_file:
                tasks = list(filter(Results(comp_file), all_tasks))
        else:
            tasks = all_tasks

        print("3 -- Running computation and/or collecting tasks for solvers :", [ s._name for s in solvers])
        print(" on files ",problem_filenames)
            
        if ask_compute:
            print(" with precision=", precision, " timeout=", utimeout, "and maxiter = ", maxiter)
            outputs = list(map(caller, tasks))
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



    if display:
        print("3 -- Running display tasks for solvers :", [ s._name for s in solvers])
        filename=None
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            comp_precision=comp_file['data']['comp'].attrs.get('precision')
            comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            # 1 n_problems
            n_problems = 0

            for solver in solvers:
                solver_name=solver.name()
                if solver_name in comp_data :
                    if file_filter == None:
                        all_filenames = comp_data[solver_name]
                    else:
                        all_filenames = list(filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name]))

                    if remove_file != None:
                        remove_file_without_ext=[]
                        for uf in remove_file:
                            remove_file_without_ext.append(uf.split('.')[:-1][0])
                        all_filenames = list(filter(lambda f: any(uf not in f for uf in remove_file_without_ext), all_filenames))

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
                        all_filenames = list(filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name]))

                    if remove_file != None:
                        remove_file_without_ext=[]
                        for uf in remove_file:
                            remove_file_without_ext.append(uf.split('.')[:-1][0])
                        all_filenames = list(filter(lambda f: any(uf not in f for uf in remove_file_without_ext), all_filenames))


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
                        all_filenames = list(filter(lambda f: any(uf in f for uf in file_filter), comp_data[solver_name]))
                    if remove_file != None:
                        remove_file_without_ext=[]
                        for uf in remove_file:
                            remove_file_without_ext.append(uf.split('.')[:-1][0])
                        all_filenames = list(filter(lambda f: any(uf not in f for uf in remove_file_without_ext), all_filenames))



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
                            print(line, file=input_file)

                out_data=np.empty([len(domain),len(comp_data)+1])
                write_report(rhos,'rhos.txt')
                write_report(solver_r,'solver_r.txt')
                def long_substr(data):
                    substr = ''
                    #print("data=",data)
                    if len(data) > 0 and len(data[0]) > 0:
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
                    #print("filenames=",filenames)
                    #print("long_substr(filenames)=", long_substr(filenames))
                    test_name = long_substr(filenames).partition('-')[0]
                    print("test_name=",test_name)

                    if test_name.endswith('_'):
                        test_name  = test_name[:-1]
                    gp.write('basename="profile-{0}"\n'.format(test_name))
                    #print filename.partition('-')[0]
                    print("test_name=",test_name)
                    gp.write('\n')
                    gp.write('term_choice_tikz=1\n')
                    gp.write('if (term_choice_tikz == 1) \\\n')
                    if (gnuplot_with_color):
                        gp.write('set term tikz standalone size 5in,3in font \'\\scriptsize\\sf\';  \\\n')
                    else:
                        gp.write('set term tikz standalone monochrome  size 5in,3in font \'\\scriptsize\\sf\';  \\\n')
                    gp.write('extension = \'.tex\'; \\\n')
                    gp.write('extension_legend = \'_legend.tex\'; \\\n')
                    gp.write('set output basename.extension; \\\n')
                    gp.write('print "output = ", basename.extension; \\\n')

                    gp.write('else \\\n')
                    gp.write('set term aqua;\\\n')
                    gp.write('\n')
                    gp.write('set title\'{0} - precision: {1} - timeout: {2} \';; \n'.format(test_name,comp_precision,comp_utimeout))


                    
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
                            gp.write('unset title; \n \n')
                            gp.write('unset tics; \n \n')
                            gp.write('unset xlabel; \n \n')
                            gp.write('unset ylabel; \n \n')
                            gp.write('set term tikz standalone  size 5in,1.5in font \'\\scriptsize\\sf\';  \\\n')
                            gp.write('set key right inside vertical maxrows {0}\n'.format(maxrows))
                            gp.write('\n plot [0:1] [0:1]')
                            gp.write(
                                ','.join([' NaN t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))



                        else:
                            gp.write(
                                ','.join(['resultfile using 1:{0} notitle w l dashtype {1} linecolor {2}'.format(index + 2,index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))

                            gp.write('\n set output basename.extension_legend; \n')
                            gp.write('print "output = ", basename.extension_legend; \n \n')
                            gp.write('unset border; \n \n')
                            gp.write('unset title; \n \n')
                            gp.write('unset tics; \n \n')
                            gp.write('unset xlabel; \n \n')
                            gp.write('unset ylabel; \n \n')
                            gp.write('set term tikz standalone monochrome  size 5in,1.5in font \'\\scriptsize\\sf\';  \\\n')
                            gp.write('set key right inside vertical maxrows {0}\n'.format(maxrows))
                            gp.write('\n plot [0:1] [0:1]')
                            gp.write(
                                ','.join([' NaN t "{1}" w l dashtype {2} linecolor {3}'.format(index + 2, solver.gnuplot_name(),index+1,8)
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
                print("Warning: no problem corresponding to the required solver")
                if (os.path.isfile('profile.gp')):
                    os.remove('profile.gp')

            if not no_matplot:
                # 5 plot
                from matplotlib.pyplot import subplot, title, plot, grid, show, get_fignums, legend, figure, xlim, ylim, xscale

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
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import show
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']

            if user_filenames == []:
                for solver_name in comp_data:

                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
            else:
                filenames = user_filenames
            
            for filename in filenames:

                fig, axs = plt.subplots(1, 1)
                ax = axs

                ax.set_title('Convergence on {0}'.format(filename))
                ax.grid(True, which="both")

                ax.set_yscale('symlog', linthreshy=0.001)
                
                for solver_name in comp_data:

                    try:
                        pfilename = os.path.splitext(filename)[0]
                        solver_problem_data = comp_data[solver_name][pfilename]
                    
                        ax.plot(np.arange(len(solver_problem_data['errors'][:])),
                                np.log(solver_problem_data['errors']),
                                label='{0}'.format(solver_name))
                        ax.legend(loc='lower left')
                        
                    except:
                        pass



    if compute_cond_rank:
        print("Tasks will be run for", problem_filenames)
        for problem_filename in problem_filenames:
            print("compute for", problem_filename,"....")
            with h5py.File(problem_filename, 'r+') as fclib_file:
                no_rank_info=True
                if (fclib_file['fclib_local']['W'].attrs.get('rank') == None) :
                    print("Rank info already not  in", problem_filename)
                else:
                    print("Rank info already in", problem_filename)
                    no_rank_info=False
            if no_rank_info:
                try:
                    [norm_lsmr, cond_lsmr, max_nz_sv, min_nz_sv, cond, rank, rank_dense, rank_svd, rank_estimate] = norm_cond(problem_filename)
                    print( problem_filename, norm_lsmr, cond_lsmr, max_nz_sv, min_nz_sv, cond,  rank_dense, rank_svd, rank_estimate)
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
                    print("-->", e)

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

    if compute_hardness:
        nc = []
        nds = []
        cond_nc = []
        max_measure = dict()

        for fileproblem in problem_filenames:
            max_measure[fileproblem] = - np.inf

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
        #print(nc)
        nc_avg = sum(nc)/float(len(nc))
        print("nc_avg", nc_avg)
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
                                max_measure[filename] = max(max_measure[filename], measure[solver_name][ip])
                            else:
                                measure[solver_name][ip] = np.inf
                        except:
                            measure[solver_name][ip] = np.nan
                        ip += 1

            #print("min_measure", min_measure)
            #print("max_measure", max_measure)

            min_measure_array=np.array([min_measure[key] for key in min_measure.keys()])
            avg_min_measure = min_measure_array.mean()
            std_min_measure = min_measure_array.std()                
            print(         "Average min resolution measure (avg fastest solver measure) = {0:12.8e}".format(avg_min_measure))
            print(         "Std min resolution measure (std fastest solver measure) = {0:12.8e}".format(std_min_measure))
            print(         "Average min resolution measure by contact = {0:12.8e}".format(avg_min_measure/nc_avg))
            

            max_measure_array=np.array([max_measure[key] for key in max_measure.keys()])
            avg_max_measure = max_measure_array.mean()
            std_max_measure = max_measure_array.std()    
            print(         "Average max resolution measure (avg slowest suceeded solver measure) = {0:12.8e}".format(avg_max_measure))
            print(         "Std max resolution measure (std fastest solver measure) = {0:12.8e}".format(std_max_measure))
            print(         "Average max resolution measure by contact = {0:12.8e}".format(avg_max_measure/nc_avg))



    if display_distrib:
        from matplotlib.pyplot import title, subplot, grid, show, get_fignums, legend, figure, hist, xlim, ylim, xscale
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
            print("nds", nds)
            print("rank_dense_W",  rank_dense_W)
            print("cond_nc", cond_nc)
            print("cond_W", cond_W)


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
                    gp.write('print \'max_nc =\', max_nc,\' min_nc =\', min_nc \n')
                    gp.write('binwidth = (max_nc-min_nc)/numberofbox\n')
                    gp.write('set boxwidth binwidth\n')

                    gp.write('print \'binwidth =\', binwidth \n')

                    gp.write('set xlabel \'number of contacts\' offset 0,1.2 \n')
                    #gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes title \'number of contacts\'  \n')
                    gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes notitle  \n')
                    gp.write('\n')
                    gp.write('print \'max_ndof =\', max_ndof,\' min_ndof =\', min_ndof \n')
                    gp.write('if ( (max_ndof-min_ndof) < numberofbox) {binwidth =1 ;} else {binwidth = (max_ndof-min_ndof)/numberofbox}\n');

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
                                    print(filename)
                                    x[filename + solver_name] = cond_problem(filename)
                                else:
                                    x[filename + solver_name] = np.nan
                figure()
                l = [x[k] for k in x]
                l.sort()
                values = array(l)
                hist(values, 100, range=(min(values), max(values)), histtype='stepfilled')
                grid()


    if display_speedup:
        from matplotlib.pyplot import subplot, title, plot, grid, show, get_fignums, legend, figure, hist, bar, xlabel, ylabel, boxplot
        print('\n display speedup is starting ...')
        with h5py.File('comp.hdf5', 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            result_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            print('solver in comp_data =',[s for s in comp_data] )
            solvers=[]
            if user_solvers != []:
                #print "user_solvers", user_solvers
                solvers.extend( list(filter(lambda s: any(us in s for us in user_solvers), comp_data)))

                if solvers == []:
                    raise RuntimeError ("Cannot find any matching solver")

            elif user_solvers_exact != []:
                #print "user_solvers_exact", user_solvers_exact
                solvers.extend(list(filter(lambda s: any(us ==  s  for us in user_solvers_exact), comp_data)))

                if solvers == []:
                    raise RuntimeError("Cannot find any solvers in specified list")

            else:
                solvers= comp_data

            print(' filtered solver in comp_data =',[s for s in solvers])

            solvers= list(filter(lambda s: ('OPENMP' in s), solvers))


            print(' solver for speedup-display =',[s for s in solvers])
            if solvers == []:
                raise RuntimeError("Cannot find any solvers in specified list")
            #---#
            # collect results by filename
            #---#
            results_filename={}
            nthread_set= set()
            n_filename=[]

            for solver_name in solvers:
                if user_filenames == []:
                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
                else:
                    filenames = user_filenames
                n_filename.append(len(filenames))
                for filename in filenames:
                    pfilename = os.path.splitext(filename)[0]
                    nthread=int(solver_name.split('-')[-1])
                    nthread_set.add(nthread)
                    measure_data =comp_data[solver_name][pfilename].attrs[measure_name]
                    nc =comp_data[solver_name][pfilename].attrs['nc']
                    n_iter =comp_data[solver_name][pfilename].attrs['iter']
                    if filename in results_filename.keys():
                        results_filename[filename].append([solver_name,nthread,measure_data,nc,n_iter])
                    else:
                        results_filename[filename] = [[solver_name,nthread,measure_data,nc,n_iter]]


            #print("n_filename", n_filename)

            results_filename_fails = {}
            for filename,solver in results_filename.items():
                for s in solver:
                    if np.isnan(s[2]):
                        if filename in results_filename.keys():
                            results_filename_fails[filename]=results_filename[filename]
                            print("\nremove failed instance for filename:",filename)
                            print(results_filename.pop(filename))




            nthread_list=list(nthread_set)
            nthread_list.sort()

            # collect results by thread

            measure_by_nthread=[]
            measure_penalized_by_nthread=[]
            measure_mean_by_nthread=[]
            measure_mean_penalized_by_nthread=[]
            index_non_failed=[]
            index_failed=[]

            for n in nthread_list:
                measure_by_nthread.append([])
                measure_mean_by_nthread.append(0.0)

            speedup_list =[]
            speedup_size_list =[]

            for n in nthread_list:
                speedup_list.append([])
                speedup_size_list.append([])

            for filename,solver in results_filename.items():
                for s in solver:
                        nthread= s[1]
                        thread_index=nthread_list.index(nthread)
                        measure_by_nthread[thread_index].append(s[2])
                        measure_mean_by_nthread[thread_index] += s[2]


            for n in nthread_list:
                thread_index    = nthread_list.index(n)
                measure_mean_by_nthread[thread_index] /= len(measure_by_nthread[thread_index])

            #raw_input()
            #print('measure_mean_by_nthread', measure_mean_penalized_by_nthread)

            speedup_list =[]
            speedup_size_list =[]
            iter_size_list=[]
            speedup_avg =[]

            for n in nthread_list:
                speedup_list.append([])
                speedup_size_list.append([])
                iter_size_list.append([])
                speedup_avg.append(0.0)
                thread_index=nthread_list.index(n)
                for i in range(len(measure_by_nthread[thread_index])):
                    speedup_list[thread_index].append(0.0)

            for n in nthread_list:
                thread_index_ref= nthread_list.index(1)
                thread_index    = nthread_list.index(n)
                for i  in range(len(measure_by_nthread[thread_index])):
                    speedup_list[thread_index][i] =  measure_by_nthread[thread_index_ref][i]/measure_by_nthread[thread_index][i]
                    speedup_avg[thread_index] += speedup_list[thread_index][i]
                speedup_avg[thread_index] /=len(measure_by_nthread[thread_index])

            print('speedup_avg', speedup_avg)
            #print('speedup_list', speedup_list)

            _cmp=0
            for filename,solver in results_filename.items():
                for s in solver:
                    thread_index=nthread_list.index(s[1])
                    #print _cmp, speedup_list[thread_index], speedup_list[thread_index][_cmp]
                    speedup_size_list[thread_index].append([s[3] , speedup_list[thread_index][_cmp]])
                    iter_size_list[thread_index].append([s[3] , s[4]])
                _cmp+=1

            count_failed=[]
            for n in nthread_list:
                count_failed.append(0)
            for filename,solver in results_filename_fails.items():
                for s in solver:
                    nthread=  s[1]
                    thread_index=nthread_list.index(nthread)
                    if np.isnan(s[2]):
                        count_failed[thread_index] +=1

            # ---------------------------- #
            # figure #
            # ---------------------------- #



            figure(figsize=(8,14))

            subplot('211')
            xlabel('problem size')
            ylabel('speedup for each thread')
            for n in nthread_list:
                index=nthread_list.index(n)
                list_size_speedup= speedup_size_list[index]
                list_size_speedup= sorted(list_size_speedup, key=lambda data: data[0])
                plot(np.array([t[0]  for t in list_size_speedup]), np.array([t[1]  for t in list_size_speedup]), label='n='+str(n))
            legend()
            subplot('212')
            xlabel('problem size')
            ylabel('iter for each thread')
            for n in nthread_list:
                index=nthread_list.index(n)
                iter_size_speedup= iter_size_list[index]
                iter_size_speedup= sorted(iter_size_speedup, key=lambda data: data[0])
                plot(np.array([t[0]  for t in iter_size_speedup]), np.array([t[1]  for t in iter_size_speedup]), label='n='+str(n))
            legend()

            figure(figsize=(8,14))

            subplot('311')
            boxplot(speedup_list,positions=nthread_list)
            ylabel('speedup distribution')
            legend()


            subplot('312')
            plot(nthread_list,speedup_avg)
            ylabel('avg. speed up')
            legend()

            subplot('313')
            bar(nthread_list,count_failed)
            ylabel('# fails')
            xlabel('number of threads')
            legend()

            figure(figsize=(16,14))

            for filename,solver in results_filename.items():
                data_tuples = []
                for s in solver:
                    data_tuples.append((s[1],s[2],s[3],s[4]))
                data_tuples=sorted(data_tuples, key=lambda data: data[0])
                try:
                    subplot('211')
                    #plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:], label =filename)
                    plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:])
                    ylabel('cpu time for each problems')
                    xlabel('number of threads')
                    legend()


                    subplot('212')
                    #plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:], label =filename)
                    plot(np.array([data[0] for data in data_tuples])[:],np.array([data_tuples[1][1]/data[1] for data in data_tuples])[:])
                    ylabel('speedup for each problems')
                    xlabel('number of threads')
                    legend()
                except:
                    pass




    display_bw=False
    if display or display_convergence or display_distrib or display_speedup:
        if not no_matplot:
            if (display_bw):
                figs = list(map(figure, get_fignums()))
                for fig in figs:
                    setFigLinesBW(fig)

            show()
