import h5py
from faf_tools import *
class Faf_postprocess():
     def __init__(self, filename, solvers, problem_filenames):
          self._result_file=filename
          self._solvers=solvers
          self._problem_filenames = problem_filenames
          
     def compute_rho(self, file_filter, remove_file, random_sample_proba, max_problems,
                     cond_nc,  measure_name, domain):
          filename=None
          solver_r = dict()
          measure = dict()
          min_measure = dict()
          for fileproblem in self._problem_filenames:
               min_measure[fileproblem] = np.inf
          with h5py.File(self._result_file, 'r') as comp_file:

               data = comp_file['data']
               comp_data = data['comp']
               comp_precision=comp_file['data']['comp'].attrs.get('precision')
               comp_utimeout=comp_file['data']['comp'].attrs.get('timeout')
               # 1 n_problems
               n_problems = 0

               for solver in self._solvers:
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

               for solver in self._solvers:
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
                              except :
                                   print('raise exception')
                                   measure[solver_name][ip] = np.nan
                              ip += 1
               # 3 solver_r
               #            for solver_name in comp_data:
               for solver in self._solvers:
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
               for solver in self._solvers:
                    solver_name=solver.name()
                    if solver_name in comp_data :
                         assert min(solver_r[solver_name]) >= 1
                         rhos[solver_name] = np.empty(len(domain))
                         for itau in range(0, len(domain)):
                              rhos[solver_name][itau] = float(len(np.where( solver_r[solver_name] <= domain[itau] )[0])) / float(n_problems)


          return rhos,solver_r,measure, min_measure, filenames, filename


     def estimate_optimal_timeout(self):
        max_rhos=dict()
        with h5py.File(self._result_file, 'r') as comp_file:
            data = comp_file['data']
            comp_data = data['comp']
            for solver in self._solvers:
                #print(solver)
                solver_name=solver.name()
                #print(rhos[solver_name])
                #print(rhos[solver_name])
                if solver_name in comp_data :
                    max_rhos[solver_name]= np.max(rhos[solver_name])
            #print(max_rhos)
            level_of_success=0.5
            nb_succeeded_solver= len(np.argwhere(np.array([max_rhos[solver_name] for solver_name in max_rhos.keys() ]) > level_of_success))
            #print(nb_succeeded_solver, "solvers has rho_max over", level_of_success)
            print("{0:2.0f} % solvers suceeded to reach rho max equal {1} ".format(nb_succeeded_solver/len(self._solvers)*100,level_of_success))
            level_of_success=0.9
            nb_succeeded_solver= len(np.argwhere(np.array([max_rhos[solver_name] for solver_name in max_rhos.keys() ]) > level_of_success))
            #print(nb_succeeded_solver, "solvers has rho_max over", level_of_success)
            print("{0:2.0f} % solvers suceeded to reach rho max equal {1} ".format(nb_succeeded_solver/len(self._solvers)*100,level_of_success))
            level_of_success=0.99
            nb_succeeded_solver= len(np.argwhere(np.array([max_rhos[solver_name] for solver_name in max_rhos.keys() ]) > level_of_success))
            #print(nb_succeeded_solver, "solvers has rho_max over", level_of_success)
            print("{0:2.0f} % solvers suceeded to reach rho max equal {1} ".format(nb_succeeded_solver/len(self._solvers)*100,level_of_success))


     def compute_cond_rank(self):
        print("Tasks will be run for", self._problem_filenames)
        for problem_filename in self._problem_filenames:
            print("compute for", self._problem_filename,"....")
            with h5py.File(problem_filename, 'r+') as fclib_file:
                no_rank_info=True
                if (not forced):
                    if (fclib_file['fclib_local']['W'].attrs.get('rank') == None) :
                        print("Rank info already not  in", problem_filename)
                    else:
                        print("Rank info already in", problem_filename)
                        print("fclib_file['fclib_local']['W'].attrs.get('rank')", fclib_file['fclib_local']['W'].attrs.get('rank'))
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

     def test_symmetry(self):
          print("Tasks will be run for", self._problem_filenames)
          symmetry_test_list=[]
          dd_test_list=[]
          for problem_filename in self._problem_filenames:
               print("compute for", problem_filename,"....")
               with h5py.File(problem_filename, 'r+') as fclib_file:
                    no_rank_info=True
                
                    try:
                         is_symmetric, symmetry_test, is_dd, dd_test = test_symmetry_W(problem_filename)
                         symmetry_test_list.append(symmetry_test)
                         dd_test_list.append(dd_test)
                         print("is_symmetric", is_symmetric)
                         print("is_dd",is_dd)
                         with h5py.File(problem_filename, 'r+') as fclib_file:
                              fclib_file['fclib_local']['W'].attrs.create('symmetry', is_symmetric)
                    except Exception as e :
                         print("-->", e)
          print("avg symmetry_test",np.array(symmetry_test_list).mean() )
          print("avg dd_test",np.array(dd_test_list).mean() )


     def compute_hardness(self):
          nc = []
          nds = []
          cond_nc = []
          max_measure = dict()
          max_measure_by_contact = dict()
          min_measure_by_contact = dict()

          for fileproblem in self._problem_filenames:
               max_measure[fileproblem] = - np.inf
               max_measure_by_contact[fileproblem] = - np.inf
               min_measure_by_contact[fileproblem] =  np.inf



          for problem_filename in self._problem_filenames:

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

                         ip = 0

                         for filename in filenames:
                              if filename not in min_measure:
                                   min_measure[filename] = np.inf
                              #try:
                              pfilename = os.path.splitext(filename)[0]
                              #print(' pfilename',pfilename)
                              #print(solver_name,pfilename,comp_data[solver_name][pfilename].attrs['info'])
                              if comp_data[solver_name][pfilename].attrs['info'] == 0:
                                   n_contact=numberOfDegreeofFreedomContacts(filename)/3
                                   #print(' filename n_contact',filename,n_contact)
                                   measure[solver_name][ip] =  comp_data[solver_name][pfilename].attrs[measure_name]
                                   min_measure[filename] = min(min_measure[filename], measure[solver_name][ip])
                                   max_measure[filename] = max(max_measure[filename], measure[solver_name][ip])
                                   min_measure_by_contact[filename] = min(min_measure_by_contact[filename], measure[solver_name][ip]/n_contact)
                                   max_measure_by_contact[filename] = max(max_measure_by_contact[filename], measure[solver_name][ip]/n_contact)
                              else:
                                   measure[solver_name][ip] = np.inf
                              #except:
                              #    measure[solver_name][ip] = np.nan
                              ip += 1

               print("min_measure", min_measure)
               #print("max_measure", max_measure)
               unsolved =0
               min_measure_list=[]
               min_measure_by_contact_list=[]
               for key in min_measure.keys():
                    if min_measure[key] ==np.inf:
                         unsolved += 1
                    else:
                         min_measure_list.append(min_measure[key])
                         min_measure_by_contact_list.append(min_measure_by_contact[key])
               print('number of unsolved problems', unsolved)
               #min_measure_array=np.array([min_measure[key] for key in min_measure.keys()])
               #min_measure_by_contact_array=np.array([min_measure_by_contact[key] for key in min_measure_by_contact.keys()])
            
               avg_min_measure = np.array(min_measure_list).mean()
               std_min_measure = np.array(min_measure_list).std()                
               avg_min_measure_by_contact = np.array(min_measure_by_contact_list).mean()
               std_min_measure_by_contact = np.array(min_measure_by_contact_list).std()                
               print(         "Average min resolution measure (avg fastest solver measure) = {0:12.8e}".format(avg_min_measure))
               print(         "Std min resolution measure (std fastest solver measure) = {0:12.8e}".format(std_min_measure))
               print(         "Average min resolution measure by contact = {0:12.8e}".format(avg_min_measure_by_contact))
               print(         "Std min resolution measure by contact = {0:12.8e}".format(std_min_measure_by_contact))


               max_measure_list=[]
               max_measure_by_contact_list=[]
               for key in min_measure.keys():
                    if max_measure[key] == -np.inf:
                         unsolved += 1
                    else:
                         max_measure_list.append(max_measure[key])
                         max_measure_by_contact_list.append(max_measure_by_contact[key])
            
               avg_max_measure = np.array(max_measure_list).mean()
               std_max_measure = np.array(max_measure_list).std()    
               print(         "Average max resolution measure (avg slowest suceeded solver measure) = {0:12.8e}".format(avg_max_measure))
               print(         "Std max resolution measure (std fastest solver measure) = {0:12.8e}".format(std_max_measure))
               print(         "Average max resolution measure by contact = {0:12.8e}".format(avg_max_measure/nc_avg))



          

