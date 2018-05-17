import h5py
from faf_tools import *
class Faf_preprocess():
     def __init__(self, filename,  problem_filenames):
          self._result_file=filename
          self._problem_filenames = problem_filenames

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





          

