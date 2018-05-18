import numpy as np

measure = 'flop'

maxiter = 1000000
precision=1e-8
maxiterls = 10
utimeout = 10

global_problem = False

with_guess = True
with_mumps = 0

domain = np.arange(1, 100, .1)

clean = False
compute = True
display = False
display_convergence = False
display_distrib = False
display_distrib_var = False
display_speedup= False
no_matplot=False
logscale=False

gnuplot_output = False
gnuplot_with_color = True
gnuplot_separate_keys = False



output_dat=False
user_filenames = []
user_solvers = []
user_solvers_exact = []

keep_files = False
output_errors = False
output_velocities = False
output_reactions = False
measure_name = 'flpops'
ask_compute = True
ask_collect = True

ref_solver_name = 'NonsmoothGaussSeidel'
random_sample_proba = None
max_problems = None
cond_nc = None
file_filter=None
remove_file=None
list_contents=False
list_contents_solver = False
compute_hardness = False
compute_cond_rank = False
adhoc= False
thread_list = []
test_symmetry= False
compute_rho=False
estimate_optimal_timeout=False
numerics_has_openmp_solvers=False


forced=False
