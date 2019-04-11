import sys
sys.path.append('../')
print(sys.path)
from SiconosSolver import *
from faf_tools import _read_fclib_format
import siconos.numerics as N
import numpy as np
import h5py


filename= 'Capsules-i97-269.hdf5'
fclib_problem, numerics_problem = _read_fclib_format(filename)



admm_br = SiconosSolver(name="ADMM-BR",
                        gnuplot_name="ADMM-BALANCING-RESIDUAL",
                        API=N.fc3d_admm,
                        TAG=N.SICONOS_FRICTION_3D_ADMM,
                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                        dparam_err=1,
                        maxiter=10,
                        precision=1e-04,
                        with_guess=False)

reactions, velocities, global_velocities = admm_br.guess(filename)

N.numerics_set_verbose(1)
result = admm_br(numerics_problem, reactions, velocities)

print('result', result)
print('reactions', reactions)




filename= 'Box_Stacks-i0035-31-0.hdf5'
fclib_problem, numerics_problem = _read_fclib_format(filename)




admm = SiconosSolver(name="ADMM",
                     gnuplot_name="ADMM",
                     API=N.gfc3d_ADMM,
                     TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                     iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                     dparam_err=1,
                     maxiter=10,
                     precision=1e-04,
                     with_guess=False,
                     global_solver=True)


reactions, velocities, global_velocities = admm.guess(filename)

N.numerics_set_verbose(1)
result = admm(numerics_problem, reactions, velocities, global_velocities)

print('result', result)
print('reactions', reactions)
