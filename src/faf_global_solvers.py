from SiconosSolver import *
from faf_tools import *
 
#
# Some solvers
#
class faf_global_solvers():
    def __init__(self, maxiter, precision, maxiterls, with_guess, with_mumps, numerics_has_openmp_solvers):
        self._maxiter= maxiter
        self._precision=precision
        self._maxiterls=maxiterls
        self._with_guess=with_guess
        self._with_mumps=with_mumps
        self._numerics_has_openmp_solvers=numerics_has_openmp_solvers

    def create_solvers(self):

        ###### NSGS Family




        nsgs_wr = SiconosSolver(name="NSGS-WR",
                             gnuplot_name="NSGS-WR",
                             API=N.gfc3d_nsgs_wr,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_NSGS_WR,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        nsgs_wr.SolverOptions().internalSolvers.iparam[0]=nsgs_wr.SolverOptions().iparam[0]
        nsgs_wr.SolverOptions().internalSolvers.dparam[0]=nsgs_wr.SolverOptions().dparam[0]
        
        admm_wr = SiconosSolver(name="ADMM-WR",
                             gnuplot_name="ADMM-WR",
                             API=N.gfc3d_admm_wr,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM_WR,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        admm_wr.SolverOptions().internalSolvers.iparam[0]=admm_wr.SolverOptions().iparam[0]
        admm_wr.SolverOptions().internalSolvers.dparam[0]=admm_wr.SolverOptions().dparam[0]

        nsn_ac_wr = SiconosSolver(name="NSN-AC-WR",
                             gnuplot_name="NSN-AC-WR",
                             API=N.gfc3d_nonsmooth_Newton_AlartCurnier_wr,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        nsn_ac_wr.SolverOptions().internalSolvers.iparam[0]=nsn_ac_wr.SolverOptions().iparam[0]
        nsn_ac_wr.SolverOptions().internalSolvers.dparam[0]=nsn_ac_wr.SolverOptions().dparam[0]
        




        nsgs = SiconosSolver(name="NSGS-AC",
                             gnuplot_name="NSGS-AC",
                             API=N.gfc3d_nsgs,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        nsn_ac = SiconosSolver(name="NSN-AC",
                             gnuplot_name="NSN-AC",
                             API=N.gfc3d_nonsmooth_Newton_AlartCurnier,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_NSN_AC,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        
 
        vi_fp = SiconosSolver(name="VI-FP",
                             gnuplot_name="VI-FP",
                             API=N.gfc3d_VI_FixedPointProjection,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_VI_FPP,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        
        vi_eg = SiconosSolver(name="VI-EG",
                             gnuplot_name="VI-EG",
                             API=N.gfc3d_VI_ExtraGradient,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_VI_EG,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)

        aclm = SiconosSolver(name="ACLM",
                             gnuplot_name="ACLM",
                             API=N.gfc3d_ACLMFixedPoint,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_ACLMFP,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)

        
        admm = SiconosSolver(name="ADMM",
                             gnuplot_name="ADMM",
                             API=N.gfc3d_ADMM,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)

        
        ipm = SiconosSolver(name="IPM",
                             gnuplot_name="IPM",
                             API=N.gfc3d_IPM,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_IPM,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        ipm.SolverOptions().iparam[N.SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING]=0
        
        ipm_nt = SiconosSolver(name="IPM_NT",
                             gnuplot_name="IPM_NT",
                             API=N.gfc3d_IPM,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_IPM,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
        ipm_nt.SolverOptions().iparam[N.SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING]=1
        

        
        nsgs_solvers = [nsgs, nsgs_wr, nsn_ac_wr, admm_wr]

        all_solvers = list(nsgs_solvers)
        all_solvers.extend([nsn_ac, admm, vi_eg, vi_fp, aclm, ipm, ipm_nt])

        all_solvers = list(filter(lambda s : s is not None, all_solvers))

        return all_solvers
