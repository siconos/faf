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
        # options= nsgs_wr.SolverOptions()
        # N.solver_options_get_internal_solver(options,0).iparam[0]=nsgs_wr.SolverOptions().iparam[0]
        # N.solver_options_get_internal_solver(options,0).dparam[0]=nsgs_wr.SolverOptions().dparam[0]

        admm_wr = SiconosSolver(name="ADMM-WR",
                             gnuplot_name="ADMM-WR",
                             API=N.gfc3d_admm_wr,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM_WR,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)

        nsn_ac_wr = SiconosSolver(name="NSN-AC-WR",
                             gnuplot_name="NSN-AC-WR",
                             API=N.gfc3d_nonsmooth_Newton_AlartCurnier_wr,
                             TAG=N.SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                             global_solver=True)
       
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
        #aclm.SolverOptions().iparam[N.SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = N.SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
        #aclm.SolverOptions().iparam[N.SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = N.SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;
        
        
        admm_constant = SiconosSolver(name="ADMM-CST",
                                      gnuplot_name="ADMM-CST",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_constant.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_constant.SolverOptions().dparam[N.SICONOS_FRICTION_3D_ADMM_RHO]= 1e+02


        admm_norm_inf = SiconosSolver(name="ADMM-NORM-INF",
                                      gnuplot_name="ADMM-CST-NORM-INF",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_norm_inf.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_norm_inf.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO]= N.SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF
        
        admm_eig = SiconosSolver(name="ADMM-EIG",
                                      gnuplot_name="ADMM-CST-EIG",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_eig.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_eig.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO]= N.SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES

        admm_constant_no = SiconosSolver(name="ADMM-CST-NO",
                                         gnuplot_name="ADMM-CST-NO",
                                         API=N.gfc3d_ADMM,
                                         TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_constant_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_constant_no.SolverOptions().dparam[N.SICONOS_FRICTION_3D_ADMM_RHO]= 1.e+02
        admm_constant_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]=N.SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION
        
        admm_norm_inf_no = SiconosSolver(name="ADMM-NORM-INF-NO",
                                         gnuplot_name="ADMM-CST-NORM-INF-NO",
                                         API=N.gfc3d_ADMM,
                                         TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_norm_inf_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_norm_inf_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO]= N.SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF
        admm_norm_inf_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]=N.SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION
        
        admm_eig_no = SiconosSolver(name="ADMM-EIG-NO",
                                      gnuplot_name="ADMM-CST-EIG-NO",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_eig_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]= N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT
        admm_eig_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO]= N.SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES
        admm_eig_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]=N.SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION
        
        admm_br = SiconosSolver(name="ADMM-BR",
                                      gnuplot_name="ADMM-BR",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_br.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]=N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING

        admm_sbr = SiconosSolver(name="ADMM-SBR",
                                      gnuplot_name="ADMM-SBR",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_sbr.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]=N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING


        admm_br_fh = SiconosSolver(name="ADMM-BR-FH",
                                      gnuplot_name="ADMM-BR-FH",
                                      API=N.gfc3d_ADMM,
                                      TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                      global_solver=True)
        admm_br_fh.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]=N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
        admm_br_fh.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H]=N.SICONOS_FRICTION_3D_ADMM_FULL_H_YES

        admm_br_scaled = SiconosSolver(name="ADMM-BR-SCALED",
                                       gnuplot_name="ADMM-BR-SCALED",
                                       API=N.gfc3d_ADMM,
                                       TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                       iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                       dparam_err=1,
                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                       global_solver=True)
        admm_br_scaled.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]=N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
        admm_br_scaled.SolverOptions().iparam[N.SICONOS_FRICTION_3D_IPARAM_RESCALING]=N.SICONOS_FRICTION_3D_RESCALING_YES
        
        admm_br_no = SiconosSolver(name="ADMM-BR-NO",
                                       gnuplot_name="ADMM-BR-NO",
                                       API=N.gfc3d_ADMM,
                                       TAG=N.SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                       iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                       dparam_err=1,
                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess,
                                       global_solver=True)
        admm_br_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]=N.SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
        admm_br_no.SolverOptions().iparam[N.SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]=N.SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION
        
        ipm_solvers=True
        try :

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

        
        except:
            ipm_solvers=False
            
        nsgs_solvers = [nsgs, nsgs_wr, nsn_ac_wr, admm_wr]

        admm_solvers = [admm_constant, admm_norm_inf, admm_eig,
                        admm_constant_no, admm_norm_inf_no, admm_eig_no,
                        admm_br, admm_sbr, admm_br_scaled, admm_br_no, admm_br_fh]
        
        all_solvers = list(nsgs_solvers)
        all_solvers.extend(admm_solvers)
        
        
        all_solvers.extend([nsn_ac, vi_eg, vi_fp, aclm])

        if ipm_solvers:
            all_solvers.extend([ipm, ipm_nt])
        

        all_solvers = list(filter(lambda s : s is not None, all_solvers))

        return all_solvers
