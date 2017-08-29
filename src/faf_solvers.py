from SiconosSolver import *
from faf_tools import *
 
#
# Some solvers
#
class faf_solvers():
    def __init__(self, maxiter, precision, maxiterls, with_guess, with_mumps, numerics_has_openmp_solvers):
        self._maxiter= maxiter
        self._precision=precision
        self._maxiterls=maxiterls
        self._with_guess=with_guess
        self._with_mumps=with_mumps
        self._numerics_has_openmp_solvers=numerics_has_openmp_solvers

    def create_solvers(self):
        ###### NSN Family
        nsn_acSTD = SiconosSolver(name="NSN-AlartCurnier",
                                   gnuplot_name="NSN-AC",
                                   API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                   TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                   iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                   dparam_err=1,
                                   maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acSTD.SolverOptions().iparam[10] = 0;
        nsn_acSTD.SolverOptions().iparam[11] = 0;
        nsn_acSTD.SolverOptions().iparam[12] = self._maxiterls;
        nsn_acSTD.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acSTD.SolverOptions().iparam[3] = 10000000

        nsn_acSTD_nls = SiconosSolver(name="NSN-AlartCurnier-NLS",
                                       gnuplot_name="NSN-AC-NLS",
                                       API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                       TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                       iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                       dparam_err=1,
                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acSTD_nls.SolverOptions().iparam[10] = 0;
        nsn_acSTD_nls.SolverOptions().iparam[11] = -1;
        nsn_acSTD_nls.SolverOptions().iparam[12] = 0;
        nsn_acSTD_nls.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acSTD_nls.SolverOptions().iparam[3] = 10000000

        nsn_acJeanMoreau = SiconosSolver(name="NSN-JeanMoreau",
                                          gnuplot_name="NSN-JM",
                                          API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                          TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                          iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                          dparam_err=1,
                                          maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acJeanMoreau.SolverOptions().iparam[10] = 1;
        nsn_acJeanMoreau.SolverOptions().iparam[11] = 0;
        nsn_acJeanMoreau.SolverOptions().iparam[12] = self._maxiterls;
        nsn_acJeanMoreau.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acJeanMoreau.SolverOptions().iparam[3] = 10000000

        nsn_acJeanMoreau_nls = SiconosSolver(name="NSN-JeanMoreau-NLS",
                                              gnuplot_name="NSN-JM-NLS",
                                              API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                              TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                              iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                              dparam_err=1,
                                              maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acJeanMoreau_nls.SolverOptions().iparam[10] = 1;
        nsn_acJeanMoreau_nls.SolverOptions().iparam[11] = -1;
        nsn_acJeanMoreau_nls.SolverOptions().iparam[12] = 0;
        nsn_acJeanMoreau_nls.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acJeanMoreau_nls.SolverOptions().iparam[3] = 10000000

        nsn_acSTDGenerated = SiconosSolver(name="NSN-AlartCurnier-Generated",
                                            gnuplot_name="NSN-AC-Generated",
                                            API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                            TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                            iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                            dparam_err=1,
                                            maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acSTDGenerated.SolverOptions().iparam[10] = 2;
        nsn_acSTDGenerated.SolverOptions().iparam[11] = 0;
        nsn_acSTDGenerated.SolverOptions().iparam[12] = self._maxiterls;
        nsn_acSTDGenerated.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acSTDGenerated.SolverOptions().iparam[3] = 10000000

        nsn_acSTDGenerated_nls = SiconosSolver(name="NSN-AlartCurnier-Generated-NLS",
                                                gnuplot_name="NSN-AC-Generated-NLS",
                                                API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                dparam_err=1,
                                                maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acSTDGenerated_nls.SolverOptions().iparam[10] = 2;
        nsn_acSTDGenerated_nls.SolverOptions().iparam[11] = -1;
        nsn_acSTDGenerated_nls.SolverOptions().iparam[12] = 0;
        nsn_acSTDGenerated_nls.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acSTDGenerated_nls.SolverOptions().iparam[3] = 10000000

        nsn_acJeanMoreauGenerated = SiconosSolver(name="NSN-JeanMoreau-Generated",
                                                   gnuplot_name="NSN-JM-Generated",
                                                   API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                   TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                   iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                   dparam_err=1,
                                                   maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acJeanMoreauGenerated.SolverOptions().iparam[10] = 3;
        nsn_acJeanMoreauGenerated.SolverOptions().iparam[11] = 0;
        nsn_acJeanMoreauGenerated.SolverOptions().iparam[12] = self._maxiterls;
        nsn_acJeanMoreauGenerated.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acJeanMoreauGenerated.SolverOptions().iparam[3] = 10000000


        nsn_acJeanMoreauGenerated_lusol = None
        if self._with_mumps:
            # reference
            nsn_acJeanMoreauGenerated_lusol = SiconosSolver(name="NSN-JeanMoreau-Generated-lusol",
                                                             gnuplot_name="NSN-JM-Generated-LUSOL",
                                                             API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                             TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                             dparam_err=1,
                                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)


            nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[10] = 3;
            nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[11] = 0;
            nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[12] = self._maxiterls;
            nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[13] = 0;
            nsn_acJeanMoreauGenerated_lusol.SolverOptions().iparam[3] = 10000000


        nsn_acJeanMoreauGenerated_nls = SiconosSolver(name="NSN-JeanMoreau-Generated-NLS",
                                                       gnuplot_name="NSN-JM-Generated-NLS",
                                                       API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                       TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                       iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                       dparam_err=1,
                                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[10] = 3;
        nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[11] = -1;
        nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[12] = self._maxiterls;
        nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acJeanMoreauGenerated_nls.SolverOptions().iparam[3] = 10000000

        nsn_acJeanMoreauGenerated_nls_lusol=None
        if self._with_mumps:
            nsn_acJeanMoreauGenerated_nls_lusol = SiconosSolver(name="NSN-JeanMoreau-Generated-NLS-lusol",
                                                                 gnuplot_name="NSN-JM-Generated-NLS-LUSOL",
                                                                 API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                                                 TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                                                 iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                                 dparam_err=1,
                                                                 maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

            nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[10] = 3;
            nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[11] = -1;
            nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[12] = self._maxiterls;
            nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[13] = 0;
            nsn_acJeanMoreauGenerated_nls_lusol.SolverOptions().iparam[3] = 10000000

        nsn_fb_gp = SiconosSolver(name="NSN-FischerBurmeister-GP",
                                   gnuplot_name="NSN-FB-GP",
                                   API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                   TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                   iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                   dparam_err=1,
                                   maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_fb_gp.SolverOptions().iparam[3] = 1000000
        nsn_fb_gp.SolverOptions().iparam[11] = 0
        nsn_fb_gp.SolverOptions().iparam[12] = self._maxiterls
        nsn_fb_gp.SolverOptions().iparam[13] = self._with_mumps

        nsn_fb_gp_lusol = None
        if self._with_mumps:
            nsn_fb_gp_lusol = SiconosSolver(name="NSN-FischerBurmeister-GP-lusol",
                                             gnuplot_name="NSN-FB-GP-LUSOL",
                                             API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                             TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

            nsn_fb_gp_lusol.SolverOptions().iparam[3] = 1000000
            nsn_fb_gp_lusol.SolverOptions().iparam[11] = 0
            nsn_fb_gp_lusol.SolverOptions().iparam[12] = self._maxiterls
            nsn_fb_gp_lusol.SolverOptions().iparam[13] = 0

        nsn_fb_fblsa = SiconosSolver(name="NSN-FischerBurmeister-FBLSA",
                                      gnuplot_name="NSN-FB-FBLSA",
                                      API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                      TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_fb_fblsa.SolverOptions().iparam[3] = 1000000
        nsn_fb_fblsa.SolverOptions().iparam[11] = 1
        nsn_fb_fblsa.SolverOptions().iparam[12] = self._maxiterls
        nsn_fb_fblsa.SolverOptions().iparam[13] = self._with_mumps

        nsn_fb_nls = SiconosSolver(name="NSN-FischerBurmeister-NLS",
                                    gnuplot_name="NSN-FB-NLS",
                                    API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                    TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                    iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                    dparam_err=1,
                                    maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_fb_nls.SolverOptions().iparam[3] = 1000000
        nsn_fb_nls.SolverOptions().iparam[11] = -1
        nsn_fb_nls.SolverOptions().iparam[12] = 0
        nsn_fb_nls.SolverOptions().iparam[13] = self._with_mumps

        nsn_fb_nls_lusol = None
        if self._with_mumps:
            nsn_fb_nls_lusol = SiconosSolver(name="NSN-FischerBurmeister-NLS-lusol",
                                              gnuplot_name="NSN-FB-NLS-LUSOL",
                                              API=N.fc3d_nonsmooth_Newton_FischerBurmeister,
                                              TAG=N.SICONOS_FRICTION_3D_NSN_FB,
                                              iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                              dparam_err=1,
                                              maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

            nsn_fb_nls_lusol.SolverOptions().iparam[3] = 1000000
            nsn_fb_nls_lusol.SolverOptions().iparam[11] = -1
            nsn_fb_nls_lusol.SolverOptions().iparam[12] = 0
            nsn_fb_nls_lusol.SolverOptions().iparam[13] = 0


        nsn_nm_gp = SiconosSolver(name="NSN-NaturalMap-GP",
                                   gnuplot_name="NSN-NM-GP",
                                   API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                   TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                   iparam_iter=1,
                                   dparam_err=1,
                                   maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_nm_gp.SolverOptions().iparam[3] = 1000000
        nsn_nm_gp.SolverOptions().iparam[11] = 0
        nsn_nm_gp.SolverOptions().iparam[12] = self._maxiterls
        nsn_nm_gp.SolverOptions().iparam[13] = self._with_mumps

        nsn_nm_gp_lusol = None
        if self._with_mumps:
            nsn_nm_gp_lusol = SiconosSolver(name="NSN-NaturalMap-GP-lusol",
                                             gnuplot_name="NSN-NM-GP-LUSOL",
                                             API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                             TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

            nsn_nm_gp_lusol.SolverOptions().iparam[3] = 1000000
            nsn_nm_gp_lusol.SolverOptions().iparam[11] = 0
            nsn_nm_gp_lusol.SolverOptions().iparam[12] = self._maxiterls
            nsn_nm_gp_lusol.SolverOptions().iparam[13] = 0

        nsn_nm_fblsa = SiconosSolver(name="NSN-NaturalMap-FBLSA",
                                      gnuplot_name="NSN-NM-FBLSA",
                                      API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                      TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                      iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                      dparam_err=1,
                                      maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_nm_fblsa.SolverOptions().iparam[3] = 1000000
        nsn_nm_fblsa.SolverOptions().iparam[11] = 1
        nsn_nm_fblsa.SolverOptions().iparam[12] = self._maxiterls
        nsn_nm_fblsa.SolverOptions().iparam[13] = self._with_mumps

        nsn_nm_nls = SiconosSolver(name="NSN-NaturalMap-NLS",
                                    gnuplot_name="NSN-NM-NLS",
                                    API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                    TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                    iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                    dparam_err=1,
                                    maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_nm_nls.SolverOptions().iparam[3] = 1000000
        nsn_nm_nls.SolverOptions().iparam[11] = -1
        nsn_nm_nls.SolverOptions().iparam[12] = 0
        nsn_nm_nls.SolverOptions().iparam[13] = self._with_mumps

        nsn_acSTD_nls_hybrid = SiconosSolver(name="NSN-AlartCurnier-NLS-HYBRID",
                                             gnuplot_name="NSN-AC-NLS-HYBRID",
                                             API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                             TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        nsn_acSTD_nls_hybrid.SolverOptions().iparam[10] = 0;
        nsn_acSTD_nls_hybrid.SolverOptions().iparam[12] = 0;
        nsn_acSTD_nls_hybrid.SolverOptions().iparam[13] = self._with_mumps;
        nsn_acSTD_nls_hybrid.SolverOptions().iparam[3] =  10000000
        nsn_acSTD_nls_hybrid.SolverOptions().iparam[N.SICONOS_FRICTION_3D_NSN_LINESEARCH] = N.SICONOS_FRICTION_3D_NSN_LINESEARCH_NO;
        nsn_acSTD_nls_hybrid.SolverOptions().iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY]=N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN;

        nsn_nm_nls_lusol = None
        if self._with_mumps:
            nsn_nm_nls_lusol = SiconosSolver(name="NSN-NaturalMap-NLS-lusol",
                                              gnuplot_name="NSN-NM-NLS-LUSOL",
                                              API=N.fc3d_nonsmooth_Newton_NaturalMap,
                                              TAG=N.SICONOS_FRICTION_3D_NSN_NM,
                                              iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                              dparam_err=1,
                                              maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

            nsn_nm_nls_lusol.SolverOptions().iparam[3] = 1000000
            nsn_nm_nls_lusol.SolverOptions().iparam[11] = -1
            nsn_nm_nls_lusol.SolverOptions().iparam[12] = 0
            nsn_nm_nls_lusol.SolverOptions().iparam[13] = 0


        hnsn_ac = SiconosHybridSolver(name = "HLocalAlartCurnier",
                                       API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                       TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                       iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                       dparam_err=1,
                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        hnsn_ac.SolverOptions().iparam[3] = 10000000

        ###### NSGS Family


        nsgs = SiconosSolver(name="NSGS-AC",
                             API=N.fc3d_nsgs,
                             TAG=N.SICONOS_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
        nsgs.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD


        nsgs_ac_gp = SiconosSolver(name="NSGS-AC-GP",
                             API=N.fc3d_nsgs,
                             TAG=N.SICONOS_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP
        nsgs_ac_gp.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD




        nsgs_jm = SiconosSolver(name="NSGS-JM",
                             API=N.fc3d_nsgs,
                             TAG=N.SICONOS_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_jm.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
        nsgs_jm.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD


        nsgs_jm_gp = SiconosSolver(name="NSGS-JM-GP",
                             API=N.fc3d_nsgs,
                             TAG=N.SICONOS_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_jm_gp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP
        nsgs_jm_gp.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD


        snsgs = SiconosSolver(name="NSGS-AC-GP-Shuffled",
                              gnuplot_name="NSGS-AC-GP Shuffled",
                              API=N.fc3d_nsgs,
                              TAG=N.SICONOS_FRICTION_3D_NSGS,
                              iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                              dparam_err=1,
                              maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        snsgs.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP
        snsgs.SolverOptions().iparam[5] = 1
        #print snsgs.SolverOptions().iparam[6]

        nsgs_sfull = SiconosSolver(name="NSGS-AC-GP-Shuffled-full",
                                   gnuplot_name="NSGS-AC-GP Fully shuffled",
                                   API=N.fc3d_nsgs,
                                   TAG=N.SICONOS_FRICTION_3D_NSGS,
                                   iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                   dparam_err=1,
                                   maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_sfull.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP
        nsgs_sfull.SolverOptions().iparam[5] = 2
        #nsgs_sfull.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
        #N.printSolverOptions(nsgs_sfull.SolverOptions())

        nsgs_pli = SiconosSolver(name="NSGS-PLI-100",
                                 gnuplot_name="NSGS-FP-VI-UPK iter=100",
                                 API=N.fc3d_nsgs,
                                 TAG=N.SICONOS_FRICTION_3D_NSGS,
                                 iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                 dparam_err=1,
                                 maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_pli.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
        nsgs_pli.SolverOptions().internalSolvers.iparam[0] = 100

        nsgs_pli_10 = SiconosSolver(name="NSGS-PLI-10",
                                    gnuplot_name="NSGS-FP-VI-UPK iter=10",
                                    API=N.fc3d_nsgs,
                                    TAG=N.SICONOS_FRICTION_3D_NSGS,
                                    iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                    dparam_err=1,
                                    maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_pli_10.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
        nsgs_pli_10.SolverOptions().internalSolvers.iparam[0] = 10


        nsgs_p = SiconosSolver(name="NSGS-P",
                               gnuplot_name="NSGS-FP-DS-One",
                               API=N.fc3d_nsgs,
                               TAG=N.SICONOS_FRICTION_3D_NSGS,
                               iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                               dparam_err=1,
                               maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_p.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone

        nsgs_pd = SiconosSolver(name="NSGS-PD",
                                gnuplot_name="NSGS-FP-DS-One  D",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                dparam_err=1,
                                maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_pd.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization

        nsgs_pr = SiconosSolver(name="NSGS-PR",
                                gnuplot_name="NSGS-FP-DS-One  R",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                dparam_err=1,
                                maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_pr.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization




        nsgs_ac_gp_hybrid_pli_nsn = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-1",
                                     API=N.fc3d_nsgs,
                                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                                     iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                     dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 1
        nsgs_ac_gp_hybrid_pli_nsn.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP

        nsgs_ac_gp_hybrid_pli_nsn_10 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10",
                                     API=N.fc3d_nsgs,
                                     TAG=N.SICONOS_FRICTION_3D_NSGS,
                                     iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                     dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_10.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_10.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_10.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 10
        nsgs_ac_gp_hybrid_pli_nsn_10.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP

        nsgs_ac_gp_hybrid_pli_nsn_100 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100",
                                         API=N.fc3d_nsgs,
                                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP

        nsgs_ac_gp_hybrid_pli_nsn_10_1 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10-1",
                                         API=N.fc3d_nsgs,
                                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_10_1.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_10_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_10_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 10
        nsgs_ac_gp_hybrid_pli_nsn_10_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1
        nsgs_ac_gp_hybrid_pli_nsn_10_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP

        nsgs_ac_gp_hybrid_pli_nsn_100_1 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-1",
                                           API=N.fc3d_nsgs,
                                           TAG=N.SICONOS_FRICTION_3D_NSGS,
                                           iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                           dparam_err=1,
                                           maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP

        nsgs_ac_gp_hybrid_pli_nsn_100 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100",
                                         API=N.fc3d_nsgs,
                                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_pli_nsn_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP



        nsgs_ac_gp_hybrid_pli_nsn_100_1 = SiconosSolver(name="NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-1",
                                           API=N.fc3d_nsgs,
                                           TAG=N.SICONOS_FRICTION_3D_NSGS,
                                           iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                           dparam_err=1,
                                           maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1
        nsgs_ac_gp_hybrid_pli_nsn_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY ] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP


        nsgs_ac_gp_hybrid_nsn_nsn_pli_100 = SiconosSolver(name="NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100",
                                         API=N.fc3d_nsgs,
                                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] = N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_NSN_PLI_LOOP 



        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1 = SiconosSolver(name="NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100-1",
                                                            API=N.fc3d_nsgs,
                                                            TAG=N.SICONOS_FRICTION_3D_NSGS,
                                                            iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                                            dparam_err=1,
                                                            maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_FORMULATION]=N.SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1
        nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1.SolverOptions().internalSolvers.iparam[N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] =  N.SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_NSN_PLI_LOOP 







        nsgs_openmp_solvers=[]

 
        if (self._numerics_has_openmp_solvers):
            error_evaluation_frequency=1

            nsgs_openmp = SiconosSolver(name="NSGS-AC-OPENMP-FOR-"+str(error_evaluation_frequency)+"-"+str(0),
                                        API=N.fc3d_nsgs,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            nsgs_openmp.SolverOptions().iparam[1] = N.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
            nsgs_openmp.SolverOptions().iparam[14]=error_evaluation_frequency
            nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
            nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0
            nsgs_openmp_solvers.append(nsgs_openmp)


            for n in thread_list:
                nsgs_openmp = SiconosSolver(name="NSGS-AC-OPENMP-FOR-"+str(error_evaluation_frequency)+"-"+str(n),
                                            API=N.fc3d_nsgs_openmp,
                                            TAG=N.SICONOS_FRICTION_3D_NSGS_OPENMP,
                                            iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                            dparam_err=1,
                                            maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                nsgs_openmp.SolverOptions().iparam[10]=n
                nsgs_openmp.SolverOptions().iparam[11]=0
                nsgs_openmp.SolverOptions().iparam[14]=error_evaluation_frequency
                nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
                nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0

                nsgs_openmp_solvers.append(nsgs_openmp)

            for n in thread_list:
                nsgs_openmp = SiconosSolver(name="NSGS-AC-OPENMP-REDBLACK-"+str(n),
                                            API=N.fc3d_nsgs_openmp,
                                            TAG=N.SICONOS_FRICTION_3D_NSGS_OPENMP,
                                            iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                            dparam_err=1,
                                            maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                nsgs_openmp.SolverOptions().iparam[10]=n
                nsgs_openmp.SolverOptions().iparam[11]=1
                nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
                nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0

                nsgs_openmp_solvers.append(nsgs_openmp)

            for n in thread_list:
                nsgs_openmp = SiconosSolver(name="NSGS-AC-OPENMP-DDM-NAIVE-"+str(n),
                                            API=N.fc3d_nsgs_openmp,
                                            TAG=N.SICONOS_FRICTION_3D_NSGS_OPENMP,
                                            iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                            dparam_err=1,
                                            maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                nsgs_openmp.SolverOptions().iparam[10]=n
                nsgs_openmp.SolverOptions().iparam[11]=2
                nsgs_openmp.SolverOptions().iparam[12]=5
                nsgs_openmp.SolverOptions().iparam[13]=10
                nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
                nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0

                nsgs_openmp_solvers.append(nsgs_openmp)


            nsgs_openmp = SiconosSolver(name="NSGS-ERROR-COMPARISON",
                                        API=N.fc3d_nsgs_openmp,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS_OPENMP,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            nsgs_openmp.SolverOptions().iparam[10]=1
            nsgs_openmp.SolverOptions().iparam[11]=10
            nsgs_openmp.SolverOptions().iparam[14]=error_evaluation_frequency
            nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
            nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0
            nsgs_openmp_solvers.append(nsgs_openmp)

            nsgs_openmp = SiconosSolver(name="NSGS-AC-OPENMP-REDBLACK-"+str(0),
                                        API=N.fc3d_nsgs_openmp,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS_OPENMP,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            nsgs_openmp.SolverOptions().iparam[10]=1
            nsgs_openmp.SolverOptions().iparam[11]=1
            nsgs_openmp.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_NSN
            nsgs_openmp.SolverOptions().internalSolvers.iparam[10]=0
            nsgs_openmp_solvers.append(nsgs_openmp)

        local_tol_values = [1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16]
        #local_tol_values = [1e-2,1e-6,1e-10,1e-16]
        nsgs_series=[]
        for local_tol in local_tol_values:
            str1 = "{0:1.0e}".format(local_tol).replace("1e","10\^{")+"}"
            nsgs_solver = SiconosSolver(name="NSGS-AC-GP-"+str(local_tol),
                                        gnuplot_name="NSGS-AC-GP \$tol\_{local}="+str1+"\$",
                                        API=N.fc3d_nsgs,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            nsgs_solver.SolverOptions().internalSolvers.dparam[0] = local_tol
            nsgs_series.append(nsgs_solver)

        for local_tol in local_tol_values:
            str1 = "{0:1.0e}".format(local_tol).replace("1e","10\^{")+"}"
            nsgs_solver = SiconosSolver(name="NSGS-PLI-"+str(local_tol),
                                        gnuplot_name="NSGS-FP-VI-UPK \$tol\_{local}="+str1+"\$",
                                        API=N.fc3d_nsgs,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            nsgs_pli.SolverOptions().internalSolvers.solverId = N.SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
            nsgs_solver.SolverOptions().internalSolvers.dparam[0] = local_tol
            nsgs_solver.SolverOptions().internalSolvers.iparam[0] = 100
            nsgs_series.append(nsgs_solver)


        snsgs_series=[]
        for i in range(10):
            snsgs_solver = SiconosSolver(name="NSGS-AC-Shuffled-"+str(i),
                                         gnuplot_name="NSGS-AC-GP Shuffled "+str(i),
                                         API=N.fc3d_nsgs,
                                         TAG=N.SICONOS_FRICTION_3D_NSGS,
                                         iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            snsgs_solver.SolverOptions().iparam[5] = 1
            snsgs_solver.SolverOptions().iparam[6] = (1237*i)*(1237*i)
            #print snsgs_solver.SolverOptions().iparam[6]
            snsgs_series.append(snsgs_solver)




        # only dense
        nsgsv = SiconosSolver(name="NSGS-Velocity",
                              API=N.fc3d_nsgs_velocity,
                              TAG=N.SICONOS_FRICTION_3D_NSGSV,
                              iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                              dparam_err=1,
                              maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)


        omega=1.5
        psor = SiconosSolver(name="PSOR-AC",
                             gnuplot_name="PSOR-AC",
                             API=N.fc3d_nsgs,
                             TAG=N.SICONOS_FRICTION_3D_NSGS,
                             iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        psor.SolverOptions().iparam[4] = 1
        psor.SolverOptions().dparam[8] = omega

        omega_values = [0.5, 0.8, 1.0, 1.1, 1.3, 1.5, 1.8]
        psor_series=[]
        for omega in omega_values:
            psor_solver = SiconosSolver(name="PSOR-AC-"+str(omega),
                                        gnuplot_name="PSOR-AC \$\\\omega="+str(omega)+"\$",
                                        API=N.fc3d_nsgs,
                                        TAG=N.SICONOS_FRICTION_3D_NSGS,
                                        iparam_iter=N.SICONOS_IPARAM_ITER_DONE,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            psor_solver.SolverOptions().iparam[4] = 1
            psor_solver.SolverOptions().dparam[8] = omega
            psor_series.append(psor_solver)

        TrescaFixedPoint = SiconosSolver(name="TrescaFixedPoint-NSGS-PLI",
                                         API=N.fc3d_TrescaFixedPoint,
                                         TAG=N.SICONOS_FRICTION_3D_TFP,
                                         iparam_iter=7,
                                         dparam_err=1,
                                         maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        ACLMFixedPoint = SiconosSolver(name="ACLMFixedPoint-SOCLCP-NSGS-PLI",
                                       API=N.fc3d_ACLMFixedPoint,
                                       TAG=N.SICONOS_FRICTION_3D_ACLMFP,
                                       iparam_iter=7,
                                       dparam_err=1,
                                       maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        SOCLCP = SiconosSolver(name="SOCLCP-NSGS-PLI",
                               API=N.fc3d_SOCLCP,
                               TAG=N.SICONOS_FRICTION_3D_SOCLCP,
                               iparam_iter=7,
                               dparam_err=1,
                               maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        DeSaxceFixedPoint = SiconosSolver(name="FixedPoint-DeSaxce",
                                          API=N.fc3d_DeSaxceFixedPoint,
                                          TAG=N.SICONOS_FRICTION_3D_DSFP,
                                          iparam_iter=7,
                                          dparam_err=1,
                                          maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        ExtraGrad = SiconosSolver(name="ExtraGradient",
                                  API=N.fc3d_ExtraGradient,
                                  TAG=N.SICONOS_FRICTION_3D_EG,
                                  iparam_iter=7,
                                  dparam_err=1,
                                  maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        FixedPointProjection = SiconosSolver(name="FixedPoint-Projection",
                                             API=N.fc3d_fixedPointProjection,
                                             TAG=N.SICONOS_FRICTION_3D_FPP,
                                             iparam_iter=7,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        VIExtraGrad = SiconosSolver(name="ExtraGradient-VI",
                                    API=N.fc3d_VI_ExtraGradient,
                                    TAG=N.SICONOS_FRICTION_3D_VI_EG,
                                    iparam_iter=7,
                                    dparam_err=1,
                                    maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        # VIExtraGrad.SolverOptions().dparam[4]=0.6
        # VIExtraGrad.SolverOptions().dparam[5]=1/0.7
        # VIExtraGrad.SolverOptions().dparam[6]=0.9
        # VIExtraGrad.SolverOptions().dparam[7]=0.3

        iparam1_values = [0,1]

        iparam2_values = [0,1]

        iparam3_values = [0]

        VIExtraGrad_series=[]
        for i1 in iparam1_values:
            for i2 in iparam2_values:
                for i3 in iparam3_values:
                    if i1 == 0 :
                        g_name="EG-VI-UPK"
                    elif i1 == 1 :
                        g_name="EG-VI-UPTS"
                    elif i1 == 2:
                        g_name="EG-VI-UPHS"

                    if i2 == 0:
                        g_name = g_name + " False"
                    elif i2 == 1:
                        g_name = g_name + " True"

                    if i3 == 0:
                        g_name = g_name
                    elif i3 == 1:
                        g_name = g_name + " min"

                    VIExtraGrad_solver= SiconosSolver(name="ExtraGrad-VI-"+str(i1)+str(i2)+str(i3),
                                                                 gnuplot_name=g_name,
                                                                 API=N.fc3d_VI_ExtraGradient,
                                                                 TAG=N.SICONOS_FRICTION_3D_VI_EG,
                                                                 iparam_iter=7,
                                                                 dparam_err=1,
                                                                 maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                    VIExtraGrad_solver.SolverOptions().iparam[1] = i1
                    VIExtraGrad_solver.SolverOptions().iparam[2] = i2
                    VIExtraGrad_solver.SolverOptions().iparam[3] = i3
                    VIExtraGrad_series.append(VIExtraGrad_solver)


        VIFixedPointProjection = SiconosSolver(name="FixedPoint-VI",
                                               API=N.fc3d_VI_FixedPointProjection,
                                               TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                               iparam_iter=7,
                                               dparam_err=1,
                                               maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        iparam1_values = [0,1,2]
        iparam1_values = [0,1]

        iparam2_values = [0]
        iparam2_values = [0,1 ]


        iparam3_values = [0,1]
        iparam3_values = [0]
        VIFixedPointProjection_series=[]
        for i1 in iparam1_values:
            for i2 in iparam2_values:
                for i3 in iparam3_values:
                    if i1 == 0 :
                        g_name="FP-VI-UPK"
                    elif i1 == 1 :
                        g_name="FP-VI-UPTS"
                    elif i1 == 2:
                        g_name="FP-VI-UPHS"

                    if i2 == 0:
                        g_name = g_name + " False"
                    elif i2 == 1:
                        g_name = g_name + " True"
                    elif i2 == 2:
                        g_name = g_name + " semi False"

                    if i3 == 0:
                        g_name = g_name
                    elif i3 == 1:
                        g_name = g_name + " min"

                    VIFixedPointProjection_solver= SiconosSolver(name="FixedPoint-VI-"+str(i1)+str(i2)+str(i3),
                                                                 gnuplot_name=g_name,
                                                                 API=N.fc3d_VI_FixedPointProjection,
                                                                 TAG=N.SICONOS_FRICTION_3D_VI_FPP,
                                                                 iparam_iter=7,
                                                                 dparam_err=1,
                                                                 maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                    VIFixedPointProjection_solver.SolverOptions().iparam[1] = i1
                    VIFixedPointProjection_solver.SolverOptions().iparam[2] = i2
                    VIFixedPointProjection_solver.SolverOptions().iparam[3] = i3
                    VIFixedPointProjection_series.append(VIFixedPointProjection_solver)

        Prox = SiconosSolver(name="PROX-NSN-AC",
                             gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        Prox.SolverOptions().internalSolvers.iparam[3] = 100


        Proxfixed = SiconosSolver(name="PROX-NSN-AC-fixed",
                             gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu=1, \\\sigma=5.0\$ fixed",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        Proxfixed.SolverOptions().dparam[3] = -1e4
        Proxfixed.SolverOptions().internalSolvers.iparam[3] = 100


        regul_value=1e03
        str1 = "{0:1.0e}".format(regul_value).replace("1e","10\^{")+"}"
        str2 = "{0:1.0e}".format(regul_value)
        Regul_variable = SiconosSolver(name="PROX-NSN-AC-regulVar-"+str2,
                             gnuplot_name="PPA-NSN-AC-GP   \$ \\\mu=1, \\\sigma=5.0\$ regulVar "+str1,
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        Regul_variable.SolverOptions().dparam[3] = regul_value
        Regul_variable.SolverOptions().iparam[9] = 1
        Regul_variable.SolverOptions().internalSolvers.iparam[3] = 100


        regul_values= [ 1e4, 1e6, 1e8, 1e10]
        regul_series =[]
        for rr in regul_values:
            str1 = "{0:1.0e}".format(rr).replace("1e","10\^{")+"}"
            str2 = "{0:1.0e}".format(rr)
            regul_solver  = SiconosSolver(name="PROX-NSN-AC-regul-"+str2,
                                          gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu="+str1+"\$",
                                          API=N.fc3d_proximal,
                                          TAG=N.SICONOS_FRICTION_3D_PROX,
                                          iparam_iter=7,
                                          dparam_err=1,
                                          maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
            regul_solver.SolverOptions().dparam[3] = - rr
            regul_solver.SolverOptions().iparam[9] = 1
            regul_solver.SolverOptions().internalSolvers.iparam[3] = 100
            regul_series.append(regul_solver)

        Prox_nls = SiconosSolver(name="PROX-NSN-AC-NLS",
                             gnuplot_name="PPA-NSN-AC  \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        Prox_nls.SolverOptions().internalSolvers.iparam[11] = -1

        ProxFB = SiconosSolver(name="PROX-NSN-FB-GP",
                             gnuplot_name="PPA-NSN-FB-GP  \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        ProxFB.SolverOptions().internalSolvers.iparam[3] = 100
        ProxFB.SolverOptions().dparam[4]=5.0 # sigma
        ProxFB.SolverOptions().dparam[5]=1.0 # nu
        nsn_fb_gp_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
        nsn_fb_gp_inprox.iparam[3] = 1000000
        nsn_fb_gp_inprox.iparam[11] = 0
        nsn_fb_gp_inprox.iparam[12] = 6 #self._maxiterls
        nsn_fb_gp_inprox.iparam[13] = self._with_mumps

        ProxFB.SolverOptions().internalSolvers = nsn_fb_gp_inprox
        ProxFB.SolverOptions().internalSolvers.iparam[3] = 100



        ProxFB_fblsa = SiconosSolver(name="PROX-NSN-FB-FBLSA",
                             gnuplot_name="PPA-NSN-FB-FBLSA  \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsn_fb_fblsa_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
        nsn_fb_fblsa_inprox.iparam[3] = 1000000
        nsn_fb_fblsa_inprox.iparam[11] = 1
        nsn_fb_fblsa_inprox.iparam[12] = 6 #self._maxiterls
        nsn_fb_fblsa_inprox.iparam[13] = self._with_mumps
        ProxFB_fblsa.SolverOptions().internalSolvers = nsn_fb_fblsa_inprox
        ProxFB_fblsa.SolverOptions().internalSolvers.iparam[3] = 1000000

        ProxFB_nls = SiconosSolver(name="PROX-NSN-FB-NLS",
                             gnuplot_name="PPA-NSN-FB-NLS  \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        nsn_fb_nls_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)
        nsn_fb_nls_inprox.iparam[3] = 1000000
        nsn_fb_nls_inprox.iparam[11] = -1
        nsn_fb_nls_inprox.iparam[13] = self._with_mumps
        ProxFB_nls.SolverOptions().internalSolvers = nsn_fb_nls_inprox
        ProxFB_nls.SolverOptions().internalSolvers.iparam[3] = 1000000

        ProxNSGS = SiconosSolver(name="PROX-NSGS-NSN-AC",
                             gnuplot_name="PROX-NSGS-NSN-AC \$ \\\mu=1, \\\sigma=5.0\$",
                             API=N.fc3d_proximal,
                             TAG=N.SICONOS_FRICTION_3D_PROX,
                             iparam_iter=7,
                             dparam_err=1,
                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        local_nsgs_inprox = N.SolverOptions(N.SICONOS_FRICTION_3D_NSGS)
        ProxNSGS.SolverOptions().internalSolvers = local_nsgs_inprox
        ProxNSGS.SolverOptions().internalSolvers.iparam[0] = 10000







        sigmavalues= [0.5, 1.0, 4.0, 5.0, 50, 100.0, 1000.0 ]
        muvalues= [0.5, 1.0, 2.0]
        prox_series =[]
        for mu in muvalues:
            for sigma in sigmavalues:
                prox_solver  = SiconosSolver(name="PROX-NSN-AC-nu"+str(mu)+"-sigma"+str(sigma),
                                             gnuplot_name="PPA-NSN-AC-GP  \$ \\\mu="+str(mu)+", \\\sigma="+str(sigma)+"\$",
                                             API=N.fc3d_proximal,
                                             TAG=N.SICONOS_FRICTION_3D_PROX,
                                             iparam_iter=7,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
                prox_solver.SolverOptions().internalSolvers.iparam[3] = 1000000
                prox_solver.SolverOptions().dparam[4]=sigma # sigma
                prox_solver.SolverOptions().dparam[5]=mu # nu
                prox_series.append(prox_solver)

        nsn_ac_wrapped = SiconosSolver(name="NSN-AlartCurnier-Wrapped",
                                        API=N.fc3d_nonsmooth_Newton_AlartCurnier,
                                        TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                        iparam_iter=1,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        def fc3d_nsn_ac_r(problem, reactions, velocities, _SO):
            SO = N.SolverOptions(N.SICONOS_FRICTION_3D_VI_FPP)
            SO.iparam[3] = 1000
            N.fc3d_VI_FixedPointProjection(problem, reactions, velocities, SO)
            #    print '->',SO.dparam[3]
            nsn_ac_wrapped.SolverOptions().dparam[3] = SO.dparam[3]
            return nsn_ac_wrapped(problem, reactions, velocities)

        # flop measure only on nsn_ac
        nsn_acr = SiconosWrappedSolver(name="NSN-AlartCurnier-R",
                                        API=fc3d_nsn_ac_r,
                                        TAG=N.SICONOS_FRICTION_3D_NSN_AC,
                                        iparam_iter=1,
                                        dparam_err=1,
                                        maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        #
        quartic = SiconosSolver(name="NSGS-Quartic",
                                API=N.fc3d_nsgs,
                                TAG=N.SICONOS_FRICTION_3D_NSGS,
                                iparam_iter=7,
                                dparam_err=1, maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        quartic3x3 = N.SolverOptions(N.SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU)

        quartic.SolverOptions().internalSolvers = quartic3x3

        # 1 contact
        #AlartCurnierNewton = SiconosSolver(name="AlartCurnierNewton",
        #                                   API=fc3d_AlartCurnierNewton,
        #                                   iparam_iter=1,
        #                                   iparam_iter=1)

        # rho estimation needed
        Prox._SO.dparam[3] = 1000

        HyperplaneProjection = SiconosSolver(name="HyperplaneProjection",
                                             API=N.fc3d_HyperplaneProjection,
                                             TAG=N.SICONOS_FRICTION_3D_HP,
                                             iparam_iter=7,
                                             dparam_err=1,
                                             maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)

        bogusPureNewton = BogusSolver(name="BogusPureNewton", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        bogusPureNewton.SolverOptions().iparam[4]=0


        bogusPureEnumerative = BogusSolver(name="BogusPureEnumerative", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        bogusPureEnumerative.SolverOptions().iparam[4]=1


        bogusHybrid = BogusSolver(name="BogusHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        bogusHybrid.SolverOptions().iparam[4]=2


        bogusRevHybrid = BogusSolver(name="BogusRevHybrid", API=wrap_bogus_solve, TAG=N.SICONOS_FRICTION_3D_NSN_FB, iparam_iter=1, dparam_err=1, maxiter=self._maxiter, precision=self._precision, with_guess=self._with_guess)
        bogusRevHybrid.SolverOptions().iparam[4]=3



        # check for correct flop
        #def pipo(*args,**kwargs):
        #    a=1.
        #    b=3.
        #    c = a+b
        #    return 0
        #
        #Pipo = SiconosSolver(name="pipo",
        #                     API=pipo,
        #                     TAG=SICONOS_FRICTION_3D_HP,
        #                     iparam_iter=1,
        #                     dparam_err=1)
        # http://code.google.com/p/mpi4py/
        #from mpi4py import MPI
        #rank = MPI.COMM_WORLD.rank
        #fc3d_sparseGlobalAlartCurnierInit(nsn_ac.SolverOptions())

        #all_solvers = [nsgs, snsgs, TrescaFixedPoint, nsn_ac, Prox, DeSaxceFixedPoint,
        #               FixedPointProjection, VIFixedPointProjection, ExtraGrad, VIExtraGrad]
        #all_solvers = [nsgs, snsgs, quartic, TrescaFixedPoint, ACLMFixedPoint, DeSaxceFixedPoint, VIFixedPointProjection, VIFixedPointProjection1, VIFixedPointProjection2, VIFixedPointProjection3, VIExtraGrad, SOCLCP, Prox, Prox2, Prox3, Prox4, Prox5, nsn_acSTD, nsn_acSTDGenerated,  nsn_acr, nsn_acJeanMoreau, nsn_acJeanMoreauGenerated, nsn_fb_gp, nsn_fb_fblsa]

        nsgs_solvers = [nsgs, nsgs_ac_gp,  nsgs_jm, nsgs_jm_gp, nsgs_sfull, snsgs,
                        nsgs_pli, nsgs_pli_10,
                        nsgs_p, nsgs_pd, nsgs_pr , quartic,
                        nsgs_ac_gp_hybrid_pli_nsn, nsgs_ac_gp_hybrid_pli_nsn_10, nsgs_ac_gp_hybrid_pli_nsn_100, nsgs_ac_gp_hybrid_pli_nsn_10_1, nsgs_ac_gp_hybrid_pli_nsn_100_1,
                        nsgs_ac_gp_hybrid_nsn_nsn_pli_100,  nsgs_ac_gp_hybrid_nsn_nsn_pli_100_1]
        # remove very nasty solver
        #nsgs_solvers.remove(nsgs_p)
        nsgs_solvers.remove(nsgs_pd)
        nsgs_solvers.remove(nsgs_pr)
        nsgs_solvers.remove(quartic)

        nsn_solvers =  [nsn_acSTD, nsn_acSTD_nls, nsn_acSTDGenerated, nsn_acSTDGenerated_nls,  nsn_acr, nsn_acJeanMoreau, nsn_acJeanMoreau_nls, nsn_acJeanMoreauGenerated, nsn_acJeanMoreauGenerated_lusol,
                        nsn_acJeanMoreauGenerated_nls, nsn_acJeanMoreauGenerated_nls_lusol,
                        nsn_fb_gp, nsn_fb_gp_lusol, nsn_fb_nls, nsn_fb_nls_lusol,
                        nsn_nm_gp, nsn_nm_gp_lusol, nsn_nm_nls, nsn_nm_nls_lusol,
                        nsn_acSTD_nls_hybrid]

        all_solvers = list(nsgs_solvers)
        all_solvers.extend(nsn_solvers)
        all_solvers.extend( [ psor,
                              TrescaFixedPoint, DeSaxceFixedPoint,
                              VIFixedPointProjection, VIExtraGrad,
                              SOCLCP,
                              Prox,  Prox_nls, ProxFB,  ProxFB_nls, ProxNSGS, Proxfixed, Regul_variable, regul_series[0],
                              ACLMFixedPoint])

        all_solver_unstable = [ProxFB_fblsa]
        all_solvers.extend(all_solver_unstable)

        all_solvers = list(filter(lambda s : s is not None, all_solvers))

        ###
        # specific studies of solvers.

        #all_solvers.extend(VIFixedPointProjection_series)
        #all_solvers.extend(VIExtraGrad_series)
        #all_solvers.extend(psor_series)
        #all_solvers.extend(prox_series)
        #all_solvers.extend(regul_series)
        #all_solvers.extend(nsgs_series)

        all_solvers.extend(nsgs_openmp_solvers)
        return all_solvers
