import os,sys

from pylmgc90.chipy import *
from numpy import *

checkDirectories()

### computation's parameters definition ###
dt = 0.002
nb_steps = 500

theta = 0.5

freq_display = 10
ref_radius   = 0.01

RBDY3_NewRotationScheme()

freq_detect = 1


# tol = 0.1666e-3
# relax = 1.0
# quad = 'Maxm '
# gs_it1 = 11
# gs_it2 = 101
# nlgs_3D_SetWithQuickScramble()


itermax = 5000
tol = 1.e-5
relax = 1. #0.25
#       123456789012345678901234567890
solver='nlgs                          '
solver_output=3
freq_output=10
SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,1,solver_output,10)



SetDimension(3,0)

utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### model reading ###
utilities_logMes('READ BODIES')
RBDY3_ReadBodies()
PLANx_LoadTactors()
SPHER_LoadTactors()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
RBDY3_LoadBehaviours()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY3_ReadIniDof()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
SPPLx_ReadIniVlocRloc()
SPSPx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY3_ReadDrivenDof()

overall_WriteBodies()
RBDY3_WriteBodies()

bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

overall_WriteDrivenDof()
RBDY3_WriteDrivenDof()

### post3D ##
OpenDisplayFiles()
OpenPostproFiles()

### compute masses ###
RBDY3_ComputeMass()

for k in range(nb_steps):
    #
    utilities_logMes('INCREMENT STEP')
    TimeEvolution_IncrementStep()
    RBDY3_IncrementStep()
    #
    utilities_logMes('DISPLAY TIMES')
    TimeEvolution_DisplayStep()
    #
    utilities_logMes('COMPUTE Fext')
    RBDY3_ComputeFext()
    #
    utilities_logMes('COMPUTE Fint')
    RBDY3_ComputeBulk()
    #
    utilities_logMes('COMPUTE Free Vlocy')
    RBDY3_ComputeFreeVelocity()
    #
    utilities_logMes('SELECT PROX TACTORS')
    overall_SelectProxTactors(freq_detect)
    SPPLx_SelectProxTactors()
    SPSPx_SelectProxTactors()

    SPPLx_RecupRloc()
    SPSPx_RecupRloc()

    #nlgs_3D_ExSolver('Exchange_Local_Global         ',quad, tol, relax, gs_it1, gs_it2)
    SiconosNumerics_ExSolver()


    SPPLx_StockRloc()
    SPSPx_StockRloc()
    #
    utilities_logMes('COMPUTE DOF')
    RBDY3_ComputeDof()
    #
    utilities_logMes('UPDATE DOF')
    TimeEvolution_UpdateStep()
    RBDY3_UpdateDof()

    ### post3D ###
    WriteDisplayFiles(freq_display,ref_radius)
    WritePostproFiles()

    ### postpro ###
    postpro_3D_PostproDuringComputation()

    ### writeout handling ###
    overall_CleanWriteOutFlags()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()

TimeEvolution_WriteLastVlocRloc()
SPSPx_WriteLastVlocRloc()
SPPLx_WriteLastVlocRloc()

CloseDisplayFiles()
ClosePostproFiles()
