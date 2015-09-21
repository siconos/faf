import os,sys

from pylmgc90.chipy import *
from numpy import *

checkDirectories()

### computation's parameters definition ### 

nb_steps = 50
dt = 1e-4
theta = 0.5

freq_display = 10
ref_radius   = 0.4

#POLYR_SkipHEBuild()
#POLYR_SkipAutomaticReorientation()

PRPRx_UseCpCundallDetection(150)
#PRPRx_UseCpF2fExplicitDetection(1.e-2)
PRPRx_ShrinkPolyrFaces(1e-2)
PRPRx_LowSizeArrayPolyr(40)

freq_detect = 1
# norm='Exchange_Local_Global         '
# tol = 1e-4
# relax = 1.0
# quad = 'Maxm '
# gs_it1 = 200
# gs_it2 = 10
# nlgs_3D_SetWithQuickScramble()
# nlgs_3D_DiagonalResolution()
freq_detect = 1  
itermax = 5000
tol = 1.e-5
relax = 1. #0.25 
#       123456789012345678901234567890
solver='nlgs                          '
solver_output=3
SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,1,solver_output,1)



##############################

SetDimension(3,0)

utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### model reading ###
utilities_logMes('READ BODIES')
RBDY3_ReadBodies()

utilities_logMes('LOAD TACTORS')
PLANx_LoadTactors()
POLYR_LoadTactors()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

RBDY3_LoadBehaviours()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY3_ReadIniDof()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
PRPLx_ReadIniVlocRloc()
PRPRx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY3_ReadDrivenDof()

overall_WriteBodies()
RBDY3_WriteBodies()

bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

overall_WriteDrivenDof()
RBDY3_WriteDrivenDof()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()

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
    PRPLx_SelectProxTactors()
    PRPRx_SelectProxTactors()
    
    PRPLx_RecupRloc()
    PRPRx_RecupRloc()

    SiconosNumerics_ExSolver()

    PRPLx_StockRloc()
    PRPRx_StockRloc()

    utilities_logMes('COMPUTE DOF')
    RBDY3_ComputeDof()
    #
    utilities_logMes('UPDATE DOF')
    TimeEvolution_UpdateStep()
    RBDY3_UpdateDof()

    ### post3D ###
    WriteDisplayFiles(freq_display,ref_radius)
    WritePostproFiles()

    ### writeout handling ###
    overall_CleanWriteOutFlags()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()

TimeEvolution_WriteLastVlocRloc()
PRPLx_WriteLastVlocRloc()
PRPRx_WriteLastVlocRloc()

CloseDisplayFiles()
ClosePostproFiles()
