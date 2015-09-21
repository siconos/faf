import os,sys

from pylmgc90.chipy import *
from numpy import *

checkDirectories()

### computation's parameters definition ### 

dt = 1e-3
nb_steps = 1000

# TT method
theta = 0.5

#DISPLAY
freq_display = 10
ref_radius   = 0.4

#SELECT
freq_detect = 1
PRPRx_SetCundallNeighbor(1e-1)
PRPRx_UseCpCundallDetection(200)
PRPRx_LowSizeArrayPolyr(10)

#GSNL
tol = 0.1666e-3
relax = 1.0
quad = 'QM/16'
gs_it1 = 51
gs_it2 = 501
nlgs_3D_DiagonalResolution()
freq_detect = 1  
itermax = 5000
tol = 1.e-5
relax = 1. #0.25 
solver='nlgs                          '
solver_output=3
SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,1,solver_output,10)
#############################
SetDimension(3)

utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

utilities_logMes('READ BODIES')
RBDY3_ReadBodies()
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

utilities_logMes('WRITE BODIES')
overall_WriteBodies()
RBDY3_WriteBodies()

utilities_logMes('WRITE BEHAVIOURS')
bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('WRITE DRIVEN DOF')
overall_WriteDrivenDof()
RBDY3_WriteDrivenDof()

### post3D ##
OpenDisplayFiles()
OpenPostproFiles()

### compute masses ###
utilities_logMes('COMPUTE MASS')
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

    #nlgs_3D_ExSolver('Exchange_Local_Global         ',quad, tol, relax, gs_it1, gs_it2)
    SiconosNumerics_ExSolver()


    PRPLx_StockRloc()
    PRPRx_StockRloc()
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

    ### writeout handling ###
    overall_CleanWriteOutFlags()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()
TimeEvolution_WriteLastVlocRloc()
PRPRx_WriteLastVlocRloc()

CloseDisplayFiles()
ClosePostproFiles()
