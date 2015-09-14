import os,sys

from pylmgc90.chipy import *
import numpy as np

checkDirectories()

#utilities_DisableLogMes()

Initialize()
SetDimension(3,0)

### computation's parameters definition ### 
utilities_logMes('INIT TIME STEPPING')
nb_iter = 10
dt = 1e-3  
theta = 0.5  

POLYR_SkipAutomaticReorientation()

PRPRx_UseCpF2fExplicitDetection(1e-3)
PRPRx_ShrinkPolyrFaces(0.05)
PRPRx_LowSizeArrayPolyr(100)

nlgs_3D_DiagonalResolution()

freq_detect = 1  
itermax = 5000
tol = 1.e-8  
relax = 1. #0.25 
#       123456789012345678901234567890
solver='nlgs                          '

SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,0,0)

freq_visu = 1
ref_radius = 0.1

TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### model reading ###
utilities_logMes('READ BODIES')
RBDY3_ReadBodies()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
RBDY3_LoadBehaviours()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY3_ReadIniDof()

POLYR_LoadTactors()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
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


OpenDisplayFiles()


### compute masses ###
utilities_logMes('COMPUTE MASS')
RBDY3_ComputeMass()

for k in range(nb_iter):
    #
    utilities_logMes('itere : '+str(k))
    #
    utilities_logMes('INCREMENT STEP')
    TimeEvolution_IncrementStep()
    RBDY3_IncrementStep()
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
    PRPRx_SelectProxTactors()
    

    PRPRx_RecupRloc()

    SiconosNumerics_ExSolver()

    PRPRx_StockRloc()

    #
    #
    utilities_logMes('COMPUTE DOF')
    RBDY3_ComputeDof()
    #

    ### post3D ###
    WriteDisplayFiles(freq_visu,ref_radius)
    # visualisation des joints
    if k % freq_visu == 0 :
      PRPRx_VisavisVTKDrawAll()



    ### postpro ###

    utilities_logMes('UPDATE DOF')
    TimeEvolution_UpdateStep()
    RBDY3_UpdateDof()

    TimeEvolution_WriteOutDof(freq_visu)
    RBDY3_WriteOutDof(-1,10000)
    TimeEvolution_WriteOutVlocRloc(freq_visu)
    PRPRx_WriteOutVlocRloc()
    ### writeout handling ###
    overall_CleanWriteOutFlags()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()


# postpro_3D_ClosePostproFiles()
CloseDisplayFiles()

