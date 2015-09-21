import os,sys

from pylmgc90.chipy import *
from numpy import *

checkDirectories()

### computation's parameters definition ### 

dt = 2e-2
nb_steps = 50000
theta = 0.5

# source point
first_object=2
radius=6.
coorx=12.
coory=12.
coorz=12.

# periodicity
xperiode = 25.0
yperiode = 25.0

freq_display = 500
ref_radius=1.

RBDY3_NewRotationScheme()

PRPRx_UseCpCundallDetection(300)
PRPRx_LowSizeArrayPolyr(70)

freq_detect = 1
tol = 0.1666e-4
relax = 1.0
quad = 'Maxm '
gs_it1 = 400
gs_it2 = 10
type='Exchange_Local_Global         '
#       123456789012345678901234567890
freq_detect = 1  
itermax = 5000
tol = 1.e-5
relax = 1. #0.25 
solver='nlgs                          '
solver_output=3
SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,1,solver_output,10)
#nlgs_3D_SetWithQuickScramble()

##########################################

SetDimension(3)

utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### model reading ###
utilities_logMes('READ BODIES')
RBDY3_ReadBodies()
PLANx_LoadTactors()
POLYR_LoadTactors()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
RBDY3_LoadBehaviours()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()

PRPLx_ReadIniVlocRloc()
PRPRx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY3_ReadDrivenDof()

overall_WriteBodies()
RBDY3_WriteBodies()
#MAILx_WriteBodies()

bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

overall_WriteDrivenDof()
RBDY3_WriteDrivenDof()

# source point
# skip,radius,dx,dy,dz
RBDY3_SetSourcePointWithIni(first_object, radius, coorx, coory, coorz)

### set periodic conditions ###
RBDY3_SetXPeriodicCondition(xperiode)
PRPRx_SetXPeriodicCondition(xperiode)
post3D_SetXPeriodicCondition(xperiode)

RBDY3_SetYPeriodicCondition(yperiode)
PRPRx_SetYPeriodicCondition(yperiode)
post3D_SetYPeriodicCondition(yperiode)


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

    #nlgs_3D_ExSolver(type, quad, tol, relax, gs_it1, gs_it2)
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

    ### postpro ###
    WritePostproFiles()

    ### writeout handling ###
    overall_CleanWriteOutFlags()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()
#
TimeEvolution_WriteLastVlocRloc()
PRPRx_WriteLastVlocRloc()

CloseDisplayFiles()
ClosePostproFiles()
