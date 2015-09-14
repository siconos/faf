import sys

from pylmgc90.chipy import *

checkDirectories()

#utilities_DisableLogMes()

Initialize()
SetDimension(3,0)

####
# info gestion du temps
dt = 1.e-4
theta = 0.505
nb_steps = 200

# bavardage de certaines fonctions
echo = 0

# info generation fichier visu
freq_display = 1
ref_radius = 0.1

# info contact
freq_detect = 1

#       123456789012345678901234567890
itermax = 5000
tol = 1.e-8  
relax = 1. #0.25 
#       123456789012345678901234567890
solver='nsgs                          '
SiconosNumerics_SetParameters(solver,tol,itermax,relax,0,0)

###

#mecaMAILx_SparseStorage()

### lecture du modele ###
utilities_logMes('READ BODIES')
MAILx_ReadBodies()
RBDY3_ReadBodies()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

utilities_logMes('READ MODELS')
models_ReadModels()

utilities_logMes('INIT MODELS')
# on dimensionne et on initie la construction des mapping
models_InitModels()
ExternalModels_InitModels()

utilities_logMes('LOADS')
# on charge les choses et on construit les mapping
mecaMAILx_LoadModels()
mecaMAILx_LoadBehaviours()
RBDY3_LoadBehaviours()

utilities_logMes('PUSH')
mecaMAILx_PushProperties()
#
utilities_logMes('STORE')
# on finalise la construction des mapping
models_StoreProperties()

utilities_logMes('CHECK')
ExternalModels_CheckProperties()

#
POLYR_LoadTactors()
CSxxx_LoadTactors()
ASpxx_LoadTactors()
#
utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
mecaMAILx_ReadIniDof()
RBDY3_ReadIniDof()

utilities_logMes('READ INI GPV')
TimeEvolution_ReadIniGPV()
mecaMAILx_ReadIniGPV()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
CSASp_ReadIniVlocRloc()
CSPRx_ReadIniVlocRloc()

#
utilities_logMes('READ DRIVEN DOF')
mecaMAILx_ReadDrivenDof()
RBDY3_ReadDrivenDof()

### ecriture paranoiaque du modele ###
utilities_logMes('WRITE BODIES')
overall_WriteBodies()
MAILx_WriteBodies()
RBDY3_WriteBodies()

utilities_logMes('WRITE MODELS')
models_WriteModels()

utilities_logMes('WRITE BEHAVIOURS')
bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('WRITE DRIVEN DOF')
overall_WriteDrivenDof()
mecaMAILx_WriteDrivenDof()
RBDY3_WriteDrivenDof()

utilities_logMes('WRITE out DOF')
TimeEvolution_WriteOutDof(1)
mecaMAILx_WriteOutDof()
RBDY3_WriteOutDof(-1,100)

utilities_logMes('WRITE LAST DOF')
TimeEvolution_WriteLastDof()
mecaMAILx_WriteLastDof()
RBDY3_WriteLastDof()

utilities_logMes('WRITE out Vloc_Rloc')
TimeEvolution_WriteOutVlocRloc(1)
CSASp_WriteOutVlocRloc()
CSPRx_WriteOutVlocRloc()

utilities_logMes('WRITE LAST Vloc_Rloc')
TimeEvolution_WriteLastVlocRloc()
CSASp_WriteLastVlocRloc()
CSPRx_WriteLastVlocRloc()

### definition des parametres du calcul ### 
utilities_logMes('INIT TIME STEPPING')

TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### post3D ##
OpenDisplayFiles()
### parameters setting ###

utilities_logMes('COMPUTE MASS')
mecaMAILx_ComputeMass()
RBDY3_ComputeMass()

utilities_logMes('COMPUTE STIFFNESS')
mecaMAILx_ComputeBulk()
mecaMAILx_AssembKT()

# precondensation
mecaMAILx_SetPreconAllBodies()
CSxxx_PushPreconNodes()
ASpxx_PushPreconNodes()
mecaMAILx_ComputePreconW()

for k in xrange(1,nb_steps+1,1):
   #
   utilities_logMes('increment : '+str(k))
   #
   utilities_logMes('INCREMENT STEP')
   TimeEvolution_IncrementStep()
   mecaMAILx_IncrementStep()
   RBDY3_IncrementStep()
   
   utilities_logMes('DISPLAY TIMES')
   TimeEvolution_DisplayStep()

   utilities_logMes('COMPUTE Fext')
   mecaMAILx_ComputeFext()
   RBDY3_ComputeFext()

   utilities_logMes('COMPUTE Fint')
   mecaMAILx_ComputeBulk()
   RBDY3_ComputeBulk()

   utilities_logMes('ASSEMBLAGE')
   mecaMAILx_AssembRHS()

   utilities_logMes('COMPUTE Free Vlocy')
   mecaMAILx_ComputeFreeVelocity()
   RBDY3_ComputeFreeVelocity()
   #
   utilities_logMes('SELECT PROX TACTORS')
   overall_SelectProxTactors(freq_detect)
   CSASp_SelectProxTactors()
   CSPRx_SelectProxTactors()
   #
   CSASp_RecupRloc()
   CSPRx_RecupRloc()

   utilities_logMes('RESOLUTION' )
    
   #utilities_logMes(type,quad,tol, relax, gs_it1, gs_it2 )
   SiconosNumerics_ExSolver()
    
   CSASp_StockRloc()
   CSPRx_StockRloc()
   #
   utilities_logMes('COMPUTE DOF, FIELDS, etc.')
   mecaMAILx_ComputeDof()
   RBDY3_ComputeDof()
   #
   utilities_logMes('UPDATE DOF, FIELDS')
   TimeEvolution_UpdateStep()
   mecaMAILx_UpdateDof()
   mecaMAILx_UpdateBulk()
   RBDY3_UpdateDof()
   #
   utilities_logMes('WRITE out DOF')
   TimeEvolution_WriteOutDof(100)
   mecaMAILx_WriteOutDof()
   RBDY3_WriteOutDof(-1,100)
   #
   utilities_logMes('WRITE out Rloc')
   TimeEvolution_WriteOutVlocRloc(100)
   CSASp_WriteOutVlocRloc()
   CSPRx_WriteOutVlocRloc()

   utilities_logMes('WRITE LAST DOF')
   TimeEvolution_WriteLastDof()
   mecaMAILx_WriteLastDof()
   RBDY3_WriteLastDof()
   #
   utilities_logMes('WRITE Vloc_Rloc')
   TimeEvolution_WriteLastVlocRloc()
   CSASp_WriteLastVlocRloc()
   CSPRx_WriteLastVlocRloc()

   ### post3D ###
   WriteDisplayFiles(freq_display,ref_radius)

   ### postpro ###
   #postpro_3D_PostproDuringComputation()

   ### gestion des writeout ###
   overall_CleanWriteOutFlags()

### postpro ###
CloseDisplayFiles()

Finalize()
