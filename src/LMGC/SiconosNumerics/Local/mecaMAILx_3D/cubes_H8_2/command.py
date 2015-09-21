
from pylmgc90.chipy import *

checkDirectories()

#utilities_DisableLogMes()

####
# info gestion du temps
dt = 1.e-2
theta = 0.505
nb_steps = 200

# bavardage de certaines fonctions
echo = 0

# info generation fichier visu
freq_display = 1
ref_radius   = 0.1

freq_write = 1

# info contact
freq_detect = 1

#CSxxx_SetQuadrature(1)

# #       123456789012345678901234567890
# type = 'Stored_Delassus_Loops         '
# norm = 'QM/16'
# tol = 1e-5
# relax = 1.0
# gs_it1 = 50
# gs_it2 = 10

freq_detect = 1  
itermax = 5000
tol = 1.e-4
relax = 1. #0.25 
#       123456789012345678901234567890
solver='nlgs                          '
solver_output=3
SiconosNumerics_SetParameters(solver,tol,10,itermax,relax,1,solver_output,1)




###
SetDimension(3,0)

### definition des parametres du calcul ### 
utilities_logMes('INIT TIME STEPPING')

TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### lecture du modele ###
utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

utilities_logMes('READ MODELS')
models_ReadModels()

utilities_logMes('READ BODIES')
MAILx_ReadBodies()

# on dimensionne et on initie la construction des mapping
utilities_logMes('INIT MODELS')
models_InitModels()
ExternalModels_InitModels()

# on charge les choses et on construit les mapping
utilities_logMes('LOADS')
mecaMAILx_LoadModels()
mecaMAILx_LoadBehaviours()

utilities_logMes('PUSH')
mecaMAILx_PushProperties()

# on finalise la construction des mapping
utilities_logMes('STORE')
models_StoreProperties()

utilities_logMes('CHECK')
ExternalModels_CheckProperties()

#
utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
mecaMAILx_ReadIniDof()

utilities_logMes('READ INI GPV')
TimeEvolution_ReadIniGPV()
mecaMAILx_ReadIniGPV()

#
utilities_logMes('READ DRIVEN DOF')
mecaMAILx_ReadDrivenDof()

#
CSxxx_LoadTactors()
ASpxx_LoadTactors()

#
utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
CSASp_ReadIniVlocRloc()

### ecriture paranoiaque du modele ###
utilities_logMes('WRITE BODIES')
overall_WriteBodies()
MAILx_WriteBodies()

utilities_logMes('WRITE MODELS')
models_WriteModels()

utilities_logMes('WRITE BEHAVIOURS')
bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('WRITE DRIVEN DOF')
overall_WriteDrivenDof()
mecaMAILx_WriteDrivenDof()

utilities_logMes('WRITE LAST DOF')
TimeEvolution_WriteLastDof()
mecaMAILx_WriteLastDof()

utilities_logMes('WRITE LAST Vloc_Rloc')
TimeEvolution_WriteLastVlocRloc()
CSASp_WriteLastVlocRloc()

### post3D ##
OpenDisplayFiles()
#OpenPostproFiles()

utilities_logMes('COMPUTE MASS')
mecaMAILx_ComputeMass()

utilities_logMes('COMPUTE STIFFNESS')
mecaMAILx_ComputeBulk()

# precondensation
mecaMAILx_SetPreconAllBodies()
CSxxx_PushPreconNodes()
ASpxx_PushPreconNodes()
mecaMAILx_ComputePreconW()

mecaMAILx_AssembKT()


for k in xrange(1,nb_steps+1,1):
   #
   utilities_logMes('INCREMENT STEP')
   TimeEvolution_IncrementStep()
   mecaMAILx_IncrementStep()
   
   utilities_logMes('DISPLAY TIMES')
   TimeEvolution_DisplayStep()

   utilities_logMes('COMPUTE Fext')
   mecaMAILx_ComputeFext()

   utilities_logMes('COMPUTE Fint')
   mecaMAILx_ComputeBulk()

   utilities_logMes('ASSEMBLAGE')
   mecaMAILx_AssembRHS()

   utilities_logMes('COMPUTE Free Vlocy')
   mecaMAILx_ComputeFreeVelocity()
   #
   utilities_logMes('SELECT PROX TACTORS')
   overall_SelectProxTactors(freq_detect)
   CSASp_SelectProxTactors()
   #
   CSASp_RecupRloc()

   utilities_logMes('RESOLUTION' )
    
   #nlgs_3D_ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
   SiconosNumerics_ExSolver()
   CSASp_StockRloc()
   #
   utilities_logMes('COMPUTE DOF, FIELDS, etc.')
   mecaMAILx_ComputeDof()
   mecaMAILx_ComputeField()
   #
   utilities_logMes('UPDATE DOF, FIELDS')
   TimeEvolution_UpdateStep()
   mecaMAILx_UpdateDof()
   mecaMAILx_UpdateBulk()
   #
   utilities_logMes('WRITE out DOF')
   TimeEvolution_WriteOutDof(freq_write)
   mecaMAILx_WriteOutDof()
   #
   utilities_logMes('WRITE out Rloc')
   TimeEvolution_WriteOutVlocRloc(freq_write)
   CSASp_WriteOutVlocRloc()

   ### post3D ###
   WriteDisplayFiles(freq_display,ref_radius)
   #WritePostproFiles()

   ### gestion des writeout ###
   overall_CleanWriteOutFlags()

### postpro ###
CloseDisplayFiles()
#ClosePostproFiles()
