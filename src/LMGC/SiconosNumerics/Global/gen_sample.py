import sys
import copy

import os

if (not os.path.isdir('./DATBOX')):
  os.mkdir('./DATBOX')

from pylmgc90.pre_lmgc import *

# definition des conteneurs:
#   * de corps
bodies = avatars()
#   * de modeles
mods = models()
#   * de materiaux
mats = materials()
#   * pour les tables de visibilite
svs = see_tables()
#   * pour les lois de contact
tacts = tact_behavs()

# on se place en 3D
dim = 3

# definition d'un modele elastique pour les cubes, mailles en hexaedres 
m3Dl = model(name='M3DH8', type='MECAx', element='H8xxx', dimension=3, 
     external_model='yes__', kinematic='small', material='elas_', 
     anisotropy='iso__', mass_storage='coher')
# on ajoute le modele dans le conteneur
mods.addModel(m3Dl)

# on definit le materiau constitutif des cubes
stone = material(name='stone', type='ELAS', density=2750., elas='standard',
   anisotropy='isotropic', young=7.e10, nu=0.2)  
# on l'ajoute dans le contenaur
mats.addMaterial(stone)

# construction des maillages :

nb_e = 100

# on contruit le maillage du cube
mesh_cube=buildMeshH8(x0=0., y0=0., z0=0., lx=1., ly=1., lz=1., nb_elem_x=nb_e, nb_elem_y=nb_e, nb_elem_z=nb_e)

# on construit un avatar deformable pour le premier cube
cube=buildMeshedAvatar(mesh=mesh_cube, model=m3Dl, material=stone)

## contacteurs :
#   * antagonistes sur la face du haut
cube.addContactors(group='up', type='ASpx4', color='BLEUx')
cube.imposeDrivenDof(group='down',component=[1,2,3],dofty='vlocy')

# ajout du cube dans le conteneur de corps
bodies += cube

# on copie le maillage pour construire le deuxieme cube
mesh_cube_2 = copy.deepcopy(mesh_cube)

# on construit un avatar deformable pour le deuxieme cube
cube_2=buildMeshedAvatar(mesh=mesh_cube_2, model=m3Dl, material=stone)

# contacteurs :
#   * candidats sur la face du bas

cube_2.addContactors(group='down', type='CSpx4', color='BLEUx')

# on place le deuxieme cube sur le premier
cube_2.translate(dz=1.)

# ajout du cube dans le conteneur de corps
bodies += cube_2


visuAvatars(bodies)

# gestion des interactions :
#   * declaration des lois
##       - entre cubes
lcsas=tact_behav('gapc0','GAP_SGR_CLB',fric=0.3)
tacts+=lcsas

##   * declaration des tables de visibilite
#       - entre cubes
svcsas = see_table(CorpsCandidat='MAILx', candidat='CSxxx',
   colorCandidat='BLEUx', behav=lcsas, CorpsAntagoniste='MAILx', 
   antagoniste='ASpxx', colorAntagoniste='BLEUx', alert=0.1, halo=0.5)
svs+=svcsas

# Ecriture des fichiers pour LMGC
writeBodies(bodies, chemin='./DATBOX/')
writeDofIni(bodies, chemin='./DATBOX/')
writeDrvDof(bodies, chemin='./DATBOX/')
writeModels(mods,chemin='./DATBOX/')
writeGPVIni(bodies, chemin='./DATBOX/')
writeBulkBehav(mats, chemin='./DATBOX/', dim=dim)
writeTactBehav(tacts,svs,chemin='DATBOX/')
writeVlocRlocIni(chemin='DATBOX/')
