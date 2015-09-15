# import des modules
import os, sys
import math, numpy, copy

from pylmgc90.pre_lmgc import *

if (not os.path.isdir('./DATBOX')):
  os.mkdir('./DATBOX')

# definition des conteneurs:
#   * de corps
bodies = avatars()
#   * de materiaux
mats = materials()
#   * de lois de contacts
tacts = tact_behavs()
#   * de tables de visibilite
svs = see_tables()

# exemple 3D
dim = 3

# definition d'un modele rigide
mR3D = model(name='rigid', type='MECAx', element='Rxx3D', dimension=dim)

# definition du materiau pour les blocs
stone = material(name='STONE', type='RIGID', density=1.)
# ajout des materiaux dans le conteneur
mats.addMaterial(stone)

# ajout des marches dans le conteneur
blk = lecture(name='gmsh/bof.msh', dim=dim)
#reorientSurfacicElements(blk)
block1 = surfacicMeshToRigid3D(surfacic_mesh=blk, model=mR3D, material=stone, color = 'BLUEx')
block1.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
bodies.addAvatar(block1)

for i in xrange(0,9):
  blk = lecture(name='gmsh/bof.msh', dim=dim)
  #reorientSurfacicElements(blk)
  block = surfacicMeshToRigid3D(surfacic_mesh=blk, model=mR3D, material=stone, color = 'BLUEx')
  #block2.rotate(type='axis', alpha=math.pi, axis=[0,1,0,], center=block2.nodes[1].coor)
  block.translate(dz=(i+1)*0.5)
  bodies.addAvatar(block)

# definition d'une loi de contact frottant, avec pre-gap
iqsg0=tact_behav(name='iqsg0', type='IQS_CLB', fric=0.75)
tacts.addBehav(iqsg0)

# definition d'une table de visibilite pour le
# contact polyedre-polyedre
sv = see_table(CorpsCandidat='RBDY3', candidat='POLYR',
colorCandidat='BLUEx', behav=iqsg0,
CorpsAntagoniste='RBDY3', antagoniste='POLYR',
colorAntagoniste='BLUEx', alert=.1)
# ajout de la table de visibilite dans le conteneur
# de tables de visibilite
svs.addSeeTable(sv)

# ecriture des fichiers de donnees pour LMGC90
#   * les materiaux
writeBulkBehav(mats, chemin='./DATBOX/', dim=dim, gravy=[0., 0., -9.81])
#   * les corps
writeBodies(bodies, chemin='./DATBOX/')
#   * les deplacements et vitesses initiaux
writeDofIni(bodies, chemin='./DATBOX/')
#   * les conditions aux limites
writeDrvDof(bodies, chemin='./DATBOX/')
#   * les contacts envisages
writeTactBehav(tacts, svs, chemin='./DATBOX/')
#   * les contacts initiaux (fichier vide)
writeVlocRlocIni(chemin='./DATBOX/')
