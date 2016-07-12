import os,sys

import numpy
import math
import pickle

import random
with_lmgc = False
if with_lmgc:
  from pylmgc90.pre_lmgc import *

from siconos.mechanics.contact_detection.tools import Contactor
from siconos.mechanics.contact_detection.convexhull import ConvexHull

from siconos.io.mechanics_io import Hdf5
#sys.path.append('../..')
#from mechanics_io import Hdf5
import siconos.numerics as Numerics

if with_lmgc:
  if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')
    
    # WARNING : in 3D by default z-axis is upward
    # this is very important to direct PLANx objects

    dim = 3
    
    bodies = avatars()
    mat = materials()
    svs = see_tables()
    tacts = tact_behavs()


    #create materials
    tdur = material(name='TDURx',type='RIGID',density=7000.)
    pdur = material(name='MOUxx',type='RIGID',density=2300.)
    mat.addMaterial(tdur,pdur)


    # create a model of rigid
    mod = model(name='rigid', type='MECAx', element='Rxx3D', dimension=dim)

unscaled_cube_size=0.05
unscaled_density=2500

scale=1.0/unscaled_cube_size*10.0

cube_size = unscaled_cube_size*scale
density =  unscaled_density/(scale**3)


box_height = 3.683*scale
box_length = 6.900*scale
right_length= 4.370*scale
box_width  = 3.430*scale
rear_height = 1.200*scale
rear_length =1.770*scale
bottom_width =1.280*scale
output_height = .6*scale

plan_thickness = 0.05*scale

body_collection={}
body_collection['scale']= scale

body_collection['plan_id']= {}
id_plan=0

def normal_plane(p1,p2,p3):

  x1=p1[0]
  y1=p1[1]
  z1=p1[2]
  x2=p2[0]
  y2=p2[1]
  z2=p2[2]
  x3=p3[0]
  y3=p3[1]
  z3=p3[2]

  vector1 = [x2 - x1, y2 - y1, z2 - z1]
  vector2 = [x3 - x1, y3 - y1, z3 - z1]

  cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1], -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]), vector1[0] * vector2[1] - vector1[1] * vector2[0]]

  a = cross_product[0]
  b = cross_product[1]
  c = cross_product[2]
  d = - (cross_product[0] * x1 + cross_product[1] * y1 + cross_product[2] * z1)

  return numpy.array([a,b,c])/numpy.linalg.norm([a,b,c])




#create some bodies

# Creation of the hdf5 file for input/output
with Hdf5() as io:

  angle= 0.0 ; math.pi/4.0
  translation_plan = [0.0, 0.0 ,0.0]
  orientation_plan = [math.cos(angle/2.0), 0.0, math.sin(angle/2.0), 0.0]
  orientation_plan = [1.0, 0.0, 0,0, 0.0]
  #orientation_plan = [0.0, 1.0, 0,0, 0.0] #--> cos(\theta/2.0) = 0, sin(\theta/2.0) =1 ==> \theta = pi
  
  ######### left_up
  id_plan=id_plan+1
  body_collection['plan_id']["left_up"]=id_plan
  #left = rigidPlan(axe1=box_length/2.0 , axe2=box_height/2.0 , axe3=plan_thickness,
  #                 center=[0.,0.,0.],material=tdur, model=mod, color='VERTx')

  v1 = numpy.array([0, 0 , box_height-rear_height])
  v2 = numpy.array([right_length-right_length*rear_height/box_height,0.0, rear_height])
  v3 = numpy.array([ box_length, 0,rear_height])

  left_up_normal = normal_plane(v1,v2,v3)
  #print('left_up_normal=', left_up_normal)
  v1_extruded = v1 + numpy.dot(plan_thickness,left_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,left_up_normal)
  v3_extruded = v3 + numpy.dot(plan_thickness,left_up_normal)

  left_up_vertices=numpy.array([v1,v2,v3,v1_extruded,v2_extruded,v3_extruded])
  #print left_up_vertices

  v1 = numpy.array([0, 0 , box_height])
  v2 = numpy.array([right_length-right_length*rear_height/box_height,0.0, rear_height])
  v3 = numpy.array([ box_length, 0,rear_height])
  v1_extruded = v1 + numpy.dot(plan_thickness,left_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,left_up_normal)
  v3_extruded = v3 + numpy.dot(plan_thickness,left_up_normal)

  left_up_vertices=numpy.array([v1,v2,v3,v1_extruded,v2_extruded,v3_extruded])
  #print left_up_vertices
  if with_lmgc:
    left_up=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                            nb_vertices=len(left_up_vertices), vertices=left_up_vertices, faces=None, radius=1., tol=0., number=None)


    bodies.addAvatar(left_up)

  io.addConvexShape('Left_up',left_up_vertices )
  io.addObject('left_up', [Contactor('Left_up')],
               translation=translation_plan, orientation=orientation_plan)


  
  ######### left_middle
  id_plan=id_plan+1
  body_collection['plan_id']["left_middle"]=id_plan

  v4 = numpy.array([right_length,bottom_width, 0.0])
  v5 = numpy.array([ box_length-rear_length, bottom_width,0.0])

  left_middle_normal = normal_plane(v2,v4,v3)
  #print('left_middle_normal=', left_middle_normal)

  v4_extruded = v4 + numpy.dot(plan_thickness,left_middle_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,left_middle_normal)

  left_middle_vertices=numpy.array([v2,v3,v4,v5,v2_extruded,v3_extruded,v4_extruded,v5_extruded])
  #print left_middle_vertices
  if with_lmgc:
    left_middle=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                                nb_vertices=len(left_middle_vertices), vertices=left_middle_vertices,
                                faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(left_middle)
    
  io.addConvexShape('Left_middle',left_middle_vertices )
  io.addObject('left_middle', [Contactor('Left_middle')],
               translation=translation_plan, orientation=orientation_plan)

  ######### left_down
  id_plan=id_plan+1
  body_collection['plan_id']["left_down"]=id_plan

  v6 = numpy.array([right_length,box_width, -output_height])
  v7 = numpy.array([box_length-rear_length, box_width,-output_height])

  left_down_normal = normal_plane(v4,v6,v5)
  #print('left_down_normal=', left_down_normal)

  v6_extruded = v6 - [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,left_down_normal)
  v7_extruded = v7 + [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,left_down_normal)

  left_down_vertices=numpy.array([v4-[plan_thickness, 0.0 ,0.],v5+[plan_thickness, 0.0 ,0.],
                                  v6-[plan_thickness, 0.0 ,0.],v7+[plan_thickness, 0.0 ,0.],
                                  v4_extruded-[plan_thickness, 0.0 ,0.],v5_extruded+[plan_thickness, 0.0 ,0.],
                                  v6_extruded,v7_extruded])
  #print left_down_vertices
  if with_lmgc:
    left_down=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                              nb_vertices=len(left_down_vertices), vertices=left_down_vertices,
                              faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(left_down)
    
  io.addConvexShape('Left_down',left_down_vertices )
  io.addObject('left_udown', [Contactor('Left_down')],
               translation=translation_plan, orientation=orientation_plan)

  ######### right_up
  id_plan=id_plan+1
  body_collection['plan_id']["right_up"]=id_plan

  v8 = numpy.array([0, box_width , box_height])
  v9 = numpy.array([ box_length, box_width,rear_height])

  v10 = numpy.array([ box_length-rear_length, box_width,0.0])
  v11 =  numpy.array([right_length,box_width, 0.0])


  right_up_normal = normal_plane(v8,v9,v10)
  #print('right_up_normal=', right_up_normal)

  v8_extruded = v8 + numpy.dot(plan_thickness,right_up_normal)
  v9_extruded = v9 + numpy.dot(plan_thickness,right_up_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,right_up_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,right_up_normal)

  right_up_vertices=numpy.array([v8,v9,v10,v11,v8_extruded,v9_extruded,v10_extruded,v11_extruded])
  #print right_up_vertices
  if with_lmgc:
    right_up=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                             nb_vertices=len(right_up_vertices), vertices=right_up_vertices,
                             faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(right_up)

  io.addConvexShape('Right_up',right_up_vertices )
  io.addObject('right_up', [Contactor('Right_up')],
               translation=translation_plan, orientation=orientation_plan)

  ######### rear_up
  id_plan=id_plan+1
  body_collection['plan_id']["rear_up"]=id_plan


  rear_up_normal = normal_plane(v1,v8,v4)
  #print('rear_up_normal=', rear_up_normal)

  v1_extruded = v1 + numpy.dot(plan_thickness,rear_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,rear_up_normal)
  v8_extruded = v8 + numpy.dot(plan_thickness,rear_up_normal)
  v4_extruded = v4 + numpy.dot(plan_thickness,rear_up_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,rear_up_normal)

  rear_up_vertices=numpy.array([v1-[0.0,plan_thickness,0.0],v2-[0.0,plan_thickness,0.0],
                                v8+[0.0,plan_thickness,0.0],v4-[0.0,plan_thickness,0.0],v11+[0.0,plan_thickness,0.0],
                                v1_extruded-[0.0,plan_thickness,0.0],v2_extruded-[0.0,plan_thickness,0.0],
                                v8_extruded+[0.0,plan_thickness,0.0],v4_extruded-[0.0,plan_thickness,0.0],
                                v11_extruded+[0.0,plan_thickness,0.0]])
  #print rear_up_vertices
  if with_lmgc:
    rear_up=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                            nb_vertices=len(rear_up_vertices), vertices=rear_up_vertices,
                            faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(rear_up)

  io.addConvexShape('Rear_up',rear_up_vertices )
  io.addObject('rear_up', [Contactor('Rear_up')],
               translation=translation_plan, orientation=orientation_plan)


  ######### rear_down
  id_plan=id_plan+1
  body_collection['plan_id']["rear_down"]=id_plan
  #v12 = numpy.array([ box_length-rear_length, box_width,-output_height])
  #v13 = numpy.array([right_length,box_width, -output_height])


  rear_down_normal = normal_plane(v4,v11,v6)
  #print('rear_down_normal=', rear_down_normal)

  v4_extruded = v4 + numpy.dot(plan_thickness,rear_down_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,rear_down_normal)
  v6_extruded = v6 + numpy.dot(plan_thickness,rear_down_normal)

  rear_down_vertices=numpy.array([v4,v11,v6,v4_extruded,v11_extruded,v6_extruded])
  #print rear_down_vertices
  if with_lmgc:
    rear_down=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                              nb_vertices=len(rear_down_vertices), vertices=rear_down_vertices,
                              faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(rear_down)

  io.addConvexShape('Rear_down',rear_down_vertices )
  io.addObject('rear_down', [Contactor('Rear_down')],
               translation=translation_plan, orientation=orientation_plan)

  ######### front_up
  id_plan=id_plan+1
  body_collection['plan_id']["front_up"]=id_plan


  front_up_normal = normal_plane(v3,v5,v9)
  #print('front_up_normal=', front_up_normal)

  v3_extruded = v3 + numpy.dot(plan_thickness,front_up_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,front_up_normal)
  v9_extruded = v9 + numpy.dot(plan_thickness,front_up_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,front_up_normal)

  front_up_vertices=numpy.array([v3-[0.0,plan_thickness,0.0],v5-[0.0,plan_thickness,0.0],
                                 v9+[0.0,plan_thickness,0.0],v10+[0.0,plan_thickness,0.0],
                                 v3_extruded-[0.0,plan_thickness,0.0],v5_extruded-[0.0,plan_thickness,0.0],
                                 v9_extruded+[0.0,plan_thickness,0.0],v10_extruded+[0.0,plan_thickness,0.0]])
  #print front_up_vertices
  if with_lmgc:
    front_up=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                             nb_vertices=len(front_up_vertices), vertices=front_up_vertices,
                             faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(front_up)
  io.addConvexShape('Front_up',front_up_vertices )
  io.addObject('front_up', [Contactor('Front_up')],
               translation=translation_plan, orientation=orientation_plan)

  ######### front_down
  id_plan=id_plan+1
  body_collection['plan_id']["front_down"]=id_plan


  front_down_normal = normal_plane(v5,v7,v10)
  #print('front_down_normal=', front_down_normal)

  v7_extruded = v7 + numpy.dot(plan_thickness,front_down_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,front_down_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,front_down_normal)

  front_down_vertices=numpy.array([v5,v7,v10,v5_extruded,v7_extruded,v10_extruded])
  #print front_down_vertices
  if with_lmgc:
    front_down=rigidPolyhedron(model=mod, material=tdur, center=numpy.zeros(3), color='VERTx', generation_type='vertices',
                               nb_vertices=len(front_down_vertices), vertices=front_down_vertices,
                               faces=None, radius=1., tol=0., number=None)
    bodies.addAvatar(front_down)
  io.addConvexShape('Front_down',front_down_vertices )
  io.addObject('front_down', [Contactor('Front_down')],
               translation=translation_plan, orientation=orientation_plan)

  if with_lmgc:
  
    sizeofparticles_min=0.05
    sizeofparticles_max=0.15
    nparticles= 1

    # random distribution  of radii
    radii=granulo_Random(nparticles, sizeofparticles_min, sizeofparticles_max)
    radius_min=min(radii)
    radius_max=max(radii)

    ratio =0.4
    lx = box_length *ratio
    ly = box_width * ratio
    lz = box_height*5

    [nb_remaining_particles, coor_particles]=depositInBox3D(radii, lx, ly, lz)

    if (nb_remaining_particles < nparticles):
       print "Warning: granulometry changed, since some particles were removed!"
       print "nparticles is now : nb_remaining_particles "
       nparticles= nb_remaining_particles

    coor_particles.resize(nb_remaining_particles*dim)
    coor_particles.shape=[nparticles,dim]

    x0 = lx /2.0 *1.2
    y0 = box_width/2.0
    z0 = box_height

    for i in range(nparticles):
      coor_particles[i,0] = coor_particles[i,0]+x0
      coor_particles[i,1] = coor_particles[i,1]+y0
      coor_particles[i,2] = coor_particles[i,2]+z0

    particles=[]
    body_collection['particles_id']= []
    for i in range(nparticles):
    #  particles.append(rigidSphere(r=radii[i],center=coor_particles[i,:],model=mod,material=pdur,color='BLEUx'))
      print("Creation particules number", i)
      particles.append(rigidPolyhedron(radius=radii[i],center=coor_particles[i,:],
                                       model=mod,material=pdur,color='BLEUx',
                                       nb_vertices=6,generation_type='random'))

      bodies += particles[i]
      body_collection['particles_id'].append(particles[i].number+1)


    


  n_cube=10
  n_row=10
  n_col=10
  x_shift=3.0

  angle =0.0
  orientation_cube = [math.cos(angle/2.0), 0.0, math.sin(angle/2.0), 0.0]
  orientation_cube = [1.0, 0.0, 0,0, 0.0]
  count =0
  for i in range(n_row):
    for j in range(n_col):
      for n in range(n_cube):
        count +=1
        print("Creation particules number", i, j, n, count)
        cube_size_rand = cube_size*(1.0 + 0.5*( random.random()-1.0))
        # Definition of a cube as a convex shape
        # cubes
        # vertices =[(-cube_size_rand, cube_size_rand, -cube_size_rand),
        #            (-cube_size_rand, -cube_size_rand, -cube_size_rand),
        #            (-cube_size_rand, -cube_size_rand, cube_size_rand),
        #            (-cube_size_rand, cube_size_rand, cube_size_rand),
        #            (cube_size_rand, cube_size_rand, cube_size_rand),
        #            (cube_size_rand, cube_size_rand, -cube_size_rand),
        #            (cube_size_rand, -cube_size_rand, -cube_size_rand),
        #            (cube_size_rand, -cube_size_rand, cube_size_rand)])
        # irregular polyhedron
        vertices=[ (-cube_size_rand, cube_size_rand, -cube_size_rand),
                   (-cube_size_rand, -cube_size_rand, -cube_size_rand),
                   (-cube_size_rand, -cube_size_rand, cube_size_rand),
                   (cube_size_rand, cube_size_rand, cube_size_rand),
                   (cube_size_rand, cube_size_rand, -cube_size_rand),
                   (cube_size_rand, -cube_size_rand, -cube_size_rand),
                   (0.0,0.0, -2.0*cube_size_rand)]
        #print('Polyhedron vertices', vertices)
 
        ch = ConvexHull(vertices)
        cm = ch.centroid()
        #print('cm', cm)
 
        # correction of vertices such that 0 is the centroid
        vertices = numpy.array(vertices)[:]-cm[:]
        #print('Corrected polyhedron vertices', vertices)
        ch = ConvexHull(vertices)
        cm = ch.centroid()
        #print('cm', cm)

        # Definition of a polyhedron as a convex shape
        io.addConvexShape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j),vertices )

        # computation of inertia and volume
        inertia,volume=ch.inertia(ch.centroid())

        # initial translation
        trans=[box_width/5.0+i*x_shift*cube_size,
               x_shift*(j+2)*cube_size,
               box_height+cube_size*x_shift*n]

        io.addObject('cube'+str(n)+'_'+str(i)+'_'+str(j),
                     [Contactor('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))],
                     translation=trans, orientation = orientation_cube,
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=volume*density, inertia = inertia*density)
        #raw_input()
  if with_lmgc:
    #b=tact_behav('rstc1','RST_CLB',fric=0.3,rstn=0.9,rstt=0.5)
    #tacts+=b
    b=tact_behav(name='iqsc0',type='IQS_CLB',fric=0.3)
    tacts+=b

    #interactions
    #sv = see_table(CorpsCandidat='RBDY3',candidat='SPHER',colorCandidat='BLEUx',behav='iqsc0',
    sv = see_table(CorpsCandidat='RBDY3',candidat='POLYR',colorCandidat='BLEUx',behav=b,
                   CorpsAntagoniste='RBDY3',antagoniste='POLYR',colorAntagoniste='VERTx',alert=.01)
    svs+=sv
    #sv = see_table(CorpsCandidat='RBDY3',candidat='SPHER',colorCandidat='BLEUx',behav='iqsc0',
    sv = see_table(CorpsCandidat='RBDY3',candidat='POLYR',colorCandidat='BLEUx',behav=b,
                   CorpsAntagoniste='RBDY3',antagoniste='POLYR',colorAntagoniste='BLEUx',alert=.01)
    svs+=sv


  # Definition of a non smooth law. As no group ids are specified it
  # is between contactors of group id 0.
  io.addNewtonImpactFrictionNSL('contact', mu=0.3)

  print body_collection

  f = open('body_collection.dict', 'w')
  pickle.dump(body_collection,f)
  f.close()

