import os,sys

import numpy
import math
import pickle


from siconos.mechanics.contact_detection.tools import Contactor
from siconos.mechanics.contact_detection.convexhull import ConvexHull

#from siconos.io.mechanics_io import Hdf5

from mechanics_io import Hdf5
import siconos.numerics as Numerics



print 'load body_collection'
f = open('body_collection_1000.dict', 'r')
body_collection = pickle.load(f)
f.close()

scale = body_collection['scale']
print ('scale', scale)

step=50000
hstep=1e-4
itermax=1000
dump_probability = .02
theta=0.50001
tolerance=1e-4


solver=Numerics.SICONOS_FRICTION_3D_NSGS

fileName = "TolvaChute_1000"
title = "TolvaChute"
description = """
Tolva Chute of 1000 particles with Bullet collision detection
Moreau TimeStepping: h={0}, theta = {1}
One Step non smooth problem: {2}, maxiter={3}, tol={4}
""".format(hstep, theta, Numerics.idToName(solver),
           itermax,
           tolerance)

mathInfo = ""

from FrictionContactTrace import FrictionContactTraceParams
friction_contact_trace_params = FrictionContactTraceParams(dump_itermax=200, dump_probability=None,
                                                           fileName=fileName, title =title,
                                                           description = description, mathInfo= mathInfo)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+', collision_margin=0.01, io_filename= 'chute_scale_scene_1000.hdf5') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
  io.run(with_timer=False,
         time_stepping=None,
         space_filter=None,
         body_class=None,
         shape_class=None,
         face_class=None,
         edge_class=None,
         gravity_scale=1.0/scale,
         t0=0,
         T=step*hstep,
         h=hstep,
         multipoints_iterations=True,
         theta=theta,
         Newton_max_iter=1,
         set_external_forces=None,
         solver=solver,
         itermax=itermax,
         tolerance=tolerance,
         numerics_verbose=False,
         violation_verbose=True,
         output_frequency=100,
         friction_contact_trace=True,
         friction_contact_trace_params=friction_contact_trace_params
  )
