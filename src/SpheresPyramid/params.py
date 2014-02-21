import Siconos.Numerics as Numerics

t0 = 0
T = 30
h = 0.0005
g = 9.81
theta = 0.50001
mu = 2.0
dump_itermax = 80
dump_probability = .05
itermax = 100000
NewtonMaxIter = 20
tolerance = 1e-8
solver = Numerics.SICONOS_FRICTION_3D_NSGS


import imp
try:
    imp.load_source('mkinput', 'mkinput.py')
except IOError as e:
    warn('I need a mkinput.py file')
    usage()
    exit(1)
import mkinput


fileName = "SpheresPyramid{0}".format(mkinput.N)
title = "SpheresPyramid with {0} levels"
description = """
Spheres pyramid under gravity on the ground with Bullet collision detection
Moreau TimeStepping: h={0}, theta = {1}
One Step non smooth problem: {2}, maxiter={3}, tol={4}
""".format(h, theta, Numerics.idToName(solver),
           itermax,
           tolerance)

mathInfo = ""

# if we want a shuffled NonsmoothGaussSeidel
#def initialize(model):
    #    model.simulation().oneStepNSProblem(0).numericsSolverOptions().iparam[9] = 1
