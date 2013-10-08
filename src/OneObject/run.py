#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2013.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY, siconos-team@lists.gforge.fr
#

import Siconos.Numerics as N
import Siconos.Kernel as K
import Siconos.FCLib as F

from Siconos.Kernel import \
     Model, Moreau, TimeDiscretisation,\
     FrictionContact, NewtonImpactFrictionNSL

from Siconos.Mechanics.ContactDetection.Bullet import IO, \
    btConvexHullShape, btCollisionObject, \
    btBoxShape, btQuaternion, btTransform, btConeShape, \
    BulletSpaceFilter, cast_BulletR, \
    BulletWeightedShape, BulletDS, BulletTimeStepping


t0 = 0       # start time
T = 10.      # end time
h = 0.0005    # time step

g = 9.81     # gravity

theta = 0.51  # theta scheme

#
# Model
#
model = Model(t0, T)

#
# Simulation
#
# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0., 0., 0.7, 3)

# (1) OneStepIntegrators
osi = Moreau(theta)

static_cobjs = []


# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
solver = N.SICONOS_FRICTION_3D_NSGS


class FrictionContactTrace(FrictionContact):

    def compute(self,time):
        info = 0
        cont = self.preCompute(time)

        if (not cont):
            return 0

        if (self.indexSetLevel() == 999):
            return 0

        self.updateMu()

        if self.getSizeOutput() != 0:

            w_backup = self.w()
            z_backup = self.z()

            info = self.solve()
            SO = self.numericsSolverOptions()
            if SO.iparam[7] > 1000:
                lnopts = self.numericsOptions()
                problem = self.frictionContactProblem()
                filename = "{0}-i{1}-{2}.hdf5".format(lnopts.fileName,
                                                      SO.iparam[7],
                                                      nopts.counter)
                lnopts.counter += 1
                K.frictionContact_fclib_write(problem,
                                              lnopts.title,
                                              lnopts.description,
                                              lnopts.mathInfo,
                                              filename)
                sol = F.fclib_solution()
                sol.u = w_backup
                sol.r = z_backup
                F.fclib_write_guesses([sol], filename)

            self.postCompute()

        return info

osnspb = FrictionContactTrace(3, solver)

osnspb.numericsSolverOptions().iparam[0] = 100000
osnspb.numericsSolverOptions().dparam[0] = 1e-8
osnspb.setMaxSize(16384)
osnspb.setMStorageType(1)
#osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)

nopts = osnspb.numericsOptions()
nopts.verboseMode = 0
nopts.outputMode = 0  # 3 | N.OUTPUT_ON_ERROR
nopts.fileName = "BoxesStack1"
nopts.title = "Boxes Stack"
nopts.description = """
Boxes (Cubes) stacking with Bullet collision detection
Moreau TimeStepping: h={0}, theta = {1}
One Step non smooth problem: {2}, maxiter={3}, tol={4}
""".format(h, theta, N.idToName(solver),
           osnspb.numericsSolverOptions().iparam[0],
           osnspb.numericsSolverOptions().dparam[0])
#nopts.mathInfo = ""

print osnspb.numericsOptions().outputMode

# (5) broadphase contact detection
broadphase = BulletSpaceFilter(model, nslaw)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = BulletTimeStepping(timedisc, broadphase)
simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)
simulation.setNewtonMaxIteration(20)

k = 1

# time loop


with IO.Dat(broadphase, osi) as io:

    model.initialize(simulation)

    io.outputStaticObjects()
    io.outputDynamicObjects()
    while(simulation.hasNextEvent()):

        print 'contact detection', k
        broadphase.buildInteractions(model.currentTime())

        print 'computeOneStep'
        simulation.computeOneStep()

        io.outputDynamicObjects()
        io.outputContactForces()
        io.outputSolverInfos()

        print 'nextStep'
        simulation.nextStep()
        k += 1
