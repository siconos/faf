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
from __future__ import print_function

import Siconos.Numerics as N
import Siconos.Kernel as K
import Siconos.FCLib as F
import imp
import sys
import time
import getopt


from Siconos.Kernel import \
     Model, Moreau, TimeDiscretisation,\
     FrictionContact, NewtonImpactFrictionNSL, BlockCSRMatrix

from Siconos.Mechanics.ContactDetection.Bullet import IO, \
    btConvexHullShape, btCollisionObject, \
    btBoxShape, btQuaternion, btTransform, btConeShape, \
    BulletSpaceFilter, cast_BulletR, \
    BulletWeightedShape, BulletDS, BulletTimeStepping

from Siconos.Kernel import GlobalFrictionContactProblem

class Timer():
    def __init__(self):
        self._t0 = time.clock()

    def elapsed(self):
        return time.clock() - self._t0

    def update(self):
        self._t0 = time.clock()

with_timer = False


def usage():
    print('{0}: Usage'.format(sys.argv[0]))


def log(fun):
    if with_timer:
        t = Timer()

        def logged(*args):
            t.update()
            print('{0} ...'.format(fun.__name__), end='')
            fun(*args)
            print('..... {0} s'.format(t.elapsed()))
        return logged
    else:
        def silent(*args):
            fun(*args)
        return silent

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'timer'])
except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)
for o, a in opts:
    if o == '--help':
        usage()
        sys.exit(0)
    elif o == '--timer':
        with_timer = True


def warn(msg):
    sys.stderr.write('{0}: {1}'.format(sys.argv[0], msg))


def usage():
    pass


try:
    imp.load_source('params', 'params.py')
except IOError as e:
    warn('I need a params.py file')
    usage()
    exit(1)
import params


t0 = params.t0       # start time
T = params.T      # end time
h = params.h    # time step
g = params.g     # gravity

theta = params.theta  # theta scheme
mu = params.mu
dump_itermax = params.dump_itermax
itermax = params.itermax
tolerance = params.tolerance

#
# Model
#
model = Model(t0, T)

#
# Simulation
#
# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0., 0., mu, 3)

# (1) OneStepIntegrators
osi = Moreau(theta)

static_cobjs = []


# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
solver = params.solver


class FrictionContactTrace(FrictionContact):

    def __init__(self, dim, solver, maxiter):
        self._maxiter = maxiter
        super(FrictionContactTrace, self).__init__(dim, solver)

    def compute(self,time):
        info = 0
        cont = self.preCompute(time)

        if (not cont):
            return 0

        if (self.indexSetLevel() == 999):
            return 0

        self.updateMu()

        if self.getSizeOutput() != 0:

            #            M = BlockCSRMatrix()
            #M.fillM(model.nonSmoothDynamicalSystem().topology().indexSet(1))
            #M.convert()

            #H = BlockCSRMatrix()
            #H.fillH(model.nonSmoothDynamicalSystem().topology().indexSet(1))
            #H.convert()

            #t = GlobalFrictionContactProblem() 

            #t.M = M.getNumericsMatSparse()

            w_backup = self.w().copy()
            z_backup = self.z().copy()

            info = self.solve()
            SO = self.numericsSolverOptions()
            if SO.iparam[7] > self._maxiter:
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
                guess = F.fclib_solution()
                guess.u = w_backup
                guess.r = z_backup
                F.fclib_write_guesses([guess], filename)

                solution = F.fclib_solution()
                solution.u = self.w()
                solution.z = self.z()
                F.fclib_write_solution(solution, filename)

            self.postCompute()

        return info

osnspb = FrictionContactTrace(3, solver, dump_itermax)

osnspb.numericsSolverOptions().iparam[0] = itermax
osnspb.numericsSolverOptions().dparam[0] = tolerance
osnspb.setMaxSize(16384)
osnspb.setMStorageType(1)
#osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)

nopts = osnspb.numericsOptions()
nopts.verboseMode = 0
nopts.outputMode = 0  # 3 | N.OUTPUT_ON_ERROR
nopts.fileName = params.fileName
nopts.title = params.title
nopts.description = params.description
nopts.mathInfo = params.mathInfo

# (5) broadphase contact detection
broadphase = BulletSpaceFilter(model, nslaw)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = BulletTimeStepping(timedisc, broadphase)
simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)
simulation.setNewtonMaxIteration(params.NewtonMaxIter)

k = 1

# time loop


with IO.Hdf5(broadphase, osi) as io:

    model.initialize(simulation)

    if 'initialize' in dir(params):
        params.initialize(model)

    io.outputStaticObjects()
    io.outputDynamicObjects()
    while(simulation.hasNextEvent()):

        print ('step', k)

        log(broadphase.buildInteractions)(model.currentTime())

        log(simulation.computeOneStep)()

        log(io.outputDynamicObjects)()

        log(io.outputContactForces)()

        log(io.outputSolverInfos)()

        log(simulation.nextStep)()

        log(io._out.flush)()

        print ('')
        k += 1
