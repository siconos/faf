#!/usr/bin/env python
from getopt import getopt, GetoptError
import sys
import shlex

import vtk
from vtk.util import numpy_support

from math import atan2, pi
import bisect
from numpy.linalg import norm

import numpy

import random


def random_color():
    r = random.uniform(0.0, 1.0)
    g = random.uniform(0.0, 1.0)
    b = random.uniform(0.0, 1.0)
    return r, g, b


def axis_angle(q):
    w, v = q[0], q[1:]
    nv = norm(v)
    theta = 2 * atan2(nv, w)

    if nv != 0.:
        v = [iv / nv for iv in v]
    else:
        v = [0., 0., 0.]
    return v, theta

axis_anglev = numpy.vectorize(axis_angle)
transforms = dict()


def set_position(instance, q0, q1, q2, q3, q4, q5, q6):

    axis, angle = axis_angle((q3, q4, q5, q6))

    transforms[instance].Identity()
    transforms[instance].Translate(q0, q1, q2)
    transforms[instance].RotateWXYZ(angle * 180 / pi,
                                    axis[0],
                                    axis[1],
                                    axis[2])    

set_positionv = numpy.vectorize(set_position)


def usage():
    print """{0}
    """.format(sys.argv[0])

try:
    opts, args = getopt(sys.argv[1:], "", [])
except GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

    #ref_filename = args[0]
    #bind_filename = args[1]
    #pos_filename = args[2]
    #cf_filename = args[3]

ref_filename = 'ref.txt'
bind_filename = 'bindings.dat'
dpos_filename = 'dpos.dat'
spos_filename = 'spos.dat'
cf_filename = 'cf.dat'

with open(ref_filename, 'r') as ref_file:
    refs = shlex.split(ref_file.read())

shape = dict()
with open(bind_filename, 'r') as bind_file:
    for line in bind_file:
        obj_id, obj_shape = shlex.split(line)
        shape[int(obj_id)] = int(obj_shape)

pos = dict()
instances = set()

import numpy
spos_data = numpy.loadtxt(spos_filename, ndmin=2)
dpos_data = numpy.loadtxt(dpos_filename, ndmin=2)
print spos_data.shape
print dpos_data.shape

cf_data = numpy.loadtxt(cf_filename)

#def contact_point_reader():
#    global dos
#    dos.GetOutput().GetFieldData()
    

#contact_field = vtk.vtkPointData()
#c1 = vtk.vtkPolyData()

#contact_data = vtk.vtkDataObject()
contact_pos = vtk.vtkDataObjectToDataSetFilter()

contact_pos_force = vtk.vtkFieldDataToAttributeDataFilter()
contact_pos_norm = vtk.vtkFieldDataToAttributeDataFilter()

#contact_data.SetFieldData(contact_field)

#id_f = numpy.where(cf_data[:,0] == min(cf_data[:,0]))[0]

#cf_provider = vtk.vtkProgrammableSource()
keeper = []


class CFprov():

    def __init__(self, data):
        if len(data) > 0:
            self._data = data
        else:
            self._data = None

        if self._data is not None:
            self._time = min(self._data[:, 0])
        else:
            self._time = 0
        self._output = None

    def method(self):
        global keeper
        self._output = vtk.vtkPolyData()
        contact_field = vtk.vtkPointData()
        ind0 = bisect.bisect(self._data[:,0], self._time)
        if self._data is not None and ind0 > 0:

            id_f = numpy.where(self._data[:, 0] == self._data[ind0,0])[0]

            cp_at_time = self._data[id_f, 1:4].copy()
            cp = numpy_support.numpy_to_vtk(cp_at_time)
            cp.SetName('contactPositions')

            cn_at_time = - self._data[id_f, 4:7]
            cn = numpy_support.numpy_to_vtk(cn_at_time)
            cn.SetName('contactNormals')
            cf_at_time = self._data[id_f, 7:10].copy()

            cf = numpy_support.numpy_to_vtk(cf_at_time)

            cf.SetName('contactForces')

            contact_field.AddArray(cp)
            contact_field.AddArray(cn)
            contact_field.AddArray(cf)

            keeper = [cp_at_time, cf_at_time, cn_at_time]
            self._output.SetFieldData(contact_field)
        else:
            pass
        

cf_prov = CFprov(cf_data)

contact_pos.SetDataSetTypeToPolyData()
contact_pos.SetPointComponent(0, "contactPositions", 0)
contact_pos.SetPointComponent(1, "contactPositions", 1)
contact_pos.SetPointComponent(2, "contactPositions", 2)

contact_pos_force.SetInputConnection(contact_pos.GetOutputPort())
contact_pos_force.SetInputFieldToDataObjectField()
contact_pos_force.SetOutputAttributeDataToPointData()
contact_pos_force.SetVectorComponent(0, "contactForces", 0)
contact_pos_force.SetVectorComponent(1, "contactForces", 1)
contact_pos_force.SetVectorComponent(2, "contactForces", 2)

contact_pos_norm.SetInputConnection(contact_pos.GetOutputPort())
contact_pos_norm.SetInputFieldToDataObjectField()
contact_pos_norm.SetOutputAttributeDataToPointData()
contact_pos_norm.SetVectorComponent(0, "contactNormals", 0)
contact_pos_norm.SetVectorComponent(1, "contactNormals", 1)
contact_pos_norm.SetVectorComponent(2, "contactNormals", 2)

times = list(set(dpos_data[:,0]))
times.sort()

ndyna = len(numpy.where(dpos_data[:, 0] == times[0]))

nstatic = len(numpy.where(spos_data[:, 0] == times[0]))

instances =  set(dpos_data[:,1]).union(set(spos_data[:,1]))

cf_prov._time = min(times[:])
cf_prov.method()
contact_pos.SetInput(cf_prov._output)
contact_pos.Update()
contact_pos_force.Update()
contact_pos_norm.Update()

arrow = vtk.vtkArrowSource()
arrow.SetTipResolution(40)
arrow.SetShaftResolution(40)
cone = vtk.vtkConeSource()
cone.SetResolution(40)
cone.SetRadius(0.3)

sphere = vtk.vtkSphereSource()

arrow_glyph = vtk.vtkGlyph3D()
arrow_glyph.SetInputConnection(contact_pos_force.GetOutputPort())
arrow_glyph.SetSourceConnection(arrow.GetOutputPort())
arrow_glyph.ScalingOn()
arrow_glyph.SetScaleModeToScaleByVector()
arrow_glyph.SetRange(-.5, 2)
arrow_glyph.ClampingOn()
arrow_glyph.SetScaleFactor(4)
arrow_glyph.SetVectorModeToUseVector()

arrow_glyph.SetInputArrayToProcess(1, 0, 0, 0, 'contactForces')
arrow_glyph.SetInputArrayToProcess(3, 0, 0, 0, 'contactForces')
arrow_glyph.OrientOn()

gmapper = vtk.vtkPolyDataMapper()
gmapper.SetInputConnection(arrow_glyph.GetOutputPort())
gmapper.SetScalarModeToUsePointFieldData()
gmapper.SetColorModeToMapScalars()
gmapper.ScalarVisibilityOn()
gmapper.SelectColorArray('contactForces')
#gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

gactor = vtk.vtkActor()
gactor.SetMapper(gmapper)

transform = vtk.vtkTransform()
transform.Translate(-0.5,0.,0.)

cone_glyph = vtk.vtkGlyph3D()
cone_glyph.SetSourceTransform(transform)

cone_glyph.SetInputConnection(contact_pos_norm.GetOutputPort())
cone_glyph.SetSourceConnection(cone.GetOutputPort())
cone_glyph.ScalingOff()
#cone_glyph.SetScaleModeToScaleByVector()
#cone_glyph.SetRange(-0.5, 2)
cone_glyph.ClampingOn()
#cone_glyph.SetScaleFactor(4)
cone_glyph.SetVectorModeToUseVector()

cone_glyph.SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
cone_glyph.OrientOn()



cmapper = vtk.vtkPolyDataMapper()
cmapper.SetInputConnection(cone_glyph.GetOutputPort())

sphere_glyph = vtk.vtkGlyph3D()
sphere_glyph.SetInputConnection(contact_pos_norm.GetOutputPort())
sphere_glyph.SetSourceConnection(sphere.GetOutputPort())
sphere_glyph.ScalingOn()
#sphere_glyph.SetScaleModeToScaleByVector()
#sphere_glyph.SetRange(-0.5, 2)
#sphere_glyph.ClampingOn()
sphere_glyph.SetScaleFactor(.1)
#sphere_glyph.SetVectorModeToUseVector()

#sphere_glyph.SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
#sphere_glyph.OrientOn()

smapper = vtk.vtkPolyDataMapper()
smapper.SetInputConnection(sphere_glyph.GetOutputPort())

#cmapper.SetScalarModeToUsePointFieldData()
#cmapper.SetColorModeToMapScalars()
#cmapper.ScalarVisibilityOn()
#cmapper.SelectColorArray('contactNormals')
#gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

cactor = vtk.vtkActor()
cactor.GetProperty().SetOpacity(0.4)
cactor.GetProperty().SetColor(0,0,1)
cactor.SetMapper(cmapper)

sactor = vtk.vtkActor()
#sactor.GetProperty().SetOpacity(0.4)
sactor.GetProperty().SetColor(1,0,0)
sactor.SetMapper(smapper)

#with open(pos_filename, 'r') as pos_file:
#    for line in pos_file:
#        data = shlex.split(line)
#        time = float(data[0])
#        instance = int(data[1])
#        q0 = float(data[2])
#        q1 = float(data[3])
#        q2 = float(data[4])        
#        q3 = float(data[5])
#        q4 = float(data[6])
#        q5 = float(data[7])
#        q6 = float(data[8])
#        instances.add(instance)
#        if time in pos:
#            pos[time].append([instance, q0, q1, q2, q3, q4, q5, q6])
#        else:
#            pos[time] = [[instance, q0, q1, q2, q3, q4, q5, q6]]

renderer = vtk.vtkRenderer()
renderer_window = vtk.vtkRenderWindow()
interactor_renderer = vtk.vtkRenderWindowInteractor()

#camera = vtk.vtkCamera()
#camera.SetViewUp(0, 0, -1)
#camera.SetPosition(0, 1, 0)
#camera.SetFocalPoint(0, 0, 0)
#camera.ComputeViewPlaneNormal()
#camera.SetRoll(180.0)
#camera.Azimuth(80.0)

renderer.SetBackground(0.85, 0.85, 0.85)
#renderer.SetActiveCamera(camera)
#renderer.ResetCamera()

readers = []
mappers = []
actors = []

for ref in refs:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(ref)
    readers.append(reader)

for instance in instances:
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInput(readers[shape[int(instance)]].GetOutput())
    mappers.append(mapper)
    actor = vtk.vtkActor()
    #    actor.GetProperty().SetOpacity(0.7)
    actor.GetProperty().SetColor(random_color())
    actor.SetMapper(mapper)
    actors.append(actor)
    renderer.AddActor(actor)
    transform = vtk.vtkTransform()
    actor.SetUserTransform(transform)
    transforms[instance] = transform

renderer.AddActor(gactor)
renderer.AddActor(cactor)
renderer.AddActor(sactor)




id_t0 = numpy.where(dpos_data[:, 0] == min(dpos_data[:, 0]))

pos_data = numpy.concatenate((spos_data, dpos_data))

set_positionv(pos_data[id_t0, 1], pos_data[id_t0, 2], pos_data[id_t0, 3],
              pos_data[id_t0, 4], pos_data[id_t0, 5], pos_data[id_t0, 6],
              pos_data[id_t0, 7], pos_data[id_t0, 8])

renderer_window.AddRenderer(renderer)
interactor_renderer.SetRenderWindow(renderer_window)
interactor_renderer.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
renderer_window.Render()

# http://www.itk.org/Wiki/VTK/Depth_Peeling ?

# #Use a render window with alpha bits (as initial value is 0 (false) ): 
# renderer_window.SetAlphaBitPlanes(1)

# # Force to not pick a framebuffer with a multisample buffer ( as initial value is 8): 
# renderer_window.SetMultiSamples(0)

# # Choose to use depth peeling (if supported) (initial value is 0 (false) ) 
# renderer.SetUseDepthPeeling(1)

# # Set depth peeling parameters. 
# renderer.SetMaximumNumberOfPeels(100)

# # Set the occlusion ratio (initial value is 0.0, exact image)
# renderer.SetOcclusionRatio(0.1)


class InputObserver():

    def __init__(self, times, slider_repres):
        self._stimes = set(times)
        self._opacity = 1.0
        self._time_step = (max(self._stimes) - min(self._stimes)) \
                           / len(self._stimes)
        self._time = min(times)
        self._slider_repres = slider_repres

    def update(self):
        index = bisect.bisect(times, self._time)
        cf_prov._time = times[index]
        cf_prov.method()
        contact_pos.SetInput(cf_prov._output)
        contact_pos.Update()
        #        contact_pos_force.Update()
        # arrow_glyph.Update()
        #gmapper.Update()

        id_t = numpy.where(pos_data[:, 0] == times[index])
        set_positionv(pos_data[id_t, 1], pos_data[id_t, 2], pos_data[id_t, 3],
                      pos_data[id_t, 4],
                      pos_data[id_t, 5], pos_data[id_t, 6], pos_data[id_t, 7],
                      pos_data[id_t, 8])

        self._slider_repres.SetValue(self._time)
        renderer_window.Render()

    def set_opacity(self):
        for actor in actors:
            actor.GetProperty().SetOpacity(self._opacity)

    def key(self, obj, event):
        key = obj.GetKeySym()
        print 'key', key

        if key == 'Up':
                self._time_step = self._time_step * 2.

        if key == 'Down':
                self._time_step = self._time_step / 2.

        if key == 'Left':
                self._time -= self._time_step

        if key == 'Right':
                self._time += self._time_step

        if key == 't':
                self._opacity -= 1
                self.set_opacity()

        if key == 'T':
                self._opacity += 1
                self.set_opacity()


        self.update()


    def time(self, obj, event):
        slider_repres = obj.GetRepresentation()
        self._time = slider_repres.GetValue()
        self.update()

slider_repres = vtk.vtkSliderRepresentation2D() 
min_time = times[0] 
max_time = times[len(times)-1]
slider_repres.SetMinimumValue(min_time)
slider_repres.SetMaximumValue(max_time)
slider_repres.SetValue(min_time)
slider_repres.SetTitleText("time")
slider_repres.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_repres.GetPoint1Coordinate().SetValue(0.4, 0.9)
slider_repres.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_repres.GetPoint2Coordinate().SetValue(0.9, 0.9)

slider_repres.SetSliderLength(0.02)
slider_repres.SetSliderWidth(0.03)
slider_repres.SetEndCapLength(0.01)
slider_repres.SetEndCapWidth(0.03)
slider_repres.SetTubeWidth(0.005)
slider_repres.SetLabelFormat("%3.4lf")
slider_repres.SetTitleHeight(0.02)
slider_repres.SetLabelHeight(0.02)

slider_widget = vtk.vtkSliderWidget()
slider_widget.SetInteractor(interactor_renderer)
slider_widget.SetRepresentation(slider_repres)
slider_widget.KeyPressActivationOff()
slider_widget.SetAnimationModeToAnimate()
slider_widget.SetEnabled(True)

input_observer = InputObserver(times, slider_repres)


slider_widget.AddObserver("InteractionEvent", input_observer.time)

interactor_renderer.AddObserver('KeyPressEvent', input_observer.key)


# Create a vtkLight, and set the light parameters.
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(0, 50, 50)

renderer.AddLight(light)

interactor_renderer.Initialize()
interactor_renderer.Start()
