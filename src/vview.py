#!/usr/bin/env python
import sys
import shlex
import vtk
from vtk.util import numpy_support
from math import atan2, pi, cos, sin
import bisect
from numpy.linalg import norm
import numpy
import random
import h5py
import getopt


def usage():
    print '{0}: Usage'.format(sys.argv[0])
    print """
    {0} filename
    """

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'h',
                                   ['help'])
except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)

for o, a in opts:
    if o == '--help':
        usage()
        exit(0)

if len(args) == 1:
    filename = args[0]
elif len(args) > 1:
    usage()
    exit(1)


def random_color():
    r = random.uniform(0.1, 0.9)
    g = random.uniform(0.1, 0.9)
    b = random.uniform(0.1, 0.9)
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

    #ref_filename = args[0]
    #bind_filename = args[1]
    #pos_filename = args[2]
    #cf_filename = args[3]

ref_filename = 'ref.txt'
cfg_filename = 'shapes.cfg'
bind_filename = 'bindings.dat'
dpos_filename = 'dpos.dat'
spos_filename = 'spos.dat'
cf_filename = 'cf.dat'

refs = []
refs_attrs = []
with open(ref_filename, 'r') as ref_file:
    for line in ref_file:
        line_tokens = shlex.split(line)
        refs.append(line_tokens[0])
        refs_attrs.append([float(x) for x in line_tokens[1:]])

shape = dict()
with open(bind_filename, 'r') as bind_file:
    for line in bind_file:
        obj_id, obj_shape = shlex.split(line)
        shape[int(obj_id)] = int(obj_shape)

pos = dict()
instances = set()

import numpy

inFile = h5py.File('out.hdf5', 'r')

#spos_data = numpy.loadtxt(spos_filename, ndmin=2)
#dpos_data = numpy.loadtxt(dpos_filename, ndmin=2)
spos_data = inFile['data']['static']
dpos_data = inFile['data']['dynamic']

print spos_data.shape
print dpos_data.shape

#cf_data = numpy.loadtxt(cf_filename)
cf_data = inFile['data']['cf'][:].copy()

#def contact_point_reader():
#    global dos
#    dos.GetOutput().GetFieldData()
    

#contact_field = vtk.vtkPointData()
#c1 = vtk.vtkPolyData()

#contact_data = vtk.vtkDataObject()
contact_posa = vtk.vtkDataObjectToDataSetFilter()
contact_posb = vtk.vtkDataObjectToDataSetFilter()

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
            self._mu_coefs = set(self._data[:, 1])
        else:
            self._data = None
            self._mu_coefs = None

        if self._data is not None:
            self._time = min(self._data[:, 0])
        else:
            self._time = 0

        self._output = vtk.vtkPolyData()
        self._contact_field = vtk.vtkPointData()


    def method(self):
        global keeper

        #        ind0 = bisect.bisect_left(self._data[:, 0], self._time)
        self._contact_field = vtk.vtkPointData()
        self._output.SetFieldData(self._contact_field)

        if self._data is not None:
            id_f = numpy.where(abs(self._data[:, 0] - self._time) < 1e-15)[0]

            if len(id_f) == 0:
                return

            self.cpa_at_time = self._data[id_f, 2:5].copy()
            self.cpa = numpy_support.numpy_to_vtk(self.cpa_at_time)
            self.cpa.SetName('contactPositionsA')

            self.cpb_at_time = self._data[id_f, 5:8].copy()
            self.cpb = numpy_support.numpy_to_vtk(self.cpb_at_time)
            self.cpb.SetName('contactPositionsB')

            self.cn_at_time = - self._data[id_f, 8:11].copy()
            self.cn = numpy_support.numpy_to_vtk(self.cn_at_time)
            self.cn.SetName('contactNormals')

            self.cf_at_time = self._data[id_f, 11:14].copy()
            self.cf = numpy_support.numpy_to_vtk(self.cf_at_time)
            self.cf.SetName('contactForces')

            self._contact_field.AddArray(self.cpa)
            self._contact_field.AddArray(self.cpb)
            self._contact_field.AddArray(self.cn)
            self._contact_field.AddArray(self.cf)

        else:
            pass
        self._output.Update()

cf_prov = CFprov(cf_data)

contact_posa.SetDataSetTypeToPolyData()
contact_posa.SetPointComponent(0, "contactPositionsA", 0)
contact_posa.SetPointComponent(1, "contactPositionsA", 1)
contact_posa.SetPointComponent(2, "contactPositionsA", 2)

contact_posb.SetDataSetTypeToPolyData()
contact_posb.SetPointComponent(0, "contactPositionsB", 0)
contact_posb.SetPointComponent(1, "contactPositionsB", 1)
contact_posb.SetPointComponent(2, "contactPositionsB", 2)

contact_pos_force.SetInputConnection(contact_posa.GetOutputPort())
contact_pos_force.SetInputFieldToDataObjectField()
contact_pos_force.SetOutputAttributeDataToPointData()
contact_pos_force.SetVectorComponent(0, "contactForces", 0)
contact_pos_force.SetVectorComponent(1, "contactForces", 1)
contact_pos_force.SetVectorComponent(2, "contactForces", 2)

contact_pos_norm.SetInputConnection(contact_posa.GetOutputPort())
contact_pos_norm.SetInputFieldToDataObjectField()
contact_pos_norm.SetOutputAttributeDataToPointData()
contact_pos_norm.SetVectorComponent(0, "contactNormals", 0)
contact_pos_norm.SetVectorComponent(1, "contactNormals", 1)
contact_pos_norm.SetVectorComponent(2, "contactNormals", 2)

times = list(set(dpos_data[:, 0]))
times.sort()

ndyna = len(numpy.where(dpos_data[:, 0] == times[0]))

nstatic = len(numpy.where(spos_data[:, 0] == times[0]))

instances = set(dpos_data[:, 1]).union(set(spos_data[:, 1]))

cf_prov._time = min(times[:])

cf_prov.method()
contact_posa.SetInput(cf_prov._output)
contact_posa.Update()
contact_posb.SetInput(cf_prov._output)
contact_posb.Update()

contact_pos_force.Update()
contact_pos_norm.Update()

arrow = vtk.vtkArrowSource()
arrow.SetTipResolution(40)
arrow.SetShaftResolution(40)

cone = vtk.vtkConeSource()
cone.SetResolution(40)

if cf_prov._mu_coefs is not None:
    cone.SetRadius(min(cf_prov._mu_coefs))  # one coef!!

cylinder = vtk.vtkCylinderSource()
cylinder.SetRadius(.01)
cylinder.SetHeight(1)

sphere = vtk.vtkSphereSource()


# 1. scale = (scalar value of that particular data index);
# 2. denominator = Range[1] - Range[0];
# 3. scale = (scale < Range[0] ? Range[0] : (scale > Range[1] ? Range[1] : scale));
# 4. scale = (scale - Range[0]) / denominator;
# 5. scale *= scaleFactor;

arrow_glyph = vtk.vtkGlyph3D()
arrow_glyph.SetInputConnection(contact_pos_force.GetOutputPort())
arrow_glyph.SetSourceConnection(arrow.GetOutputPort())
arrow_glyph.ScalingOn()
arrow_glyph.SetScaleModeToScaleByVector()
arrow_glyph.SetRange(0, .01)
arrow_glyph.ClampingOn()
arrow_glyph.SetScaleFactor(5)
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
transform.Translate(-0.5, 0., 0.)

cone_glyph = vtk.vtkGlyph3D()
cone_glyph.SetSourceTransform(transform)

cone_glyph.SetInputConnection(contact_pos_norm.GetOutputPort())
cone_glyph.SetSourceConnection(cone.GetOutputPort())
#cone_glyph.ScalingOn()
#cone_glyph.SetScaleModeToScaleByVector()
#cone_glyph.SetRange(0, 100)
#cone_glyph.ClampingOn()
#cone_glyph.SetScaleFactor(1000000)
cone_glyph.SetVectorModeToUseVector()

cone_glyph.SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
cone_glyph.OrientOn()

ctransform = vtk.vtkTransform()
ctransform.Translate(-0.5, 0, 0)
ctransform.RotateWXYZ(90, 0, 0, 1)
cylinder_glyph = vtk.vtkGlyph3D()
cylinder_glyph.SetSourceTransform(ctransform)

cylinder_glyph.SetInputConnection(contact_pos_norm.GetOutputPort())
cylinder_glyph.SetSourceConnection(cylinder.GetOutputPort())
cylinder_glyph.ScalingOff()
cylinder_glyph.SetVectorModeToUseVector()

cylinder_glyph.SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
cylinder_glyph.OrientOn()

cmapper = vtk.vtkPolyDataMapper()
cmapper.SetInputConnection(cone_glyph.GetOutputPort())

clmapper = vtk.vtkPolyDataMapper()
clmapper.SetInputConnection(cylinder_glyph.GetOutputPort())


sphere_glypha = vtk.vtkGlyph3D()
sphere_glypha.SetInputConnection(contact_posa.GetOutputPort())
sphere_glypha.SetSourceConnection(sphere.GetOutputPort())
sphere_glypha.ScalingOn()
#sphere_glypha.SetScaleModeToScaleByVector()
#sphere_glypha.SetRange(-0.5, 2)
#sphere_glypha.ClampingOn()
sphere_glypha.SetScaleFactor(.1)
#sphere_glypha.SetVectorModeToUseVector()

sphere_glyphb = vtk.vtkGlyph3D()
sphere_glyphb.SetInputConnection(contact_posb.GetOutputPort())
sphere_glyphb.SetSourceConnection(sphere.GetOutputPort())
sphere_glyphb.ScalingOn()
#sphere_glyphb.SetScaleModeToScaleByVector()
#sphere_glyphb.SetRange(-0.5, 2)
#sphere_glyphb.ClampingOn()
sphere_glyphb.SetScaleFactor(.1)
#sphere_glyphb.SetVectorModeToUseVector()

#sphere_glyphb.SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
#sphere_glyph.OrientOn()

smappera = vtk.vtkPolyDataMapper()
smappera.SetInputConnection(sphere_glypha.GetOutputPort())
smapperb = vtk.vtkPolyDataMapper()
smapperb.SetInputConnection(sphere_glyphb.GetOutputPort())

#cmapper.SetScalarModeToUsePointFieldData()
#cmapper.SetColorModeToMapScalars()
#cmapper.ScalarVisibilityOn()
#cmapper.SelectColorArray('contactNormals')
#gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

cactor = vtk.vtkActor()
cactor.GetProperty().SetOpacity(0.4)
cactor.GetProperty().SetColor(0, 0, 1)
cactor.SetMapper(cmapper)

clactor = vtk.vtkActor()
#cactor.GetProperty().SetOpacity(0.4)
clactor.GetProperty().SetColor(1, 0, 0)
clactor.SetMapper(clmapper)

sactora = vtk.vtkActor()
sactora.GetProperty().SetColor(1, 0, 0)
sactora.SetMapper(smappera)

sactorb = vtk.vtkActor()
sactorb.GetProperty().SetColor(0, 1, 0)
sactorb.SetMapper(smapperb)

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
#camera.SetPosition(-221, 40, 204)
#camera.SetFocalPoint(0, 0, 0)
#camera.ComputeViewPlaneNormal()
#camera.SetRoll(180.0)
#camera.Azimuth(80.0)

#renderer.SetBackground(0.85, 0.85, 0.85)
#renderer.SetActiveCamera(camera)
#renderer.ResetCamera()

readers = []
mappers = []
actors = []

for ref, attrs in zip(refs, refs_attrs):
    if '.vtp' in ref:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(ref)
        readers.append(reader)
    else:
        if ref == 'Sphere':
            source = vtk.vtkSphereSource()
            source.SetRadius(attrs[0])

        elif ref == 'Cone':
            source = vtk.vtkConeSource()
            source.SetRadius(attrs[0])
            source.SetHeight(attrs[1])
            source.SetResolution(15)
            source.SetDirection(0, 1, 0) # needed

        elif ref == 'Cylinder':
            source = vtk.vtkCylinderSource()
            source.SetResolution(15)
            source.SetRadius(attrs[0])
            source.SetHeight(attrs[1])
            #           source.SetDirection(0,1,0)

        elif ref == 'Box':
            source = vtk.vtkCubeSource()
            source.SetXLength(attrs[0])
            source.SetYLength(attrs[1])
            source.SetZLength(attrs[2])

        elif ref == 'Capsule':
            sphere1 = vtk.vtkSphereSource()
            sphere1.SetRadius(attrs[0])
            sphere1.SetCenter(0, attrs[1] / 2, 0)
            sphere1.SetThetaResolution(15)
            sphere1.SetPhiResolution(15)
            sphere1.Update()

            sphere2 = vtk.vtkSphereSource()
            sphere2.SetRadius(attrs[0])
            sphere2.SetCenter(0, -attrs[1] / 2, 0)
            sphere2.SetThetaResolution(15)
            sphere2.SetPhiResolution(15)
            sphere2.Update()

            cylinder = vtk.vtkCylinderSource()
            cylinder.SetRadius(attrs[0])
            cylinder.SetHeight(attrs[1])
            cylinder.SetResolution(15)
            cylinder.Update()

            data = vtk.vtkMultiBlockDataSet()
            data.SetNumberOfBlocks(3)
            data.SetBlock(0, sphere1.GetOutput())
            data.SetBlock(1, sphere2.GetOutput())
            data.SetBlock(2, cylinder.GetOutput())
            source = vtk.vtkMultiBlockDataGroupFilter()
            source.AddInput(data)

        readers.append(source)


for instance in instances:
    mapper = vtk.vtkCompositePolyDataMapper()
    mapper.SetInputConnection(readers[shape[int(instance)]].GetOutputPort())
    mappers.append(mapper)
    actor = vtk.vtkActor()
    if int(instance) >= 0:
        actor.GetProperty().SetOpacity(0.7)
    actor.GetProperty().SetColor(random_color())
    actor.SetMapper(mapper)
    actors.append(actor)
    renderer.AddActor(actor)
    transform = vtk.vtkTransform()
    actor.SetUserTransform(transform)
    transforms[instance] = transform

renderer.AddActor(gactor)
renderer.AddActor(cactor)
renderer.AddActor(clactor)
renderer.AddActor(sactora)
renderer.AddActor(sactorb)

import imp
try:
    imp.load_source('myview', 'myview.py')
    import myview
    this_view = myview.MyView(renderer)
except IOError as e:
    pass

id_t0 = numpy.where(dpos_data[:, 0] == min(dpos_data[:, 0]))

pos_data = numpy.concatenate((spos_data, dpos_data))

set_positionv(pos_data[id_t0, 1], pos_data[id_t0, 2], pos_data[id_t0, 3],
              pos_data[id_t0, 4], pos_data[id_t0, 5], pos_data[id_t0, 6],
              pos_data[id_t0, 7], pos_data[id_t0, 8])

renderer_window.AddRenderer(renderer)
interactor_renderer.SetRenderWindow(renderer_window)
interactor_renderer.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

# http://www.itk.org/Wiki/VTK/Depth_Peeling ?

# #Use a render window with alpha bits (as initial value is 0 (false) ): 
renderer_window.SetAlphaBitPlanes(1)

# # Force to not pick a framebuffer with a multisample buffer ( as initial value is 8): 
renderer_window.SetMultiSamples(0)

# # Choose to use depth peeling (if supported) (initial value is 0 (false) ) 
renderer.SetUseDepthPeeling(1)

# # Set depth peeling parameters. 
renderer.SetMaximumNumberOfPeels(100)

# # Set the occlusion ratio (initial value is 0.0, exact image)
renderer.SetOcclusionRatio(0.1)


class InputObserver():

    def __init__(self, times, slider_repres):
        self._stimes = set(times)
        self._opacity = 1.0
        self._time_step = (max(self._stimes) - min(self._stimes)) \
                           / len(self._stimes)
        self._time = min(times)
        self._slider_repres = slider_repres
        self._current_id = vtk.vtkIdTypeArray()
        self._renderer = renderer
        self._renderer_window = renderer_window
        self._times = times

    def update(self):
        index = bisect.bisect_left(self._times, self._time)
        index = max(0, index)
        index = min(index, len(self._times)-1)        
        cf_prov._time = self._times[index]
        cf_prov.method()

        # contact_posa.SetInput(cf_prov._output)
        contact_posa.Update()

        # contact_posb.SetInput(cf_prov._output)
        contact_posb.Update()

        #contact_pos_force.Update()
        # arrow_glyph.Update()
        #gmapper.Update()

        id_t = numpy.where(pos_data[:, 0] == self._times[index])
        set_positionv(pos_data[id_t, 1], pos_data[id_t, 2], pos_data[id_t, 3],
                      pos_data[id_t, 4],
                      pos_data[id_t, 5], pos_data[id_t, 6], pos_data[id_t, 7],
                      pos_data[id_t, 8])

        self._slider_repres.SetValue(self._time)
        renderer_window.Render()

        self._current_id.SetNumberOfValues(1)
        self._current_id.SetValue(0,index)

        self._iter_plot.SetSelection(self._current_id)
        self._prec_plot.SetSelection(self._current_id)
        self._iter_plot_view.Update()
        self._prec_plot_view.Update()
        self._iter_plot_view.GetRenderer().GetRenderWindow().Render()
        self._prec_plot_view.GetRenderer().GetRenderWindow().Render()

    def set_opacity(self):
        for instance, actor in zip(instances, actors):
            if instance >= 0:
                actor.GetProperty().SetOpacity(self._opacity)

    def key(self, obj, event):
        key = obj.GetKeySym()
        print 'key', key

        if key == 'Up':
                self._time_step = self._time_step * 2.
                self._time += self._time_step

        if key == 'Down':
                self._time_step = self._time_step / 2.
                self._time -= self._time_step

        if key == 'Left':
                self._time -= self._time_step

        if key == 'Right':
                self._time += self._time_step

        if key == 't':
                self._opacity -= .1
                self.set_opacity()

        if key == 'T':
                self._opacity += .1
                self.set_opacity()

        if key == 'c':
                print 'camera position:',self._renderer.GetActiveCamera().GetPosition()
                print 'camera focal point', self._renderer.GetActiveCamera().GetFocalPoint()
                print 'camera clipping plane', self._renderer.GetActiveCamera().GetClippingRange()

                
        if key == 'C':
                this_view.action(self)
        self.update()

    def time(self, obj, event):
        slider_repres = obj.GetRepresentation()
        self._time = slider_repres.GetValue()
        self.update()

        # observer on 2D chart
    def iter_plot_observer(self, obj, event):

        if self._iter_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self._iter_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[self._iter_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

    def prec_plot_observer(self, obj, event):
        if self._prec_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self._prec_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[self._prec_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

slider_repres = vtk.vtkSliderRepresentation2D()
min_time = times[0]
max_time = times[len(times) - 1]
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
interactor_renderer.AddObserver('KeyPressEvent', input_observer.key)

# Create a vtkLight, and set the light parameters.
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(0, 0, 500)
#light.SetLightTypeToHeadlight()
renderer.AddLight(light)

hlight = vtk.vtkLight()
hlight.SetFocalPoint(0, 0, 0)
#hlight.SetPosition(0, 0, 500)
hlight.SetLightTypeToHeadlight()
renderer.AddLight(hlight)


#solv_data = numpy.loadtxt('solv.dat', ndmin=2)
solv_data = inFile['data']['solv']

#import numpy as np
#import matplotlib.pyplot as plt

#plt.plot(solv_data[:,0], solv_data[:,1])
#plt.show()


# Warning! numpy support offer a view on numpy array
# the numpy array must not be garbage collected!
nxtime = solv_data[:, 0].copy()
nxiters = solv_data[:, 1].copy()
nprecs = solv_data[:, 2].copy()
xtime = numpy_support.numpy_to_vtk(nxtime)
xiters = numpy_support.numpy_to_vtk(nxiters)
xprecs = numpy_support.numpy_to_vtk(nprecs)

xtime.SetName('time')
xiters.SetName('iterations')
xprecs.SetName('precisions')

table = vtk.vtkTable()
table.AddColumn(xtime)
table.AddColumn(xiters)
table.AddColumn(xprecs)
#table.Dump()

tview_iter = vtk.vtkContextView()
tview_prec = vtk.vtkContextView()

chart_iter = vtk.vtkChartXY()
chart_prec = vtk.vtkChartXY()
tview_iter.GetScene().AddItem(chart_iter)
tview_prec.GetScene().AddItem(chart_prec)
iter_plot = chart_iter.AddPlot(vtk.vtkChart.LINE)
iter_plot.SetLabel('Solver iterations')
iter_plot.GetXAxis().SetTitle('time')
iter_plot.GetYAxis().SetTitle('iterations')

prec_plot = chart_prec.AddPlot(vtk.vtkChart.LINE)
prec_plot.SetLabel('Solver precisions')
prec_plot.GetXAxis().SetTitle('time')
prec_plot.GetYAxis().SetTitle('precisions')

iter_plot.SetInput(table, 'time', 'iterations')
prec_plot.SetInput(table, 'time', 'precisions')
iter_plot.SetWidth(5.0)
prec_plot.SetWidth(5.0)
iter_plot.SetColor(0, 255, 0, 255)
prec_plot.SetColor(0, 255, 0, 255)

input_observer._iter_plot = iter_plot
input_observer._prec_plot = prec_plot
input_observer._iter_plot_view = tview_iter
input_observer._prec_plot_view = tview_prec

tview_iter.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                       input_observer.iter_plot_observer)

tview_prec.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                       input_observer.prec_plot_observer)


#screen_size = renderer_window.GetScreenSize()
renderer_window.SetSize(600, 600)
tview_iter.GetRenderer().GetRenderWindow().SetSize(600, 200)
tview_prec.GetRenderer().GetRenderWindow().SetSize(600, 200)


tview_iter.GetInteractor().Initialize()
#tview_iter.GetInteractor().Start()
tview_iter.GetRenderer().SetBackground(.9, .9, .9)
tview_iter.GetRenderer().Render()

tview_prec.GetInteractor().Initialize()
#tview_prec.GetInteractor().Start()
tview_prec.GetRenderer().SetBackground(.9, .9, .9)
tview_prec.GetRenderer().Render()

renderer_window.Render()
interactor_renderer.Initialize()
interactor_renderer.Start()
