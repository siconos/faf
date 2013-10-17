import vtk
from math import cos, sin, pi
class MyView():
    def __init__(self, renderer):
        self._renderer = renderer
        axes2 = vtk.vtkCubeAxesActor()
        axes2.SetBounds(
            -50, 50, -50, 50, 0.6, 0.6)
        axes2.SetDrawXGridlines(True)
        axes2.SetDrawYGridlines(True)
        #axes2.SetDrawZGridlines(True)
        camera = vtk.vtkCamera()
        camera.SetPosition(0,0,100)
        axes2.SetCamera(camera)
        axes2.SetXLabelFormat("%6.1f")
        axes2.SetYLabelFormat("%6.1f")
        axes2.SetZLabelFormat("%6.1f")
        try:
            axes2.GetXAxesGridlinesProperty().SetColor(0.1,0.1,0.1)
            axes2.GetYAxesGridlinesProperty().SetColor(0.1,0.1,0.1)
            axes2.GetLabelTextProperty().SetColor(0.1,0.1,0.1)
        except:
            pass
        #axes2.SetFlyModeToClosestTriad()
        axes2.SetScreenSize(20.0)
        self._renderer.AddViewProp(axes2)

    def action(self, observer):
        image_maker = vtk.vtkWindowToImageFilter()
        image_maker.SetInput(observer._renderer_window)
        image_maker.DebugOn()

        recorder = vtk.vtkOggTheoraWriter()
        recorder.SetQuality(2)  # best
        recorder.DebugOn()

        recorder.SetFileName('out.avi')
        recorder.SetInputConnection(image_maker.GetOutputPort())

        recorder.Start()
        lt = len(observer._times)
        observer._time = 0.
        camera = observer._renderer.GetActiveCamera()
        max_time = max(observer._times)
        i = 0
        steps = 100.0
        # camera position: (-2.5511092948434606, -61.51574013040903, 79.7316663804111)
        # camera focal point (1.157025654521458, 6.282142233587582, -0.0182814735287522)
        while observer._time < max_time:
            print 'camera position : ', camera.GetPosition()
            camera.SetPosition(100.*cos(2.*i*pi/lt), 100.*sin(2.*i*pi/lt), 100.*i/lt)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetRoll(-90)
            observer.update()
            image_maker.Modified()
            recorder.Write()
            observer._time += 10*observer._time_step
            i += 1
            
            #print 'step 2'
            #fx, fy, fz = camera.GetFocalPoint()
            #for i in range(0,int(steps)):
            #camera.SetFocalPoint(fx*(steps-i-1)/steps - 20*(i+1)/steps, fy*(steps-i-1)/steps, fz*(steps-i-1)/steps + (i+1)/steps)
            #observer.update()
            #image_maker.Modified()
            #recorder.Write()

        print 'step 2'
        fx, fy, fz = camera.GetFocalPoint()
        px, py, pz = camera.GetPosition()
        for i in range(0,int(steps)):
            camera.SetPosition(px / (i+1) + 40 * i/steps, py / (i+1), 10 * i/steps + pz / (i+1)) 
            camera.SetFocalPoint(fx*(steps-i-1)/steps - 20*(i+1)/steps, fy*(steps-i-1)/steps, fz*(steps-i-1)/steps + (i+1)/steps)
            print 'camera position : ', camera.GetPosition()
            camera.SetRoll(-90)
            observer.update()
            image_maker.Modified()
            recorder.Write()

        recorder.End()
