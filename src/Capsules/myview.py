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
        # camera.SetPosition(655, -148, 1073)
        # camera.SetFocalPoint(0,0,0)
        # camera.SetClippingRange(0.1,1000)
        # for i in range(0, int(steps)+1):
        #     camera.SetPosition(100*i/steps + 655*(steps-i)/steps,
        #                        0*i/steps -148*(steps-i)/steps,
        #                        70*i/steps + 1073*(steps-i)/steps)
        #     camera.SetClippingRange(1,2000)
        #     observer.update()
        #     image_maker.Modified()
        #     recorder.Write()

        i=0
        while observer._time < max_time:
            #            print 'camera position : ', camera.GetPosition()
            #camera.SetPosition(100.*cos(2.*i*pi/steps), 100.*sin(2.*i*pi/steps), 70)
            #camera.SetFocalPoint(0, 0, 0)
            #camera.SetRoll(-90)
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

        recorder.End()
