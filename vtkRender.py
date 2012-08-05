# An example from scipy cookbook demonstrating the use of numpy arrys in vtk
# Borrowed from Kitware (http://www.vtk.org/Wiki/VTK/Examples/Python/vtkWithNumpy)
 
import vtk
import cPickle as pkl
from numpy import *

def read_data( fname ):
    """
    Load the 3D array and try to make it a little smaller in memory. 
    """
    with open( fname ) as fh:
        data = pkl.load( fh )
        #data = load( fname )
    data = asarray( data, dtype=uint32 )
    return data


def render( data, height ):
    """
    Assume 3D data array with integer-valued input.
    """
    # For VTK to be able to use the data, it must be stored as a
    # VTK-image. This can be done by the vtkImageImport-class which
    # imports raw data and stores it.
    dataImporter = vtk.vtkImageImport()
    # The previously created array is converted to a string of chars and imported.
    data_string = data.tostring()
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToUnsignedChar()
    # Because the data that is imported only contains an intensity
    # value (it isnt RGB-coded or someting similar), the importer must
    # be told this is the case.
    dataImporter.SetNumberOfScalarComponents(1)
    # The following two functions describe how the data is stored and
    # the dimensions of the array it is stored in. For this simple
    # case, all axes are of length 75 and begins with the first
    # element. For other data, this is probably not the case.  I have
    # to admit however, that I honestly dont know the difference
    # between SetDataExtent() and SetWholeExtent() although VTK
    # complains if not both are used.
    nz, ny, nx = data.shape
    dataImporter.SetDataExtent(0, nx-1, 0, ny-1, 0, nz-1)
    dataImporter.SetWholeExtent(0, nx-1, 0, ny-1, 0, nz-1)

    # The following class is used to store transparencyv-values for
    # later retrival. In our case, we want the value 0 to be completly
    # opaque whereas the three different cubes are given different
    # transperancy-values to show how it works.
    alphaChannelFunc = vtk.vtkPiecewiseFunction()
    alphaChannelFunc.AddPoint(0, 0.0)
    alphaChannelFunc.AddPoint(255, 0.0)
    # alphaChannelFunc.AddPoint(dataAvg, 0.1)
    # alphaChannelFunc.AddPoint(dataMax, 0.2)

    # This class stores color data and can create color tables from a
    # the intensity values.
    lut = vtk.vtkLookupTable()
    lut.Build()
    lutNum = data.max()
    lut.SetNumberOfTableValues(lutNum)
    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToDiverging()
    ctf.AddRGBPoint(0.0, 0, 0, 1.0)
    ctf.AddRGBPoint(data.max(), 1.0, 0, 0 )
    # Conversion to RGB tuples based on intensity values -- coarsen
    # the number of height values to only 256
    for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
        cc = ctf.GetColor( ss )
        lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
   
    # The preavius two classes stored properties. Because we want to
    # apply these properties to the volume we want to render, we have
    # to store them in a class that stores volume properties.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor( ctf ) #colorFunc)
    volumeProperty.SetScalarOpacity(alphaChannelFunc)

    # This class describes how the volume is rendered (through ray tracing).
    compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
    # We can finally create our volume. We also have to specify the
    # data for it, as well as how the data will be rendered.
    volumeMapper = vtk.vtkVolumeRayCastMapper()
    volumeMapper.SetVolumeRayCastFunction(compositeFunction)
    volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

    # The class vtkVolume is used to pair the preaviusly declared
    # volume as well as the properties to be used when rendering that
    # volume.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    #create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
    # plane=vtk.vtkPlane()
    # plane.SetOrigin(0,100,100)
    # plane.SetNormal(1,0,0)
 
    # #create cutter
    # cutter=vtk.vtkCutter()
    # cutter.SetCutFunction(plane)
    # cutter.SetInputConnection(dataImporter.GetOutputPort())
    # cutter.Update()
    # cutterMapper=vtk.vtkPolyDataMapper()
    # cutterMapper.SetInputConnection( cutter.GetOutputPort())

    # #create plane actor
    # planeActor=vtk.vtkActor()
    # planeActor.GetProperty().SetColor(1.0,1,0)
    # planeActor.GetProperty().SetLineWidth(2)
    # planeActor.SetMapper(cutterMapper)

    #create cube actor
    # cubeActor=vtk.vtkActor()
    # cubeActor.GetProperty().SetColor(0.5,1,0.5)
    # cubeActor.GetProperty().SetOpacity(0.5)
    # cubeActor.SetMapper(cubeMapper)

    # With almost everything else ready, its time to initialize the
    # renderer and window, as well as creating a method for exiting
    # the application
    renderer = vtk.vtkRenderer()
    #renderer.AddActor( planeActor )
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
 
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)

    # We add the volume to the renderer ...
    renderer.AddVolume(volume)
    # ... set background color to white ...
    renderer.SetBackground(0,0,0)
    # ... and set window size.
    renderWin.SetSize(600, 600)

    # add a scalar color bar
    sb = vtk.vtkScalarBarActor() 
    sb.SetTitle("Elevation") 
    # If the orientation is vertical there is a problem. 
    sb.SetOrientationToHorizontal() 
    # Vertical is OK. 
    # sb.SetOrientationToVertical() 
    sb.SetWidth(0.6) 
    sb.SetHeight(0.17) 
    sb.SetPosition(0.1, 0.05) 
    sb.SetLookupTable(ctf) 

    sbw = vtk.vtkScalarBarWidget() 
    sbw.SetInteractor(renderInteractor) 
    sbw.SetScalarBarActor(sb) 
    sbw.On() 

    # A simple function to be called when the user decides to quit the application.
    def exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)

    # Tell the application to use the function as an exit check.
    renderWin.AddObserver("AbortCheckEvent", exitCheck)

    renderInteractor.Initialize()
    # Because nothing will be rendered without any input, we order the
    # first render manually before control is handed over to the
    # main-loop.
    renderWin.Render()
    renderInteractor.Start()

    # write a PNG image to disk
    writer = vtk.vtkPNGWriter()
    writer.SetFileName("rbc_stackVTK_height_"+str(height)+".png")
    writer.SetInput(dataImporter.GetOutput())
    writer.Write()


if __name__ == "__main__":

    #fname = "stacktest_height_1500.npy"
    fname = './data/new_110125_height_'
    #fname = "stack.pkl"
    # for i in range( 1000, 5000, 1000 ):
    ht = 2000
    datafile = fname + str( ht )
    arr = read_data( datafile )
    render( arr, height=ht )
        
        
