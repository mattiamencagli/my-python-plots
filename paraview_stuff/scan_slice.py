# state file generated using paraview version 5.11.2
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [3463, 1718]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [383.5, 255.5, 255.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2065.994861610943, 255.5, 255.5]
renderView1.CameraFocalPoint = [383.5, 255.5, 255.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 526.9086733011709
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = 'Z Axis'
renderView1.AxesGrid.ZTitle = 'X Axis'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(3463, 1718)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
density_t300000xmf = XDMFReader(registrationName='density_t.300000.xmf', FileNames=['/home/matti/DATA/Mayonese/staggered_run_1.5Z/density_t.300000.xmf'])
density_t300000xmf.CellArrayStatus = ['rho1']
density_t300000xmf.GridStatus = ['Grid_2']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=density_t300000xmf)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [766.0, 255.5, 255.5]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [383.5, 255.5, 255.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'rho1'
rho1TF2D = GetTransferFunction2D('rho1')
rho1TF2D.ScalarRangeInitialized = 1
rho1TF2D.Range = [0.18, 1.18, 0.0, 1.0]

# get color transfer function/color map for 'rho1'
rho1LUT = GetColorTransferFunction('rho1')
rho1LUT.AutomaticRescaleRangeMode = 'Never'
rho1LUT.TransferFunction2D = rho1TF2D
rho1LUT.RGBPoints = [0.18, 0.231373, 0.298039, 0.752941, 0.6799999999999999, 0.865003, 0.865003, 0.865003, 1.18, 0.705882, 0.0156863, 0.14902]
rho1LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'rho1']
slice1Display.LookupTable = rho1LUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 51.1
slice1Display.SelectScaleArray = 'rho1'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'rho1'
slice1Display.GaussianRadius = 2.555
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.SelectInputVectors = [None, '']
slice1Display.WriteLog = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for rho1LUT in view renderView1
rho1LUTColorBar = GetScalarBar(rho1LUT, renderView1)
rho1LUTColorBar.WindowLocation = 'Any Location'
rho1LUTColorBar.Position = [0.7170083742419868, 0.32188591385331783]
rho1LUTColorBar.Title = 'rho1'
rho1LUTColorBar.ComponentTitle = ''
rho1LUTColorBar.ScalarBarLength = 0.33000000000000035

# set color bar visibility
rho1LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'rho1'
rho1PWF = GetOpacityTransferFunction('rho1')
rho1PWF.Points = [0.18, 0.0, 0.5, 0.0, 1.18, 1.0, 0.5, 0.0]
rho1PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')