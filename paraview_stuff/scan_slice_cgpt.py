# state file generated using paraview version 5.11.2
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# Enable offscreen rendering
paraview.options.batch = True
paraview.options.symmetric = True

# Setup views used in the visualization
materialLibrary1 = GetMaterialLibrary()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]  # Set a standard size for images
renderView1.StereoType = 'Crystal Eyes'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1  

renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = 'Z Axis'
renderView1.AxesGrid.ZTitle = 'X Axis'
renderView1.AxesGrid.YTitle = 'Y Axis'

# Increase CameraParallelScale for better visibility of axes
renderView1.CameraParallelScale = 550.0

SetActiveView(None)

# Setup view layouts
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1920, 1080)

SetActiveView(renderView1)

# Load the data
data = XDMFReader(
    registrationName='density_t.600000.xmf', 
    FileNames=['/home/matti/DATA/Mayonese/staggered_run_1.5Z/density_t.600000.xmf']
)
data.CellArrayStatus = ['rho1']
data.GridStatus = ['Grid_2']

# Create a slice object
slice1 = Slice(registrationName='Slice1', Input=data)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# Visualization settings
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'rho1']

# Color map setup
rho1LUT = GetColorTransferFunction('rho1')
rho1PWF = GetOpacityTransferFunction('rho1')

slice1Display.LookupTable = rho1LUT
slice1Display.OpacityTransferFunction = rho1PWF

# Set up color legend
rho1LUTColorBar = GetScalarBar(rho1LUT, renderView1)
rho1LUTColorBar.Visibility = 1
slice1Display.SetScalarBarVisibility(renderView1, True)

# Loop through slice positions and save images
output_directory = "/home/matti/DATA/Mayonese/staggered_run_1.5Z/images_slice_600000"
import os
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for i in range(1, 768):  # Loop from 1 to 767
    slice1.SliceType.Origin = [i, 255.5, 255.5]  # Update slice position
    
    # Adjust camera position to follow the slice
    renderView1.CameraPosition = [i + 1200, 255.5, 255.5]
    renderView1.CameraFocalPoint = [i, 255.5, 255.5]

    renderView1.AxesGrid.Visibility = 0  # Temporarily hide axes
    Render()
    renderView1.AxesGrid.Visibility = 1  # Show them again
    Render() 

    # Render and save the image
    filename = os.path.join(output_directory, f"slice_{i:04d}.png")
    SaveScreenshot(filename, renderView1, ImageResolution=[1920, 1080])
    print(f"Saved {filename}")

print("Image sequence saved successfully.")

# ffmpeg -framerate 15 -i slice_0%03d.png -c:v libx264 -r 15 -pix_fmt yuv420p output.mp4
