# trace generated using paraview version 5.8.0-RC2
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`
import os
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


#INPUT =============================
t_0 = 0
t_f = 0
timestep = 2000

datalocation = '/mnt/extra_storage/ben/hemocell/examples/oneCellShear/tmp/'
CELLname = "RBC"
#INPUT =============================


def convert_CELL_time(time,CELLname,datalocation):

    #first create new x3d path to save all created x3d files.
    x3dpath = os.path.join(datalocation,"x3d")
    if not os.path.isdir(x3dpath):
        os.mkdir(x3dpath)


    print("Converting "+CELLname+'.'+str(time).zfill(12)+'.xmf to x3d')

    # create a new 'XDMF Reader'
    CELL = XDMFReader(FileNames=[datalocation+CELLname+'.'+str(time).zfill(12)+'.xmf'])
    CELL.PointArrayStatus = ['Area force', 'Bending force', 'Link force', 'Position', 'Total force', 'Viscous force', 'Volume force']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [903, 796]
    # get layout
    layout1 = GetLayout()
    # show data in view
    CELLDisplay = Show(CELL, renderView1, 'UnstructuredGridRepresentation')
    # reset view to fit data
    renderView1.ResetCamera()
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    # update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(Input=CELL)
    # create a new 'Smooth'
    smooth1 = Smooth(Input=extractSurface1)
    # Properties modified on smooth1
    smooth1.NumberofIterations = 300
    # show data in view
    smooth1Display = Show(smooth1, renderView1, 'GeometryRepresentation')
    # hide data in view
    Hide(extractSurface1, renderView1)
    # update the view to ensure updated data information
    renderView1.Update()
    # export view
    ExportView(datalocation+'x3d/'+CELLname+'.'+str(time).zfill(12)+'.x3d', view=renderView1)
    
    Delete(layout1)
    Delete(renderView1)
    Delete(smooth1)
    Delete(smooth1Display)
    Delete(extractSurface1)
    Delete(CELL)
    Delete(CELLDisplay)


for t in range(t_0,t_f+timestep,timestep):
    convert_CELL_time(time = t,CELLname=CELLname,datalocation=datalocation)
    
    #If you have another celltype in the simultion, you can just run the same function, just change the cellname
    #convert_RBC_time(time = t,RBCname="RBC_stiff",datalocation=datalocation)
    
