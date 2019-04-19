import sys
sys.path.append("C:/Users/jdabo/MYSTUFF/STAGE_local/hemocell/scripts/measure/")  # 'default'
# Command-line arguments
args = sys.argv
for i, arg in enumerate(args):
    if arg == "-l":      # Readhdf5 library
        sys.path.append(args[i+1])
from matplotlib import pyplot as plt        
import HCELL_readhdf5_boundary as HCELL_readhdf5
import numpy as np
import numpy.ma as ma

# Where does the data sit?
#datapath = "L:/hemocell/examples/short_AR1/output_0/hdf5/"
#datapath = "C:/Users/jdabo/MYSTUFF/STAGE_local/hemocell/examples/with_aneurysm/output/hdf5/"
datapath = "C:/Users/jdabo/Desktop/long_run/AR2_R0.049/output/hdf5/" # WINDOWS

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 44           # How many procs did you run on?
read_procs = 4        # Number of cores you want to read in with
timestep =    100000  # What was the time step?
TIME_begin = 7000000   # When do you want to start reading the data
TIME_end =   8100000   # When do you want to end reading the data.
DX = 5e-7;   DT = 1e-7

SAVE_FIGS = False
for i, arg in enumerate(args):
    if arg == "-d":      # data
        datapath = args[i+1]
    if arg == "-tmin":   # time first step (inclusive)
        TIME_begin = int(args[i+1])
    if arg == "-tmax":   # time last step (exclusive I think)
        TIME_end = int(args[i+1])        
    if arg == "-tstep":  # time step
        timestep = int(args[i+1])
    if arg == "-n":      # Number of processors the sim was run on
        nprocs = int(args[i+1])
    if arg == "save_figs":
        SAVE_FIGS = True

#%%

Z_SLICE_GEOMETRY = 67.5

#%%

# Read t=0 Fluid box, only for the geometry boundaries
fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, begin=0, \
    end=0+timestep, timestep=timestep, datapath=datapath, fluidname='Fluid')
fluid=fluid[0]

#%%

# Read RBC data
_, RBC, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
    p=False, read_procs=read_procs, nprocs=nprocs, begin=TIME_begin, \
    end=TIME_end, timestep=timestep, datapath=datapath, rbcname='RBC')

#%%

fluid_position = np.array(fluid.position)

#%%

# Get the x,y,z boundaries of geometry bounding box
def getSpatialBoundaries(fluid_position):
    xmin, ymin, zmin =  np.infty,  np.infty,  np.infty
    xmax, ymax, zmax = -np.infty, -np.infty, -np.infty
    for pos_triplet in fluid_position:        
        # Update min/max boundaries
        if pos_triplet[0] < xmin:
            xmin = pos_triplet[0]
        if pos_triplet[0] > xmax:
            xmax = pos_triplet[0]
        if pos_triplet[1] < ymin:
            ymin = pos_triplet[1]
        if pos_triplet[1] > ymax:
            ymax = pos_triplet[1]
        if pos_triplet[2] < zmin:
            zmin = pos_triplet[2]
        if pos_triplet[2] > zmax:
            zmax = pos_triplet[2]
            
    return xmin, xmax, ymin, ymax, zmin, zmax

# Get a 2d-array of geometry boundary cells
def getGeometrySlice(fluid_position, fluid_boundary, xmin, xmax, ymin, ymax, zmin, zmax, Z_SLICE_GEOMETRY):
    x_range = int(round(xmax - xmin)) + 1
    y_range = int(round(ymax - ymin)) + 1    
    geo_slice = np.zeros((x_range, y_range), dtype=int)
    EPS = 1e-5  # crude compare floats threshold
    for pos_triplet, bdry in zip(fluid_position, fluid_boundary):
        if abs(pos_triplet[2] - Z_SLICE_GEOMETRY) < EPS and int(bdry) == 1:
            x,y = pos_triplet[0:2]
            i = int(round(x - xmin))
            j = int(round(y - ymin))
            geo_slice[i,j] = 1
    return geo_slice

def computeTimeAveragedSideview(RBC, xmin, xmax, ymin, ymax, zmin, zmax):
    x_range = int(round(xmax - xmin)) + 1
    y_range = int(round(ymax - ymin)) + 1
    
    sideview = np.zeros((x_range, y_range))
    N_t = len(RBC)  # number of time steps
    
    for t_idx in range(N_t):
        rbc = RBC[t_idx]
        # 2d array of lists, dimensions x/y
        unique_ids_per_dxdy = [[ [] for j in range(y_range)] for i in range(x_range)]  
        for pos_triplet, cid in zip(np.array(rbc.position), np.array(rbc.cid)):
            x,y = pos_triplet[0:2]
            # LSP positions can be outside normal xmin/xmax range; exclude
            if x > xmin and x < xmax and y > ymin and y < ymax:  
                i = int(round(x - xmin))
                j = int(round(y - ymin))
                if not cid in unique_ids_per_dxdy[i][j]:
                    unique_ids_per_dxdy[i][j].append(cid)
            
        sideview += np.array( [ [len(id_list) for id_list in list_of_id_lists] \
                              for list_of_id_lists in unique_ids_per_dxdy]   )
    sideview /= N_t
    return sideview

def plotSideview(sideview, geo_slice):
    # First create two color maps. One for the geometry contour and one with 
    # transparency for the overlaid RBC heat map
    from matplotlib.colors import LinearSegmentedColormap

    # Color map for geometry
    N = 256
    geometry_color = (125.0, 137.0, 209.0)
    geometry_color = [c/N for c in geometry_color]
    geometry_color.append(0.3)  # half transparent
    
    geometry_color_array = np.array([geometry_color, [1.0, 1.0, 1.0, 1.0]])
    map_object = LinearSegmentedColormap.from_list(name='geo_map', colors=geometry_color_array)
    plt.register_cmap(cmap=map_object)
    
    # Color map for RBCs, interpolating between zero-value = geometry color and 1.0 = red
    # WRITE EXPLANATION FOR THR
    THR = int(np.floor(N/2))
    c1 = np.ones(N)
    c1[:THR] = np.linspace(geometry_color[0], 1.0, THR)
    c2 = np.zeros(N)
    c2[:THR] = np.linspace(geometry_color[1], 0.0, THR)
    c3 = np.zeros(N)
    c3[:THR] = np.linspace(geometry_color[2], 0.0, THR)
    c4 = np.ones(N)
    c4[:THR] = np.linspace(0.0, 1.0, THR)
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)
    
    # Plot
    fig, ax = plt.subplots(1, 1)
    ax.imshow(geo_slice.transpose(), cmap='geo_map')
    ax.imshow(sideview.transpose(), cmap='RBC_colors')
    ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    spatial_conv_factor = 0.5
    xlabels = [str(tick * spatial_conv_factor) for tick in xticks]
    ylabels = [str(tick * spatial_conv_factor) for tick in yticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)    
    #ax.set_xlim([0, 800])
    #ax.set_ylim([0, 300])
    plt.grid(True)
    plt.title("RBC location averaged over time")
    plt.xlabel(r"Position in $\mu$m")
    plt.ylabel(r"Position in $\mu$m")
    if SAVE_FIGS:
        plt.savefig("RBCs_pos_over_time")
    else:
        plt.show()
        
# %%
xmin, xmax, ymin, ymax, zmin, zmax = getSpatialBoundaries(fluid)

#%%
#bleh = np.array(fluid.position)
geo_slice = getGeometrySlice(fluid_position, np.array(fluid.boundary), xmin, \
                             xmax, ymin, ymax, zmin, zmax, Z_SLICE_GEOMETRY)

#%%
sideview = computeTimeAveragedSideview(RBC, xmin, xmax, ymin, ymax, zmin, zmax)

#%%
#plotSideview(sideview, geo_slice)
plotSideview(sideview, geo_slice)

