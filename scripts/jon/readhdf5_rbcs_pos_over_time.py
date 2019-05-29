import sys
#sys.path.append("/home/jdbouter/hemocell/scripts/measure/")  # 'default'
sys.path.append("L:/hemocell/scripts/measure/")  # 'default'
import matplotlib
from matplotlib import pyplot as plt        
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import csv
import pickle

# Where does the data sit?
#datapath = "/home/jdbouter/no_backup/AR2_stiff/output/hdf5/"
datapath = "L:/no_backup/AR2/output/hdf5/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 4          # Number of cores you want to read in with
timestep =     100000
TIME_begin =   100000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.
DX = 5e-7;   DT = 1e-7

SAVE_FIGS = False

#%%

Z_SLICES_GEOMETRY = list(np.arange(93.5, 110.5))

#%%

timesteps = list(range(TIME_begin, TIME_end + timestep, timestep))
N_t = len(timesteps)

#%%

# Read t=0 Fluid box, only for the geometry boundaries
fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=[timesteps[0]], \
    datapath=datapath, fluidname='Fluid')
fluid=fluid[0]

#%%

# Read RBC data
#_, RBC, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
#    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=timesteps, \
#    datapath=datapath, rbcname='RBC')

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
def getGeometrySlice(fluid_position, fluid_boundary, xmin, xmax, ymin, ymax, zmin, zmax, Z_SLICES_GEOMETRY):
    x_range = int(round(xmax - xmin)) + 1
    y_range = int(round(ymax - ymin)) + 1    
    geo_slice = np.ones((x_range, y_range), dtype=int)
    EPS = 1e-5  # crude compare floats threshold
    for pos_triplet, bdry in zip(fluid_position, fluid_boundary):
        # Check if the z coordinate is in one of the specified slices
        z_coordinate_in_slice = False
        z = pos_triplet[2]
        for SLICE in Z_SLICES_GEOMETRY:
            if abs(z - SLICE) < EPS:
                z_coordinate_in_slice = True
        # Assign 1 in geo_slice if it's boundary (i.e. outside of geometry fluid field)
        if z_coordinate_in_slice and int(bdry) == 0:
            x,y = pos_triplet[0:2]
            i = int(round(x - xmin))
            j = int(round(y - ymin))
            geo_slice[i,j] = 0
    return geo_slice

def computeTimeAveragedSideview(RBC, timesteps, xmin, xmax, ymin, ymax, zmin, zmax):
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
    
        with open("L:/no_backup/AR2_stiff/output/sideview.csv", "w+") as my_csv:
            csvWriter = csv.writer(my_csv, delimiter=' ')
            try:
                csvWriter.writerow(["time n/N", str(timesteps[t_idx]), str(t_idx) + "/" + str(N_t)])
            except:
                pass
            csvWriter.writerows(sideview)
    
    sideview /= N_t
    return sideview

def plotSideview(sideview, geo_slice, N_t):
    # First create two color maps. One for the geometry contour and one with 
    # transparency for the overlaid RBC heat map
    from matplotlib.colors import LinearSegmentedColormap

    # Color map for geometry
    N = 256
    geometry_color = (200.0, 200.0, 200.0)
    geometry_color = [c/N for c in geometry_color]
    geometry_color.append(1.0)  # no transparency
    white_c = [1.0, 1.0, 1.0, 1.0]
    
    geometry_color_array = np.array([geometry_color, white_c])
    map_object = LinearSegmentedColormap.from_list(name='geo_map', colors=geometry_color_array)
    plt.register_cmap(cmap=map_object)
    
    # Comment here
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], 1.0, N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], 0.0, N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], 0.0, N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)
    
    # Changed sideview so that it's a fraction of timesteps
    sideview /= N_t
    
    # Mirror both matrices
    geo_slice = np.array(geo_slice[::-1])
    sideview  = np.array(sideview[::-1])
    
    # Transform sideview data with log for log scaling
    TwoPlusLog = lambda v : 2 + np.log10(v) if v > 0.0 else 0.0
    sideview = np.array([[TwoPlusLog(v) for v in row] for row in sideview])
    
    # Plot
    FONTSIZE = 24
    
    matplotlib.rc('xtick', labelsize=FONTSIZE) 
    matplotlib.rc('ytick', labelsize=FONTSIZE) 
    fig, ax = plt.subplots(1, 1)
    ax.imshow(geo_slice.transpose(), cmap='geo_map')
    fig2 = ax.imshow(sideview.transpose(), cmap='RBC_colors')
    colorbar = plt.colorbar(fig2)
    colorbar.ax.set_ylabel('Avg # of RBCs present over time', rotation=270, labelpad=40, fontsize=FONTSIZE)
    colorbar.set_ticks([0.0, 0.5, 1.0, 1.5])
    colorbar.set_ticklabels(["0.01", r"10$^{-1.5}$", "0.1", r"10$^{-0.5}$"])
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
    ax.set_xlim([0, 400])
    ax.set_ylim([210, 0])
    plt.grid(False)
    plt.title("AR2 RBC location averaged over time", fontsize=FONTSIZE)
    plt.xlabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    plt.ylabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    if SAVE_FIGS:
        matplotlib.use('Agg')
        plt.savefig("RBCs_pos_over_time")
    else:
        matplotlib.use('Qt5Agg')
        plt.show()

#%%        
#xmin, xmax, ymin, ymax, zmin, zmax = getSpatialBoundaries(fluid_position)
#
#geo_slice = getGeometrySlice(fluid_position, np.array(fluid.boundary), xmin, \
#                             xmax, ymin, ymax, zmin, zmax, Z_SLICES_GEOMETRY)
#
#
#pickle.dump(geo_slice, open("L:/hemocell/scripts/jon/AR2_geometry_slice.pickle", "wb"), protocol = 2)
#pickle.dump([xmin, xmax, ymin, ymax, zmin, zmax], \
#            open("L:/hemocell/scripts/jon/AR2_spatial_boundaries.pickle", "wb"), protocol = 2)
xmin, xmax, ymin, ymax, zmin, zmax = \
    pickle.load(open("L:/hemocell/scripts/jon/AR2_spatial_boundaries.pickle", "rb"))
geo_slice = pickle.load(open("L:/hemocell/scripts/jon/AR2_geometry_slice.pickle", "rb"))

#%%
#sideview = computeTimeAveragedSideview(RBC, timesteps, xmin, xmax, ymin, ymax, zmin, zmax)
sideview = np.loadtxt(open("L:/no_backup/AR2/output/sideview.csv", "rb"), delimiter=" ", skiprows=1)

plotSideview(sideview, geo_slice, N_t)

#%%