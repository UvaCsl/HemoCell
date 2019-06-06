import sys
#sys.path.append("/home/jdbouter/hemocell/scripts/measure/")  # 'default'
sys.path.append("L:/hemocell/scripts/jon/")  # 'default'    
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import csv
import pickle

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 15          # Number of cores you want to read in with
timestep =     100000
TIME_begin =   100000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.
DX = 5e-7;   DT = 1e-7

#%%

CASE = "AR1"

# Where does the data sit?
#datapath = "/home/jdbouter/no_backup/" + CASE + "/output/hdf5/"
DATAPATH = "L:/no_backup/" + CASE + "/output/hdf5/"
RESULTS_PATH = "L:/no_backup/" + CASE + "/output/heatmap/"

Z_SLICES_GEOMETRY = list(np.arange(93.5, 110.5))

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 16          # Number of cores you want to read in with
read_procs = 4          # Number of cores you want to read in with
timestep =     100000
TIME_begin =   100000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.

#timesteps = list(range(TIME_begin, TIME_end + timestep, timestep))
#N_t = len(timesteps)

timesteps = [100000]
N_t = 1

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

# Get a 2d-array of geometry boundary cells (i.e. 1 for outside of geo and 0 for inside)
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
        # Assign 0 in geo_slice if this triplet is in the slice and it's not boundary
        if z_coordinate_in_slice and int(bdry) == 0:
            x,y = pos_triplet[0:2]
            i = int(round(x - xmin))
            j = int(round(y - ymin))
            geo_slice[i,j] = 0
    return geo_slice


def writeSideviewCounts(RBC, timesteps, xmin, xmax, ymin, ymax, zmin, zmax):
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

    csv_path = RESULTS_PATH + CASE + "_sideview.csv"
    with open(csv_path, "w+") as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=' ')
        csvWriter.writerow(["N_t", str(N_t)])
        csvWriter.writerows(sideview)

#%%
# Read t=0 Fluid box, only for the geometry boundaries
fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=[timesteps[0]], \
    datapath=DATAPATH, fluidname='Fluid')

fluid_position = np.array(fluid[0].position)
fluid_boundary = np.array(fluid[0].boundary)

#%%

# Find the spatial boundaries of the fluid data
xmin, xmax, ymin, ymax, zmin, zmax = getSpatialBoundaries(fluid_position)
pickle.dump([xmin, xmax, ymin, ymax, zmin, zmax], \
            open(RESULTS_PATH + CASE + "_spatial_boundaries.pickle", "wb"), protocol = 2)

# Compute a 2D boolean array that represents whether there was geometry there
# in the given range of Z slices
geo_slice = getGeometrySlice(fluid_position, fluid_boundary, xmin, xmax, \
                                  ymin, ymax, zmin, zmax, Z_SLICES_GEOMETRY)

pickle.dump(geo_slice, open(RESULTS_PATH + CASE + "_geometry_slice.pickle", "wb"), protocol = 2)

#%%

# Read RBC data
_, RBC, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=timesteps, \
    datapath=DATAPATH, rbcname='RBC')

#%%

# Write counts of RBC per x/y coordinate to file
writeSideviewCounts(RBC, timesteps, xmin, xmax, ymin, ymax, zmin, zmax)