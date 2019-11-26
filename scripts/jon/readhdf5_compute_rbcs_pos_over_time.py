import sys
sys.path.append("L:/hemocell/scripts/jon/")  # 'default'    
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import csv
import pickle

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 4          # Number of cores you want to read in with
timestep =      10000
TIME_begin =    10000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.
DX = 5e-7;   DT = 1e-7

#%%

CASE = "AR3"

ROOT = "L:/no_backup/"

# Where does the data sit?
DATAPATH = ROOT + CASE + "/output/hdf5/"
RESULTS_PATH = ROOT + CASE + "/output/heatmap/"

Z_SLICES_GEOMETRY = list(np.arange(93.5, 110.5))

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
#read_procs = 16          # Number of cores you want to read in with
read_procs = 4          # Number of cores you want to read in with
timestep =      10000
TIME_begin =    10000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.

timesteps = list(range(TIME_begin, TIME_end + timestep, timestep))
N_t = len(timesteps)
#timesteps = [100000]
#N_t = 1

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


def getSideviewCounts(rbc, sideview, xmin, xmax, ymin, ymax, zmin, zmax):
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
    return sideview


#%%
boundaries_path = RESULTS_PATH + CASE + "_spatial_boundaries.pickle"
geo_slice_path = RESULTS_PATH + CASE + "_geometry_slice.pickle"

        
xmin, xmax, ymin, ymax, zmin, zmax = pickle.load(open(boundaries_path, "rb"), encoding = "latin1")
#xmin, xmax, ymin, ymax, zmin, zmax = -1.5, 410.5, -1.5, 251.5, -1.5, 410.5
geo_slice = pickle.load(open(geo_slice_path, "rb"))

#%%
# Read RBC data
x_range = int(round(xmax - xmin)) + 1
y_range = int(round(ymax - ymin)) + 1
sideview = np.zeros((x_range, y_range))
for t in timesteps:
    print(t, end=" ")
    _, rbc, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
        p=False, read_procs=read_procs, nprocs=nprocs, timesteps=[t], \
        datapath=DATAPATH, rbcname='RBC')
    sideview = getSideviewCounts(rbc[0], sideview, xmin, xmax, ymin, ymax, zmin, zmax)

#%%
csv_path = RESULTS_PATH + CASE + "_sideview.csv"
with open(csv_path, "w+") as my_csv:
    csvWriter = csv.writer(my_csv, delimiter=' ')
    csvWriter.writerow(["N_t", str(N_t)])
    csvWriter.writerows(sideview)