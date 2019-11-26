import sys
sys.path.append("L:/hemocell/scripts/jon/")  # 'default'    
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import csv
import pickle

#%%

CASE = "AR2"

# Where does the data sit?
DATAPATH = "L:/no_backup/" + CASE + "/output/hdf5/"
RESULTS_PATH = "L:/no_backup/" + CASE + "/output/heatmap_PLT/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 4          # Number of cores you want to read in with
timestep =      10000
TIME_begin =    10000   # When do you want to start reading the data
TIME_end =   10000000   # When do you want to end reading the data.

timesteps = list(range(TIME_begin, TIME_end + timestep, timestep))

#%%

def getSideviewCounts(plt, sideview, xmin, xmax, ymin, ymax, zmin, zmax):
    # 2d array of lists, dimensions x/y
    unique_ids_per_dxdy = [[ [] for j in range(y_range)] for i in range(x_range)]  
    for pos_triplet, cid in zip(np.array(plt.position), np.array(plt.cid)):
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

xmin, xmax, ymin, ymax, zmin, zmax = pickle.load(open(boundaries_path, "rb"), encoding='latin1')
#xmin, xmax, ymin, ymax, zmin, zmax = -1.5, 410.5, -1.5, 251.5, -1.5, 410.5
geo_slice = pickle.load(open(geo_slice_path, "rb"))

#%%
# Read RBC data
x_range = int(round(xmax - xmin)) + 1
y_range = int(round(ymax - ymin)) + 1
sideview = np.zeros((x_range, y_range))
for t in timesteps:
    print(t, end=" ")
    _, _, plt, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=False, \
        p=True, read_procs=read_procs, nprocs=nprocs, timesteps=[t], \
        datapath=DATAPATH, pltname='PLT')
    sideview = getSideviewCounts(plt[0], sideview, xmin, xmax, ymin, ymax, zmin, zmax)

sideview /= len(timesteps)

#%%
csv_path = RESULTS_PATH + CASE + "_PLT_sideview.csv"
with open(csv_path, "w+") as my_csv:
    csvWriter = csv.writer(my_csv, delimiter=' ')
    csvWriter.writerow(["N_t", str(len(timesteps))])
    csvWriter.writerows(sideview)