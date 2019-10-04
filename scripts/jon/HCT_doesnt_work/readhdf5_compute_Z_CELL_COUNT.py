import sys
sys.path.append("L:/hemocell/scripts/jon/")  # 'default'    
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import pickle

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 4          # Number of cores you want to read in with

#%%


CASE = "AR2"

# Where does the data sit?
ROOT = "F:/"
DATAPATH = ROOT + CASE + "/output/hdf5/"
RESULTS_PATH = ROOT + CASE + "/output/HCT/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40             # How many procs did you run on?
read_procs = 4          # Number of cores you want to read in with

# spatial boundaries
xmin, xmax, ymin, ymax, zmin, zmax = -1.5, 410.5, -1.5, 251.5, -1.5, 126.5

# function for empty sideview
def emptySideview(xmin, xmax, ymin, ymax):
    x_range = int(round(xmax - xmin)) + 1
    y_range = int(round(ymax - ymin)) + 1
    return np.zeros((x_range, y_range))    

def get_z_counts(fluid_position, fluid_boundary, xmin, xmax, ymin, ymax, zmin, zmax):
    
    z_counts = emptySideview(xmin, xmax, ymin, ymax)
    
    for pos_triplet, bdry in zip(fluid_position, fluid_boundary):
        x,y = pos_triplet[[0,1]]
        i = int(round(x - xmin))
        j = int(round(y - ymin))
        if int(bdry) == 0:
            z_counts[i,j] += 1
    
    return z_counts

#%%
# Read t=0 Fluid box, only for z-dimension cell counts that are IN geometry
fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=[0], datapath=DATAPATH, fluidname='Fluid')

#%%
fluid_position = np.array(fluid[0].position)
fluid_boundary = np.array(fluid[0].boundary)    

z_cell_counts = get_z_counts(fluid_position, fluid_boundary, xmin, xmax, ymin, ymax, zmin, zmax)

#%%

pickle.dump(z_cell_counts, open(RESULTS_PATH + CASE + "_z_cell_count.pickle", "wb"), protocol = 2)
