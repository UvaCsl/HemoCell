import sys
sys.path.append("C:/Users/jdabo/Dropbox/STAGE/code/jon/")  # 'default'    
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import csv
import pickle

#%%

# function for empty sideview
def emptySideview(xmin, xmax, ymin, ymax):
    x_range = int(round(xmax - xmin)) + 1
    y_range = int(round(ymax - ymin)) + 1
    return np.zeros((x_range, y_range))    

def getHCT(rbc, Z_CELL_COUNTS, xmin, xmax, ymin, ymax, zmin, zmax):
    
    # First, sideview is going to keep track of 2D map of LSP-counts on that x,y coordinate
    sideview = emptySideview(xmin, xmax, ymin, ymax)
    LSP_count = 0
    
    rbc_position = rbc.position
        
    # iterate through LSPs
    for pos_triplet in np.array(rbc_position):
        x,y,z = pos_triplet

        # LSP positions can be outside spatial boundaries: exclude
        if x > xmin and x < xmax and y > ymin and y < ymax and z > zmin and z < zmax:
            i = int(round(x - xmin))
            j = int(round(y - ymin))
            sideview[i,j] += 1
            LSP_count += 1
        
    # Divide by z-volume
    for i, row in enumerate(sideview):
        for j, value in enumerate(row):
            # Don't bother with empty space
            if value != 0.0:
                # This check if for possible LSPs slightly outside of geometry boundaries
                if Z_CELL_COUNTS[i,j] != 0:    
                    sideview[i,j] = value / float(Z_CELL_COUNTS[i,j])
                else:
                    LSP_count -= sideview[i,j]
                    sideview[i,j] = 0
                    
    # Count number of unique RBCs
    RBC_count = len(np.unique(rbc.cid))

    # ugly corrective step: This method overestimates the RBC volume, but I know 
    # what the RBC volume should be, so I'm just going to divide all the HCT values 
    # by overestimated_volume / true_volume
    volume_overestimation = LSP_count * LATTICE_VOLUME
    fraction_of_LSPS_to_be_plotted = LSP_count / len(rbc_position)
    volume_true = RBC_count * RBC_VOLUME * fraction_of_LSPS_to_be_plotted
    corrective_factor = volume_true / volume_overestimation
    
    # CORRECT
    for i, row in enumerate(sideview):
        for j, value in enumerate(row):
            value = sideview[i,j]
            sideview[i,j] = value * corrective_factor

    return sideview

#%%

#CASES = ["AR1", "AR2", "AR2_stiff", "AR3"]
CASES = ["AR2"]
CASE = "AR2"
# spatial boundaries
xmin, xmax, ymin, ymax, zmin, zmax = -1.5, 410.5, -1.5, 251.5, -1.5, 126.5
sideviews = [emptySideview(xmin, xmax, ymin, ymax) for i in range(4)]

for k, CASE in enumerate(CASES):
    
    # Where does the data sit?
    ROOT = "L:/no_backup/"
    DATAPATH = ROOT + CASE + "/output/hdf5/"
    RESULTS_PATH = ROOT + CASE + "/output/HCT/"
    
    #FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
    nprocs = 40             # How many procs did you run on?
    read_procs = 4          # Number of cores you want to read in with
    timestep =      10000
    TIME_begin =  1010000   # When do you want to start reading the data
    TIME_end =    1310000   # When do you want to end reading the data.
    
    timesteps = list(range(TIME_begin, TIME_end + timestep, timestep))
    N_t = len(timesteps)
    
    RBC_VOLUME = 90000.0    # in attoliters or 10^-21 m^3 
    LATTICE_VOLUME = 125.0      # in attoliters or 10^-21 m^3 
        
    
    # number of intra-geometry fluid cells in z direction per x/y coordinate
    Z_CELL_COUNTS_PATH = RESULTS_PATH + CASE + "_z_cell_count.pickle"
    Z_CELL_COUNTS = pickle.load(open(Z_CELL_COUNTS_PATH, "rb"), encoding='latin1')

    # Read RBC data
    sideview = emptySideview(xmin, xmax, ymin, ymax)

    for t in timesteps:
        print(t, end=" ")
        # Load RBCs for timestep t
        _, rbc, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
            p=False, read_procs=read_procs, nprocs=nprocs, timesteps=[t], \
            datapath=DATAPATH, rbcname='RBC')
        
        # Add HCT values to sideview. (average over time later)
        sideview += getHCT(rbc[0], Z_CELL_COUNTS, xmin, xmax, ymin, ymax, zmin, zmax)
        
        # Flush print stdout
        sys.stdout.flush()
    
    # Average over time
    sideviews[k] = np.array([[v / float(N_t) for v in row] for row in sideview])
    
    
    #%%
    RESULTS_PATH = "C:/Users/jdabo/Desktop/"
    
    for CASE, sideview in zip(CASES, sideviews):
        csv_path = RESULTS_PATH + CASE + "_sideview.csv"
        with open(csv_path, "w+") as my_csv:
            csvWriter = csv.writer(my_csv, delimiter=' ')
            csvWriter.writerow(["N_t", str(N_t)])
            csvWriter.writerows(Z_CELL_COUNTS)
    