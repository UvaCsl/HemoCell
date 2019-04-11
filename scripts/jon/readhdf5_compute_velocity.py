from matplotlib import pyplot as plt
import sys
from pathlib import Path
sys.path.append("C:/Users/jdabo/Desktop/")
import HCELL_readhdf5_boundary as HCELL_readhdf5
#import HCELL_readhdf5
#import HCELL_measure as measure
import numpy as np
import numpy.ma as ma

# Where does the data sit?
#datapath = "L:/hemocell/examples/short_AR1/output_0/hdf5/"
#datapath = "C:/Users/jdabo/MYSTUFF/STAGE_local/hemocell/examples/with_aneurysm/output/hdf5/"
datapath = "C:/Users/jdabo/Desktop/" # WINDOWS

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 1           # How many procs did you run on?
nprocs_preinlet = 5   # How many of those were devoted to the preinlet?
read_procs = 1        # Number of cores you want to read in with
timestep =     2500   # What was the time step?
TIME_begin = 227500   # When do you want to start reading the data
TIME_end =   230000   # When do you want to end reading the data.
DX = 5e-7;   DT = 1e-7

if len(sys.argv) > 1:
    datapath = sys.argv[1]
    for arg in str(sys.argv):
        if arg == '-r' or arg == '-R':
            REVERSE_X = True

REVERSE_X = True     # Measure in -x direction instead of x direction

#%%

# Read data
#fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, p=False, read_procs=read_procs, nprocs=nprocs, begin=TIME_begin, end=TIME_end, timestep=timestep, datapath=datapath)

fluid_pre, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, p=False, read_procs=read_procs, nprocs=nprocs_preinlet, begin=TIME_begin, end=TIME_end, timestep=timestep, datapath=datapath, fluidname='Fluid_PRE')


#%%


# Only look at the last time instance
fluid_last_time_instance = fluid_pre[-1]

# Compute average and maximum FLUID velocity
max_vel = -np.infty
running_sum = 0.0;
count = 0
for vel_triplet, bdry in zip(fluid_last_time_instance.velocity, \
                             fluid_last_time_instance.boundary):
    # If this cell is WITHIN boundary
    if bdry[0] == 0.0:
        vel_x = -vel_triplet[0] if REVERSE_X else vel_triplet[0]            
        running_sum += vel_x
        count += 1
        if vel_x > max_vel:
            max_vel = vel_x
conversion_factor = DX / DT
print("")
print("Max velocity  : " + str(max_vel * conversion_factor))
print("Mean velocity : " + str((running_sum / count) * conversion_factor))
