import sys
sys.path.append("L:/hemocell/scripts/measure/")  # 'default'
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import seaborn as sb;    sb.set()

# Where does the data sit?
datapath = "L:/hemocell/cases/AR1/output/hdf5/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 5              # How many procs did you run on?
read_procs = 4           # Number of cores you want to read in with
DX = 5e-7;   DT = 1e-7

timesteps = [10000000, 10000100]

#%%

# Read RBC data
_, RBC, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=timesteps, \
    datapath=datapath, rbcname='RBC_PRE')

#%%

prev_avg_positions = []
avg_vel_per_t = []
rbc_counts = []

first_iteration = True


for i_t, t in enumerate(timesteps):
    # This instance in time
    rbc = RBC[i_t]
    
    # Dictionary with lsp positions per RBC id
    lsp_pos_per_id = {}
    N_lsps_per_id = {}

    # positions and ids per LSP in this time instance
    positions = np.array(rbc.position)
    cids = np.array(rbc.cid)
    
    # Add together all LSP positions to a single sum of x,y,z coordinates (for averaging later)
    for pos_triplet, cid in zip(positions, cids):
        cid = cid[0]
        # Add entry to dictionary if it's not there yet
        if not cid in lsp_pos_per_id:
            lsp_pos_per_id[cid] = [0.0, 0.0, 0.0]
            N_lsps_per_id[cid]  = 0
            
        # Add LSP position and keep count of number of LSPs
        lsp_pos_per_id[cid] += pos_triplet
        N_lsps_per_id[cid]  += 1

    # Average LSP positions to get average position per RBC
    avg_pos_per_id = {}
    for cid, pos in lsp_pos_per_id.items():
        avg_pos_per_id[cid] = pos / N_lsps_per_id[cid]

    # Don't compute velocity for first time iteration
    if not first_iteration:
        # Velocity per RBC
        velocities = []
        
        rbc_count = 0
        
        # Itera1te through every RBC's position and compute velocity
        skipped = 0
        for cid, pos in avg_pos_per_id.items():
            try:
                delta_pos = pos - prev_avg_positions[cid]
                distance = np.linalg.norm(delta_pos)
                duration = t - timesteps[i_t - 1]
                velocities.append(distance / duration)                
                rbc_count += 1
            # RBCs need not be in both time instances. That's okay, just skip 'em
            except KeyError:
                pass

        # Compute average velocity over all RBCs
        rbc_counts.append(rbc_count)
        avg_vel_per_t.append(np.mean(velocities))

    # Remember average positions of last iteration
    prev_avg_positions = avg_pos_per_id
    prev_t = t
    
    if first_iteration:
        first_iteration = False
    

# Average over time, convert to mm/s, print
avg_rbc_count = np.mean(rbc_counts)
avg_vel_over_time = np.mean(avg_vel_per_t)
conversion_factor = DX / DT * 1000   #  * 1000 for m/s -> mm/s
print("Average number of RBCs measured per iteration:")
print(avg_rbc_count)
print("Average speed:")
print(avg_vel_over_time * conversion_factor, "mm/s")