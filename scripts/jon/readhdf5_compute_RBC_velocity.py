import sys
sys.path.append("C:/Users/jdabo/MYSTUFF/STAGE_local/hemocell/scripts/measure/")  # 'default'
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import seaborn as sb
sb.set()

# Where does the data sit?
datapath = r"L:/no_backup/AR2/output/hdf5/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 40              # How many procs did you run on?
read_procs = 4           # Number of cores you want to read in with
DX = 5e-7;   DT = 1e-7

# Define time steps. Add t+1 at every t for computing velocity from position
#T_STEP  =         1
#T_START =  10000000
#T_END   =  10000010
#N_steps =        11
#timesteps = list(np.linspace(T_START, T_END, N_steps, endpoint = True))
#timesteps = sorted(timesteps.append([t+1 for t in timesteps]))

timesteps = [5000000, 7500000, 10000000]

#%%

# Read RBC data
_, RBC, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, f=False, r=True, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=timesteps, \
    datapath=datapath, rbcname='RBC')

#%%

prev_avg_positions = []
avg_vel_per_t = []

first_iteration = True

for i_t, t in enumerate(timesteps[:-1]):    # Don't compute velocity for last time step
    # This instance in time
    rbc = RBC[i_t]
    
    # Dictionary with lsp positions per RBC id
    lsp_pos_per_id = {}
    N_lsps_per_id = {}

    # positions and ids per LSP in this time instance
    positions = np.array(rbc.position)
    cids = np.array(rbc.cid)
    
    # Add together all LSP positions to a single sum of x,y,z coordinates (for averaging later)
    print("c1")
    length = len(rbc.position)
    count = 0
    printed = []
    for pos_triplet, cid in zip(positions, cids):
        cid = cid[0]
        count += 1
        bleh = int((count / length) * 10.0)
        if not bleh in printed:
            print(bleh, end=" ")
            printed.append(bleh)
        # Add entry to dictionary if it's not there yet
        if not cid in lsp_pos_per_id:
            lsp_pos_per_id[cid] = [0.0, 0.0, 0.0]
            N_lsps_per_id[cid]  = 0
            
        # Add LSP position and keep count of number of LSPs
        lsp_pos_per_id[cid] += pos_triplet
        N_lsps_per_id[cid]  += 1

    # Average LSP positions to get average position per RBC
    print("\nc2")
    avg_pos_per_id = {}
    for cid, pos in lsp_pos_per_id.items():
        avg_pos_per_id[cid] = pos / N_lsps_per_id[cid]

    # Remember previous positions or compute velocity depending on even/odd iteration
    if not first_iteration:
        # Velocity per RBC
        velocities = []
        
        # Iterate through every RBC's position and compute velocity
        for cid, pos in avg_pos_per_id.items():
            try:
                delta_pos = pos - prev_avg_positions[cid]
                distance = np.linalg.norm(delta_pos)
                duration = t - timesteps[i_t - 1]
                velocities.append(distance / duration)                
            # RBCs need not be in both time instances. That's okay, just skip 'em
            except KeyError:
                pass

        # Compute average velocity over all RBCs
        avg_vel_per_t.append(np.mean(velocities))

    print("c3")
    # Remember average positions of last iteration
    prev_avg_positions = avg_pos_per_id
    prev_t = t
    
    if first_iteration:
        first_iteration = False
    

# Average over time, convert to mm/s, print
avg_vel_over_time = np.mean(avg_vel_per_t)
conversion_factor = DX / DT * 1000   #  * 1000 for m/s -> mm/s
print("test!")
print(avg_vel_over_time * conversion_factor)