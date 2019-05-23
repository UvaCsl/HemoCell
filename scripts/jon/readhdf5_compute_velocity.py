import sys
sys.path.append("L:/hemocell/scripts/measure/")  # 'default'
import HCELL_readhdf5_jon as HCELL_readhdf5
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb;    sb.set()

# Where does the data sit?
datapath = "L:/hemocell/cases/AR2_stiff/output/hdf5/"

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 5              # How many procs did you run on?
read_procs = 4           # Number of cores you want to read in with
DX = 5e-7;   DT = 1e-7

timesteps = [10000000, 10000100]

REVERSE_X = True
SAVE_FIGS = False

    #%%

# Read data
#fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, p=False, read_procs=read_procs, nprocs=nprocs, begin=TIME_begin, end=TIME_end, timestep=timestep, datapath=datapath)

fluid_pre, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, timesteps=timesteps, \
    datapath=datapath, fluidname='Fluid_PRE')

#%%
fluid = fluid_pre[-1]  # (look at last time instance)

#%%
pos = np.array(fluid.position)
vel = np.array(fluid.velocity)
bdrs = np.array(fluid.boundary)

#%%

def removeDuplicates(pos, vel, bdrs):
    pos, indices = np.unique(pos, return_index = True, axis=0)
    vel = vel[indices]
    bdrs = bdrs[indices]
    return pos, vel, bdrs

def printFluidVelocity(pos, vel, bdrs):
    # Compute average and maximum velocity
    max_vels = {}
    running_sum = 0.0;
    count = 0
    for pos_triplet, vel_triplet, bdry in zip(pos, vel, bdrs):
        # If this cell is WITHIN geometry boundary
        if bdry[0] == 0.0:
            # Index [0] is for gettin velocity in x-direction only
            vel_x = -vel_triplet[0] if REVERSE_X else vel_triplet[0]            
            running_sum += vel_x
            count += 1
            # See if this max is highest for this slice
            x = pos_triplet[0]
            if not x in max_vels:
                max_vels[x] = vel_x
            elif vel_x > max_vels[x]:
                max_vels[x] = vel_x
            
    conversion_factor = (DX / DT) * 1000   # *1000 for m/s -> mm/s
    print("")
    print("Median max fluid velocity mm/s: " + str(np.median(list(max_vels.values())) * conversion_factor))
    print("Mean       fluid velocity mm/s: " + str((running_sum / count) * conversion_factor))
    
# Get the x,y,z boundaries of geometry bounding box
def getSpatialBoundaries(pos):
    xmin, ymin, zmin =  np.infty,  np.infty,  np.infty
    xmax, ymax, zmax = -np.infty, -np.infty, -np.infty
    for pos_triplet in np.array(pos):
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
   
def computeVelocityXaverage(pos, vel, xmin, xmax, ymin, ymax, zmin, zmax):
    # We're going to average over x dimension
    x_range = int(round(xmax - xmin)) + 1;
    y_range = int(round(ymax - ymin)) + 1;
    z_range = int(round(zmax - zmin)) + 1;   
    average_slice = np.zeros((y_range, z_range))
    
    # Make a dictionary that translates y,z coordinates to corresponding indices
    pos_to_index = {}
    y = ymin
    for i in range(average_slice.shape[0]):
        z = zmin
        for j in range(average_slice.shape[1]):
            pos_to_index[(y,z)] = (i,j)
            z += 1.0
        y += 1.0    
    
    # Compute average x-velocity slice where the slice has y/z dimensions (so average over x dimension)
    print("Computing x-averaged velocity slice...")
    print("Percentage: ", end = ' ')
    length = float(len(fluid.velocity))
    count = 0;    already_printed = []
    for pos_triplet, vel_triplet in zip(pos, vel):
        percentage = int(round((count / length) * 100))
        if not percentage in already_printed and percentage % 10 == 0:
            print(percentage, end=' ')
            already_printed.append(percentage)

        (y,z) = pos_triplet[1:]
        average_slice[pos_to_index[(y,z)]] += -vel_triplet[0] if REVERSE_X else vel_triplet[0]
            
        count += 1
    average_slice /= x_range  # Average
    
    print("")
    return average_slice


def plotFluidVelocityHeatmap(average_slice):
    # Convert values from LU to SI for plot
    conversion_factor = DX / DT * 1000   # *1000 for m/s -> mm/s
    average_slice_convert = average_slice * conversion_factor
    
    # Plot
    plt.figure(0)
    plt.imshow(average_slice_convert, cmap="hot")
    plt.title("x-averaged fluid velocity in mm/s.")
    plt.xlabel(r"Y position in $\mu$m.")
    plt.ylabel(r"Z position in $\mu$m.")
    ax = plt.gca()
    ticks = [0, 10, 20, 30]
    spatial_conv_factor = 0.5
    labels = [str(tick * spatial_conv_factor) for tick in ticks]
    ax.set_xticks(ticks[:-1])
    ax.set_xticklabels(labels[:-1])
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    plt.colorbar()
    if not SAVE_FIGS:
        plt.show()
    else:
        plt.savefig('vel_heatmap.png')

# Plot fluid velocity versus radius profile
def plotFluidVelocityVSradius(average_slice, ymin, ymax, zmin, zmax):    

    # This function returns distance to center
    def compute_distance(Y,Z,middle):
        dist = np.linalg.norm(np.array([[Y - middle[0], 0],[0, Z - middle[1]]]))
        if Y < middle[0]:
            dist *= -1
        return dist
    
    # Compute y/z center of slice    
    middle = (ymin + (ymax - ymin) / 2.0, zmin + (zmax - zmin) / 2.0)
    fluid_vel_radius = {}    
    for y, row in enumerate(average_slice):
        for z, vel in enumerate(row):
            if vel > 0.0:  # We don't care about 0-velocity cells outside of boundary
                Y = y + ymin;   Z = z + zmin
                dist = compute_distance(Y,Z,middle)
                fluid_vel_radius[dist] = vel

    Xplot = sorted(fluid_vel_radius.keys())
    Yplot = [fluid_vel_radius[key] for key in Xplot]

    # Convert values from LU to SI for plot
    conversion_factor = DX / DT * 1000   # *1000 for m/s -> mm/s
    Yplot = [vel * conversion_factor for vel in Yplot]
    
    # 2 LU is one micron
    Xplot = [length * 0.5 for length in Xplot]

    plt.figure(1)
    #sb.scatterplot(Xplot, Yplot)
    plt.plot(Xplot, Yplot, "r.")
    plt.title("Fluid velocity versus distance from center")
    plt.ylabel("Velocity in mm/s")
    plt.xlabel(r"Distance from center of pipe in $\mu$m.")
    plt.grid(True)
    if not SAVE_FIGS:
        plt.show()
    else:
        plt.savefig('vel_vs_radius.png')

#%%
pos, vel, bdrs = removeDuplicates(pos, vel, bdrs)
#%%
printFluidVelocity(pos, vel, bdrs)
#%%
doTheRest = False

# Check if we want to do all the other plotting stuff	
if doTheRest:
	xmin, xmax, ymin, ymax, zmin, zmax = getSpatialBoundaries(pos)

	average_slice = computeVelocityXaverage(pos, vel, xmin, xmax, ymin, ymax, zmin, zmax)

#%%
	plotFluidVelocityHeatmap(average_slice)

	plotFluidVelocityVSradius(average_slice, ymin, ymax, zmin, zmax)