import sys
sys.path.append("C:/Users/jdabo/Desktop/")  # 'default'
# Command-line arguments
args = sys.argv
for i, arg in enumerate(args):
    if arg == "-l":      # Readhdf5 library
        sys.path.append(args[i+1])
from matplotlib import pyplot as plt        
import HCELL_readhdf5_boundary as HCELL_readhdf5
import numpy as np
import numpy.ma as ma

# Where does the data sit?
datapath = "C:/Users/jdabo/Desktop/hdf5/"

DX = 5e-7;   DT = 1e-7;

#FIRST SPECIFY SOME STUFF ABOUT YOUR DATA
nprocs = 10           # -n How many procs did you run on?
read_procs = 4        # -d Number of cores you want to read in with
timestep =     2500   # -tstep What was the time step?
TIME_begin = 97500   # -tmin When do you want to start reading the data
TIME_end =   100000   # -tmax When do you want to end reading the data.
REVERSE_X = True     # Simulation is going from x -> -x or reverse

SAVE_FIGS = False
REVERSE_X = True
for i, arg in enumerate(args):
    if arg == "-d":      # data
        datapath = args[i+1]
    if arg == "-tmin":   # time first step (inclusive)
        TIME_begin = int(args[i+1])
    if arg == "-tmax":   # time last step (exclusive I think)
        print(args[i+1])        
        TIME_end = int(args[i+1])        
    if arg == "-tstep":  # time step
        timestep = int(args[i+1])
    if arg == "-n":      # Number of processors the sim was run on
        nprocs = int(args[i+1])
    if arg == "R":
        REVERSE_X = True
    if arg == "save_figs":
        SAVE_FIGS = True

#%%

# Read data
#fluid, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, p=False, read_procs=read_procs, nprocs=nprocs, begin=TIME_begin, end=TIME_end, timestep=timestep, datapath=datapath)

fluid_pre, _, _, _ = HCELL_readhdf5.open_hdf5_files(ct3=False, half=True, r=False, \
    p=False, read_procs=read_procs, nprocs=nprocs, begin=TIME_begin, \
    end=TIME_end, timestep=timestep, datapath=datapath, fluidname='Fluid_PRE')

#%%
fluid = fluid_pre[-1]  # (look at last time instance)

#%%

def removeDuplicates(fluid):
    pos, indices = np.unique(fluid.position, return_index = True, axis=0)
    vel = np.array(fluid.velocity)[indices]
    bdry = np.array(fluid.boundary)[indices]
    fluid.position = pos
    fluid.velocity = vel
    fluid.boundary = bdry
    return fluid

def printFluidVelocity(fluid):
    # Compute average and maximum velocity
    max_vel = -np.infty
    running_sum = 0.0;
    count = 0
    for pos_triplet, vel_triplet, bdry in zip(fluid.position, fluid.velocity, fluid.boundary):
        # If this cell is WITHIN geometry boundary
        if bdry[0] == 0.0:
            # Index [0] is for gettin velocity in x-direction only
            vel_x = -vel_triplet[0] if REVERSE_X else vel_triplet[0]            
            running_sum += vel_x
            count += 1
            if vel_x > max_vel:
                max_vel = vel_x
    conversion_factor = DX / DT
    print("")
    print("Max  fluid velocity mm/s: " + str(max_vel * conversion_factor * 1000))  # *1000 for m/s -> mm/s
    print("Mean fluid velocity mm/s: " + str((running_sum / count) * conversion_factor * 1000))
    
# Get the x,y,z boundaries of geometry bounding box
def getSpatialBoundaries(fluid):
    xmin, ymin, zmin =  np.infty,  np.infty,  np.infty
    xmax, ymax, zmax = -np.infty, -np.infty, -np.infty
    for pos_triplet in np.array(fluid.position):
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
   
def computeVelocityXaverage(fluid, xmin, xmax, ymin, ymax, zmin, zmax):
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
    
    # Keeps track of coordinates already done; skip double values
    already_done = []
    
    # Compute average x-velocity slice where the slice has y/z dimensions (so average over x dimension)
    print("Computing x-averaged velocity slice...")
    print("Percentage: ", end = ' ')
    length = float(len(fluid.velocity))
    count = 0;    already_printed = []
    for pos_triplet, vel_triplet, bdry in zip(fluid.position, fluid.velocity, fluid.boundary):
        percentage = int(round((count / length) * 100))
        if percentage % 10 == 0:
            if not percentage in already_printed:
                print(percentage, end=' ')
                already_printed.append(percentage)

        (y,z) = pos_triplet[1:]
        average_slice[pos_to_index[(y,z)]] += -vel_triplet[0] if REVERSE_X else -vel_triplet[0]
        already_done.append(list(pos_triplet))
            
        count += 1
    average_slice /= x_range  # Average
    
    print(" ")
    return average_slice


def plotFluidVelocityHeatmap(average_slice):
    # Convert values from LU to SI for plot
    conversion_factor = DX / DT * 1000   # *1000 for m/s -> mm/s
    average_slice_convert = average_slice * conversion_factor
    
    # Plot
    plt.figure(0)
    plt.imshow(average_slice_convert, cmap='hot')
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
fluid = removeDuplicates(fluid)
#%%
printFluidVelocity(fluid)
#%%
# Check if we want to do all the other plotting stuff
doTheRest = True
for arg in args:
	if arg == "stats":
		doTheRest = False
	
if doTheRest:
	xmin, xmax, ymin, ymax, zmin, zmax = getSpatialBoundaries(fluid)

	average_slice = computeVelocityXaverage(fluid, xmin, xmax, ymin, ymax, zmin, zmax)

	plotFluidVelocityHeatmap(average_slice)

	plotFluidVelocityVSradius(average_slice, ymin, ymax, zmin, zmax)