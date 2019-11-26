import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb;   sb.set();   sb.set_style("whitegrid")
from matplotlib import lines as mlines


# Read data and sort
data = pd.read_csv("C:/Users/jdabo/Dropbox/STAGE/results/measurements.csv", delimiter = " ", index_col = False) 
data.sort_values(by=["region", "celltype"], axis=0, inplace=True)

# ===============================================
# NECK WIDTH
# ===============================================

# Use list comprehension to make data subset, only RBCS/PLTS without stiff
RBC_subset = [True if b1 == b2 == True else False for b1,b2 in \
              zip(data["celltype"] == "RBC", data["case"] != "AR2_stiff")]
PLT_subset = [True if b1 == b2 == True else False for b1,b2 in \
              zip(data["celltype"] == "PLT", data["case"] != "AR2_stiff")]
RBCS_no_stiffs = data[RBC_subset]
PLTS_no_stiffs = data[PLT_subset]

# Some bookkeeping variables
cases = ["AR1", "AR2", "AR3"];                            nCases   = len(cases)
regions = ["an. sac", "outlet", "vessel", "preinlet"];    nRegions = len(regions)
AR_colors = ["#003F5C", "#7A5195", "#EF5675"]
cellNames = ["RBC", "platelet"]
FONTSIZE = 32
FONTSIZE_LEGEND = 28
FONTSIZE_ANNOTATIONS = 28

# Make a plot separately for RBCs and platelets
for df, cellName in zip([RBCS_no_stiffs, PLTS_no_stiffs], cellNames):
    fig, ax = plt.subplots(figsize=(13,10))
    
    # Plot the bars
    count = 0
    for i, row in df.iterrows():
        # set x position so that some space is left after each group of  (=nCases)
        pos = count + int(count / nCases)
        vel = row['velocity']
        avgSD = row['avgSD']
        SD_t = row['SD_time']
        nCells = row['nCells']
        # Label the first 4 (for legend
        if count < nCases:
            ax.bar(pos, height = vel, align='center', \
                   color=AR_colors[count], label=cases[count])
        # After that just plot without labeling
        else:
            ax.bar(pos, height = vel, align='center', \
                   color=AR_colors[count % nCases])
        # Errorbars
        ax.errorbar(pos, vel, yerr=[[0], [avgSD]], capsize=10, \
                    barsabove=True, color='black')
        ax.errorbar(pos, vel, yerr=[[avgSD], [0]], capsize=10, \
                    barsabove=True, color='black')
        # Annotations
        if vel > 1.0:
            annotationStartsAtThisY, color, weight = 0.16, 'white', 'bold'
        else:
            annotationStartsAtThisY, color, weight = 1.25, 'black', 'bold'
        ax.annotate("n = " + str(nCells), (pos, annotationStartsAtThisY), weight=weight,
                    rotation=90, color=color, fontsize = FONTSIZE_ANNOTATIONS,  ha="center", va="bottom")
                    
        count += 1
    
    # set ticks, set their fontsize
    xticks = [(nCases+1)*n + (nCases / 2 - 0.5) for n in range(nRegions)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(regions, fontsize = FONTSIZE)
    yticks = ax.get_yticks()
    ax.set_yticklabels(yticks, fontsize = FONTSIZE)

    # Shrink current axis's width by 20% on the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width * 0.8, box.height])
    
    # Put a bar legend to the right of axis
    legend1 = ax.legend(loc='upper left',
              fancybox=True, shadow=True, ncol=1, fontsize = FONTSIZE_LEGEND)

    plt.gca().add_artist(legend1)

    # Title and show
    #plt.title("Average " + cellName + " velocity by region and Aspect Ratio", fontsize = FONTSIZE)
    plt.xlabel("region", fontsize = FONTSIZE)
    plt.ylabel("Velocity (mm/s)", fontsize = FONTSIZE)
    plt.show()

#%%
# ===============================================
# REGULAR VERSUS STIFF
# ===============================================
    
cases = ["AR2", "AR2 stiff"];    nCases = len(cases)
normal_stiff_colors = ["crimson", "darkblue"]
         
# Use list comprehension to make data subset, only RBC/PLT of AR2 / AR2_stiff cases
RBC_subset = [True if b1 == True and (b2 == True or b3 == True) else False for b1,b2,b3 in \
              zip(data["celltype"] == "RBC", data["case"] == "AR2", data["case"] == "AR2_stiff")]
PLT_subset = [True if b1 == True and (b2 == True or b3 == True) else False for b1,b2,b3 in \
              zip(data["celltype"] == "PLT", data["case"] == "AR2", data["case"] == "AR2_stiff")]
RBCS = data[RBC_subset]
PLTS = data[PLT_subset]

# Some bookkeeping variables
regions = ["pre-inlet", "vessel", "an. sac", "outlet"];    nRegions = len(regions)
cellNames = ["RBC", "platelet"]

# hacky
pos_to_new_pos = {0 : 6, 1 : 7, 3: 9, 4: 10, 6 : 3, 7 : 4, 9 : 0, 10 : 1}

# Make a plot separately for RBCs and platelets
for df, cellName in zip([RBCS, PLTS], cellNames):
    fig, ax = plt.subplots(figsize=(14,8))
    
    # Plot the bars
    count = 0
    for i, row in df.iterrows():
        # set x position so that some space is left after each group of  (=nCases)
        pos = count + int(count / nCases)
        pos = pos_to_new_pos[pos]
        vel = row['velocity']
        avgSD = row['avgSD']
        SD_t = row['SD_time']
        nCells = row['nCells']
        # Label the first 4 (for legend
        if count < nCases:
            ax.bar(pos, height = vel, align='center', \
                   color=normal_stiff_colors[count], label=cases[count])
        # After that just plot without labeling
        else:
            ax.bar(pos, height = vel, align='center', \
                   color=normal_stiff_colors[count % nCases])
        # Errorbars
        ax.errorbar(pos, vel, yerr=[[0], [avgSD]], capsize=10, \
                    barsabove=True, color='black')
        ax.errorbar(pos, vel, yerr=[[avgSD], [0]], capsize=10, \
                    barsabove=True, color='black')
        # Annotations
        if vel > 1.0:
            annotationStartsAtThisY, color, weight = 0.18, 'white', 'bold'
        else:
            annotationStartsAtThisY, color, weight = 1.25, 'black', 'bold'
        ax.annotate("n = " + str(nCells), (pos, annotationStartsAtThisY), weight=weight,
                    rotation=90, color=color, fontsize = FONTSIZE_ANNOTATIONS,  ha="center", va="bottom")
                    
        count += 1
    
    # set ticks, set their fontsize
    xticks = [(nCases+1)*n + (nCases / 2 - 0.5) for n in range(nRegions)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(regions, fontsize = FONTSIZE)
    yticks = ax.get_yticks()
    ax.set_yticklabels(yticks, fontsize = FONTSIZE)

    # Shrink current axis's width by 20% on the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width * 0.8, box.height])
    
    # Put a bar legend to the right of axis
    legend1 = ax.legend(loc='upper right',
              fancybox=True, shadow=True, ncol=1, fontsize = FONTSIZE_LEGEND)

    plt.gca().add_artist(legend1)

    # Title and show
    #plt.title("Average " + cellName + " velocity by region and healty / stiff", fontsize = FONTSIZE)
    plt.xlabel("region", fontsize = FONTSIZE)
    plt.ylabel("Velocity (mm/s)", fontsize = FONTSIZE)
    plt.show()