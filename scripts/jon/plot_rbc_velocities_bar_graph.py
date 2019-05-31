import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sb;   sb.set()

# Function for annotating each bar
def annotateBars(row, ax):
    ax.text()

# Read data and sort
data = pd.read_csv("measurements.csv", delimiter = " ", index_col = False) 
data.sort_values(by=["region", "celltype"], axis=0, inplace=True)

# Data subsets based on celltype
RBC = data[data["celltype"] == "RBC"]
PLT = data[data["celltype"] == "PLT"]

# Some bookkeeping variables
cases = ["AR1", "AR2", "AR2 stiff", "AR3"];    nCases = len(cases)
regions = ["geometry", "inlet", "outlet", "preinlet"]
colors = ["#003F5C", "#7A5195", "#EF5675", "#ffA600"]
cellNames = ["RBC", "platelet"]
FONTSIZE = 24

# Make a plot separately for RBCs and platelets
for df, cellName in zip([RBC, PLT], cellNames):
    fig, ax = plt.subplots(figsize=(20,10))
    
    # Plot the bars
    count = 0
    for i, row in df.iterrows():
        # set x position so that some space is left after each group of 4 (=nCases)
        pos = count + int(count / nCases)
        vel = row['velocity']
        avgSD = row['avgSD']
        SD_t = row['SD_time']
        nCells = row['nCells']
        # Label the first 4 (for legend)
        if count < nCases:
            ax.bar(pos, height = vel, align='center', \
                   color=colors[count], label=cases[count])
        # After that just plot without labeling
        else:
            ax.bar(pos, height = vel, align='center', \
                   color=colors[count % nCases])
        # Errorbars
        ax.errorbar(pos, vel, yerr=[[0], [avgSD]], capsize=10, \
                    barsabove=True, color='red')
        ax.errorbar(pos, vel, yerr=[[SD_t], [0]], capsize=10, \
                    barsabove=True, color='black')
        # Annotations
        annotationEndsAtThisY = 0.2 if cellName == "RBC" else 0.05
        ax.annotate("n = " + str(nCells), (pos, annotationEndsAtThisY), \
                    rotation=90, color='white', fontsize = 18,  ha="center", va="bottom")
                
        count += 1
    
    # set ticks, set their fontsize
    xticks = [5*n + 1.5 for n in range(4)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(regions, fontsize = FONTSIZE)
    yticks = ax.get_yticks()
    ax.set_yticklabels(yticks, fontsize = FONTSIZE)

    # Apply annotations
    #data.apply(annotateBars, ax=ax, axis=1)

    # Shrink current axis's width by 20% on the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width * 0.8, box.height])
    
    
    #legend1 = pyplot.legend(plot_lines[0], ["algo1", "algo2", "algo3"], loc=1)
    #pyplot.legend([l[0] for l in plot_lines], parameters, loc=4)
    
    # Put a bar legend to the right of axis
    legend1 = ax.legend(loc='upper center', bbox_to_anchor=(1.11, 0.725),
              fancybox=True, shadow=True, ncol=1, fontsize = FONTSIZE)
    
    # Second error bar legend
    from matplotlib import lines as mlines
    red_line   = mlines.Line2D([], [], color='red', markersize=15, label='avg SD')
    black_line = mlines.Line2D([], [], color='black', markersize=15, label='SD over time')
    
    ax.legend(handles=[red_line, black_line], loc='upper center', bbox_to_anchor=(1.139, 0.45),
              fancybox=True, shadow=True, ncol=1, fontsize = FONTSIZE)

    plt.gca().add_artist(legend1)

    # Title and show
    plt.title("Average " + cellName + " velocity by region and case", fontsize = FONTSIZE)
    plt.show()
