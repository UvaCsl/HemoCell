import sys
sys.path.append("L:/hemocell/scripts/jon/")  # 'default'
import matplotlib
from matplotlib import pyplot as plt        
import numpy as np
import pickle
import csv

#%%

SAVE_FIGS = False

def plotSideview(sideview, geo_slice, N_t):
    # First create two color maps. One for the geometry contour and one with 
    # transparency for the overlaid RBC heat map
    from matplotlib.colors import LinearSegmentedColormap

    # Color map for geometry
    N = 256
    geometry_color = (200.0, 200.0, 200.0)
    geometry_color = [c/N for c in geometry_color]
    geometry_color.append(1.0)  # no transparency
    white_c = [1.0, 1.0, 1.0, 1.0]
    
    geometry_color_array = np.array([geometry_color, white_c])
    map_object = LinearSegmentedColormap.from_list(name='geo_map', colors=geometry_color_array)
    plt.register_cmap(cmap=map_object)
    
    # Comment here
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], 1.0, N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], 0.0, N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], 0.0, N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)
    
    # Changed sideview so that it's a fraction of timesteps
    sideview /= N_t
    
    # Mirror both matrices
    geo_slice = np.array(geo_slice[::-1])
    sideview  = np.array(sideview[::-1])
    
    # Transform sideview data with log for log scaling
    TwoPlusLog = lambda v : 2 + np.log10(v) if v > 0.0 else 0.0
    sideview = np.array([[TwoPlusLog(v) for v in row] for row in sideview])
    
    # Plot
    FONTSIZE = 24
    
    matplotlib.rc('xtick', labelsize=FONTSIZE) 
    matplotlib.rc('ytick', labelsize=FONTSIZE) 
    fig, ax = plt.subplots(1, 1)
    ax.imshow(geo_slice.transpose(), cmap='geo_map')
    fig2 = ax.imshow(sideview.transpose(), cmap='RBC_colors')
    colorbar = plt.colorbar(fig2)
    colorbar.ax.set_ylabel('Avg # of RBCs present over time', rotation=270, labelpad=40, fontsize=FONTSIZE)
    colorbar.set_ticks([0.0, 0.5, 1.0, 1.5])
    colorbar.set_ticklabels(["0.01", r"10$^{-1.5}$", "0.1", r"10$^{-0.5}$"])
    ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    spatial_conv_factor = 0.5
    xlabels = [str(tick * spatial_conv_factor) for tick in xticks]
    ylabels = [str(tick * spatial_conv_factor) for tick in yticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)    
    ax.set_xlim([10, 400])
    ax.set_ylim([210, 0])
    plt.grid(False)
    plt.title("AR2 RBC location averaged over time", fontsize=FONTSIZE)
    plt.xlabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    plt.ylabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    if SAVE_FIGS:
        matplotlib.use('Agg')
        plt.savefig("RBCs_pos_over_time")
    else:
        matplotlib.use('Qt5Agg')
        plt.show()

def plotDifference(diff1, diff2, geo_slice, N_t):
    # First create two color maps. One for the geometry contour and one with 
    # transparency for the overlaid RBC heat map
    from matplotlib.colors import LinearSegmentedColormap

    # Color map for geometry
    N = 256
    geometry_color = (200.0, 200.0, 200.0)
    geometry_color = [c/N for c in geometry_color]
    geometry_color.append(1.0)  # no transparency
    white_c = [1.0, 1.0, 1.0, 1.0]
    #red_c = [1.0, 0.0, 0.0, 1.0]
    green_c = [0.0, 1.0, 0.0, 1.0]
    blue_c = [0.0, 0.0, 1.0, 1.0]
    
    
    geometry_color_array = np.array([geometry_color, white_c])
    map_object = LinearSegmentedColormap.from_list(name='geo_map', colors=geometry_color_array)
    plt.register_cmap(cmap=map_object)
    
    # Create color map for normal RBCs
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], green_c[0], N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], green_c[1], N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], green_c[2], N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)   # transparency
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)

    # Create color map for normal RBCs
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], blue_c[0], N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], blue_c[1], N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], blue_c[2], N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)   # transparency
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_stiff_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)
    
    # Mirror both matrices
    geo_slice = np.array(geo_slice[::-1])
    diff1  = np.array(diff1[::-1])
    diff2  = np.array(diff2[::-1])

    # Divide through time to get fraction of time
    diff1 /= N_t;   diff2 /= N_t
    
    # Transform sideview data with log for log scaling
    TwoPlusLog = lambda v : 2 + np.log10(v) if v > 0.0 else 0.0
    diff1 = np.array([[TwoPlusLog(v) for v in row] for row in diff1])
    diff2 = np.array([[TwoPlusLog(v) for v in row] for row in diff2])
    
    # Plot
    FONTSIZE = 24
    
    matplotlib.rc('xtick', labelsize=FONTSIZE) 
    matplotlib.rc('ytick', labelsize=FONTSIZE) 
    fig, ax = plt.subplots(1, 1)
    ax.imshow(geo_slice.transpose(), cmap='geo_map')
    fig2 = ax.imshow(diff1.transpose(), cmap='RBC_stiff_colors')
    fig3 = ax.imshow(diff2.transpose(), cmap='RBC_colors')
    colorbar = plt.colorbar(fig3)
    colorbar.ax.set_ylabel('Avg # of RBCs present over time', rotation=270, labelpad=40, fontsize=FONTSIZE)
    colorbar.set_ticks([0.0, 0.5, 1.0, 1.5])
    colorbar.set_ticklabels(["0.01", r"10$^{-1.5}$", "0.1", r"10$^{-0.5}$"])
    ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    spatial_conv_factor = 0.5
    xlabels = [str(tick * spatial_conv_factor) for tick in xticks]
    ylabels = [str(tick * spatial_conv_factor) for tick in yticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)    
    ax.set_xlim([10, 400])
    ax.set_ylim([210, 0])
    plt.grid(False)
    plt.title("AR2 RBC location averaged over time", fontsize=FONTSIZE)
    plt.xlabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    plt.ylabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    if SAVE_FIGS:
        matplotlib.use('Agg')
        plt.savefig("RBCs_pos_over_time")
    else:
        matplotlib.use('Qt5Agg')
        plt.show()

#%%        
CASE = "AR1"
RESULTS_PATH = "L:/no_backup/" + CASE + "/output/heatmap/"

boundaries_path = RESULTS_PATH + CASE + "_spatial_boundaries.pickle"
geo_slice_path = RESULTS_PATH + CASE + "_geometry_slice.pickle"
        
xmin, xmax, ymin, ymax, zmin, zmax = pickle.load(open(boundaries_path, "rb"))
geo_slice = pickle.load(open(geo_slice_path, "rb"))

#%%
sideview_path = RESULTS_PATH + CASE + "_sideview.csv"

with open(sideview_path, newline='\r\n') as f:
  reader = csv.reader(f, delimiter = ' ')
  row1 = next(reader)  # gets the first line
  N_t = int(row1[1])

sideview = np.loadtxt(open(sideview_path, "rb"), delimiter=" ", skiprows=1)

#absdiff = np.zeros(sideview1.shape)
#diff1 = np.zeros(sideview1.shape)
#diff2 = np.zeros(sideview1.shape)
#I,J = sideview1.shape
#for i in range(I):
#    for j in range(J):
#        absdiff[i][j] = abs(sideview2[i][j] - sideview1[i][j])
#        diff1[i][j] = sideview2[i][j] - sideview1[i][j]
#        diff2[i][j] = sideview1[i][j] - sideview2[i][j]

plotSideview(sideview, geo_slice, N_t)
#plotSideview(diff1, geo_slice, N_t)
#plotSideview(diff2, geo_slice, N_t)
#plotDifference(diff1, diff2, geo_slice, N_t)

#%%