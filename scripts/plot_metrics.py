#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:07:39 2016

@author: gzavo
"""

import numpy as np
import matplotlib.pyplot as plt
import process_out


names = ['iteration', 'wall-time (s)', 'largest force (lN)', 'mean velocity (m/s)', 'apparent rel. viscosity']
fnames = ['wall_time.png', 'largest_force.png', 'mean_vel.png', 'app_rel_visc.png']

def plotit(data, column):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(data[:,0], data[:,column], label=names[column])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.08))
    plt.xlabel(names[0])
    plt.ylabel(names[column])
    plt.savefig(fnames[column-1], dpi=200)
    


# Columns: 
#   0 iteration
#   1 wall time / iter
#   2 largest force
#   3 mean velocity
#   4 app. rel. visc.

if __name__ == "__main__":
    
    directory = './'
    
    process_out.process(directory)    
    data = np.loadtxt(directory + '/' + "metrics.dat")
    for i in range(len(fnames)):
        plotit(data, i+1)