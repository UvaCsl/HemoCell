# coding: utf-8

# In[ ]:

import numpy as np
import pylab as pl

# Load experimental data
DA = np.loadtxt('./Stretching-FedosovCaswell2010b.DA.dat')
DA_std = np.loadtxt('./Stretching-FedosovCaswell2010b.DA-STD.dat')

DT = np.loadtxt('./Stretching-FedosovCaswell2010b.DT.dat')
DT_std = np.loadtxt('./Stretching-FedosovCaswell2010b.DT-STD.dat')

# Load ficsion data
data = np.loadtxt('./cases_rbcModel-0.dat', delimiter=',')


# In[ ]:

# Plot results
pl.errorbar(DA[:,0], DA[:,1], yerr=DA_std[:,1]-DA[:,1], fmt='rd')
pl.errorbar(DT[:,0], DT[:,1], yerr=DT_std[:,1]-DT[:,1], fmt='rd')

pl.plot(data[:,0], data[:,1], 'bo-')
pl.plot(data[:,0], data[:,3], 'bo-')


# In[ ]:



