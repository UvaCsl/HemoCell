#!/usr/bin/env python
import numpy as np
import pylab as pl
import h5py as h5
from glob import glob
import sys 
import matplotlib.pyplot as plot
import pylab


filename= sys.argv[1:]

for i in range(0, len(sys.argv[1:])):
        lab = str(filename[i])[:len(str(filename[i]))-4]
	#data = pylab.loadtxt(filename, skiprows=1)
	data=np.genfromtxt(filename[i], skip_header=1, skip_footer=1)
        bnddata= np.genfromtxt("bnd-"+filename[i], skip_header=1, skip_footer=1)
	datalen= len(data[:,0])


	Dist= ['None']*datalen
	Time= ['None']*datalen
	ForceZ = ['None']*datalen
	dX = ['None']*datalen

	Time1= ['None']*len(bnddata)
	Nbnd=['None']*len(bnddata)

	for i in range(datalen):
	    Dist[i]= -data[i,3]*1e6
	    Time[i]= data[i,0]*1e6
	    ForceZ[i] = data[i,4]*1e12
	    dX[i]= abs(data[i,1]-data[i,2])*1e6	
	  #  print Time[i], '\n'

   
        for i in range(len(bnddata)):
	    Nbnd[i]= bnddata[i,1]
	    Time1[i]= bnddata[i,0]*1e6

	

	plot.subplot(2, 1, 1)
	pylab.plot( Time, Dist,  label=lab)
	#pylab.plot( Time, dX,  label="(x0-x1)" )

	pylab.ylim((0,2))
	pylab.xlim((0,20000))
	pylab.title("Platelet disaggregation dynamics")
	pylab.xlabel("Time (us)")
	pylab.ylabel("Distance (um)")

	plot.subplot(2, 1, 2)
	pylab.plot( Time1, Nbnd,  label=lab )
	#pylab.legend()
	pylab.xlim((0,20000))
	#pylab.title("Number on bonds in time")
	pylab.xlabel("Time (us)")
	pylab.ylabel("Number of bonds")

#pylab.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plot.show()






