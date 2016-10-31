#!/usr/bin/env python
import numpy as np
import pylab as pl
import h5py as h5
from glob import glob
import sys 
import matplotlib.pyplot as plot
import pylab

from optparse import OptionParser
""" Parse Options """
parser = OptionParser() 
parser.add_option("-I", "--inputdir", dest="inputdir", default='./tmp',
                  help="Directory where input files are stored (without '/hdf5/')[default:'./tmp']")            
parser.add_option("-O", "--outputfile", dest="outputfile", default='./test.npy',
                help="Directory where the output will be stored [default: './test.npy' ]")
parser.add_option("-N", "--numberOfParticles", dest="Np", default=2,
                help="Number of platelets [default: '2' ]")


(options, args) = parser.parse_args()
caseFolder = options.inputdir + '/hdf5/'
saveFile = options.outputfile
Np = int(options.Np)


allFiles = sorted(glob(caseFolder + '/*PLT_*.h5'))
fnames = filter(lambda x: '.h5' in x, allFiles)
Nt = len(fnames)



keys2get = ['Force', 'Max positions', 'Min positions', 'Position']
dtype = np.dtype([('Force', '<f4', (3,)), ('Max positions', '<f4', (3,)), ('Min positions', '<f4', (3,)), ('Position', '<f4', (3,)), ('t', '<f4', (1,))   ])

ra = np.zeros((Nt, Np), dtype=dtype)

h = h5.File(fnames[0])
dt = h.attrs['dt'][0]
dx = h.attrs['dx'][0]
df = 1025 * dx**4/dt**2
h.close()

text_file = open(saveFile, "w")
text_file.write('t; x0; x1; min{dx}; fx0; fy0; fz0 \n')


print 't; x0; x1; min{dx}; fx0; fy0; fz0'
for it, fname in enumerate(fnames):
    timestep=fname.split('.')[-4]
    proc = int(fname.split('.')[-2])
    parts = filter(lambda x: timestep + '.' in x, allFiles)
    print timestep, it*100.0/len(fnames)
    storageDict = {}
    for part in parts: # Check each processor output
        pid_file = int( part.split('.')[-2] )
        h = h5.File(part)
        for ic, cellId in enumerate(h['cellIds'].value):
            cId = int(cellId)
            for key in keys2get:
                ra[key][it, cId] = h[key].value[ic]
            ra['t'][it, cId] = np.int(timestep)
        h.close()
    text_file.write(str(ra['t'][it,0][0]*dt) +
                   " " + str(ra['Position'][it,0][0]*dx) +
                   " " + str(ra['Position'][it,1][0]*dx) + 
                   " " + str((ra['Max positions'][it,0] - ra['Min positions'][it,1])[0]*dx) +
                   " " + str(ra['Force'][it,1][0]*df) +
                   " " + str(ra['Force'][it,1][1]*df) +
                   " " + str(ra['Force'][it,1][2]*df) + "\n")

text_file.close()


ra['t'] *= dt
ra['Position'] *= dx
ra['Max positions'] *= dx
ra['Min positions'] *= dx
ra['Force'] *= df

Dist= ['None']*len(ra)
Time= ['None']*len(ra)
ForceZ = ['None']*len(ra)


for i in range(len(ra)):
    print ra['t'][i,0][0],
    print ";", ra['Position'][i,0][0],
    print ";", ra['Position'][i,1][0], 
    print ";", (ra['Max positions'][i,0] - ra['Min positions'][i,1])[0],
    print ";", ra['Force'][i,1][0],
    print ";", ra['Force'][i,1][1],
    print ";", ra['Force'][i,1][2]
    Dist[i]= (ra['Max positions'][i,0] - ra['Min positions'][i,1])[0]
    Time[i]= ra['t'][i,0][0]
    ForceZ[i] = ra['Force'][i,1][2]


#np.save(saveFile, ra)


plot.subplot(2, 1, 1)
pylab.plot( Time, Dist,  label="distance" )
#pylab.legend()
pylab.title("Platelet Distance in time")
pylab.xlabel("Time (s)")
pylab.ylabel("Distance (m)")

plot.subplot(2, 1, 2)
pylab.plot( Time, ForceZ,  label="forceZ" )
#pylab.legend()
pylab.title("Platelet Adhesion Force (X) in time")
pylab.xlabel("Time (s)")
pylab.ylabel("ForceX (N)")


#plot.show()
#savefig(saveFile.split[-4]+".png")






