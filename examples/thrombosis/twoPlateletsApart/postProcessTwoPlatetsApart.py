#!/usr/bin/env python
import numpy as np
import pylab as pl
import h5py as h5
from glob import glob
import sys 

from optparse import OptionParser
""" Parse Options """
parser = OptionParser() 
parser.add_option("-I", "--inputdir", dest="inputdir", default='./',
                  help="Directory where input files are stored (without './tmp/hdf5/')[default:'./']")            
parser.add_option("-O", "--outputfile", dest="outputfile", default='./test.npy',
                help="Directory where the output will be stored [default: './test.npy' ]")
parser.add_option("-N", "--numberOfParticles", dest="Np", default=2,
                help="Directory where the output will be stored [default: '2' ]")


(options, args) = parser.parse_args()
caseFolder = options.inputdir + '/tmp/hdf5/'
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


print 't; x0; x1; min{dx}; fx0; fy0; fz0'
for it, fname in enumerate(fnames):
    timestep=fname.split('.')[-4]
    proc = int(fname.split('.')[-2])
    parts = filter(lambda x: timestep + '.' in x, allFiles)
#    print timestep, it*100.0/len(fnames)
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


ra['t'] *= dt
ra['Position'] *= dx
ra['Max positions'] *= dx
ra['Min positions'] *= dx
ra['Force'] *= df

for i in range(len(ra)):
    print ra['t'][i,0][0],
    print ";", ra['Position'][i,0][0],
    print ";", ra['Position'][i,1][0], 
    print ";", (ra['Max positions'][i,0] - ra['Min positions'][i,1])[0],
    print ";", ra['Force'][i,1][0],
    print ";", ra['Force'][i,1][1],
    print ";", ra['Force'][i,1][2]

# np.save(saveFile, ra)



