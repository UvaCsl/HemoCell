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
saveFile =  options.outputfile
Np = int(options.Np)


allFiles = sorted(glob(caseFolder + '/*BondFieldParticles*.h5'))
fnames = filter(lambda x: '.h5' in x, allFiles)
Nt = len(fnames)


keys2get = ['numberOfParticles']
dtype = np.dtype([('numberOfParticles', '<i4', (1,)),  ('t', '<f4', (1,))  ])

ra = np.zeros((Nt), dtype=dtype)

h = h5.File(fnames[0])
dt = h.attrs['dt'][0]
dx = h.attrs['dx'][0]
df = 1025 * dx**4/dt**2
h.close()

text_file = open(saveFile, "w")
text_file.write('t; BndNum \n')


print 't; BndNum '
for it, fname in enumerate(fnames):
    timestep=fname.split('.')[-4]
    proc = int(fname.split('.')[-2])
    parts = filter(lambda x: timestep + '.' in x, allFiles)
   # print timestep, it*100.0/len(fnames)
    storageDict = {}
    for part in parts: # Check each processor output
        pid_file = int( part.split('.')[-2] )
        h = h5.File(part)
        ra['numberOfParticles'][it] += h.attrs['numberOfParticles']
        ra['t'][it] = np.int(timestep)
        h.close()
    print ra['t'][it][0]*dt, ra['numberOfParticles'][it][0]
    text_file.write(str(ra['t'][it][0]*dt) +
                   " " + str(ra['numberOfParticles'][it][0])  + "\n")

text_file.close()


ra['t'] *= dt


BndNum= ['None']*len(ra)
Time= ['None']*len(ra)


for i in range(len(ra)):
    print ra['t'][i][0],
    print ";", ra['numberOfParticles'][i][0]
    Time[i]= ra['t'][i][0]
    BndNum[i] = ra['numberOfParticles'][i][0]


#np.save(saveFile, ra)


plot.subplot(2, 1, 1)
pylab.plot( Time, BndNum,  label="bndnum" )
#pylab.legend()
pylab.title("Number of bonds in time")
pylab.xlabel("Time (s)")
pylab.ylabel("Number of bonds")




#plot.show()
#savefig(saveFile.split[-4]+".png")






