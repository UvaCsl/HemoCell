#!/usr/bin/env python
import numpy as np
import pylab as pl
import h5py as h5
from glob import glob
import sys


def writeProgressBar(string2write):
  """ Quick and dirty progress bar. 
  If string2write is an integer or float, writeProgressBar takes care of the
  format (has to be [0,100]).
  """
  if type(string2write) == type(1) or type(string2write) == type(1.2):
    string2write = "|"+  "-"*int(string2write/10)  + "| %.02f"%(string2write)
  sys.stdout.write("\r")
  sys.stdout.write(string2write)
  sys.stdout.flush()
  return

def dtypeFromH5(h5file):
    """ Extracts a numpy dtype from an h5 file.

    Keyword arguments:
        h5file -- opened h5 file.

    Returns:
        numpy dtype file.

    """
    dtypeNames, dtypeSize = [], []
    
    for i,v in h5file.iteritems():
        if ' ' not in i:
            dtypeNames += [str(i)]
            dtypeSize += [ str(v.shape[-1]) + '<' + str(v.dtype)[0] + '4' ]
    for i,v in h5file.attrs.iteritems():
        if ' ' not in i:
            dtypeNames += [str(i)]
            dtypeSize += [ str(v.shape[-1]) + '<' + str(v.dtype)[0] + '4' ]

    dtype = np.dtype({'names':dtypeNames, 'formats':dtypeSize} )
    return dtype

def readFicsionCellResults(hdf5_dir, cellType='RBC', progressBar=False):
    """ Reads all the CellType h5 files inside a dir. 

    Keyword arguments:
        hdf5_dir -- string containing the output of a ficsion run.
        cellType -- The CellType to convert to a numpy type (eg. RBC, PLT, WBC). 
            Reads {cellType}_Cell3D.* files.

    Returns:
        numpy array containing all the data from all the timesteps.
            Dimensions are [timestep, cell, dimension]


    """

    dtype = np.dtype([('Position', '<f4', (3,)), ('Velocity', '<f4', (3,)), ('cellIds', '<f4', (1,)), 
                      ('Surface', '<f4', (1,)),  ('Volume', '<f4', (1,)), ('t', '<f4', (1,))])
    keys2get = list(dtype.names)
    keys2get.remove('t')

    identifier = hdf5_dir + '/'+cellType+'_Cell3D.*'
    hdf5_files = glob(identifier)

    # Extract the sorted list of timesteps in strings. Sorted strings can be different than integers.
    timesteps = sorted(set(map(lambda x: x.split('.')[-4], hdf5_files)))
    timesteps_int = map(int, timesteps)
    timesteps = zip( *sorted(zip(timesteps_int,timesteps)) )[1]


    Nt = len(timesteps)
    ra = None; 
    for it, timestep in enumerate(timesteps): # Loop through all the timesteps and processors.
        if progressBar: writeProgressBar( it*100.0/Nt )
        parts = filter(lambda x: timestep + '.' in x, hdf5_files)
        storageDict = {}
        for part in parts: # Collect output from each processor
            pid_file = int( part.split('.')[-2] )
            h = h5.File(part)
            for ic, cellId in enumerate(h['cellIds'].value):
                cId = int(cellId)
                sD = storageDict[cId] = storageDict.get(cId, {})
                for key in keys2get:
                    sD[key] = h[key].value[ic]
                sD['t'] = h.attrs['dt'][0] * np.int(timestep)
            h.close()
        if type(ra) == type(None): # Initialize array if not initialized.
            Np = len(storageDict.keys())
            ra = np.empty((Nt, Np), dtype=dtype)
            # There are cases where certain cellIds are missing.
            cellId2Id = dict( zip( storageDict.keys(),np.arange(Np) ) )
        for cellId, sD in storageDict.items():
            for key in sD.keys():
                ra[key][it, cellId2Id[cellId] ] = sD[key]
    if progressBar: 
        writeProgressBar( 100.0 )
        print 
    return ra



# ############################ #
#       Helper functions       #
# ############################ #

def rotateVector(vector, thetaDegrees):
    """ vector has dimensions [Nt, Np, 3], when Nt are the timesteps, Np are the number of particles and 3 the dimensions. """
    thetaRadians = thetaDegrees * np.pi/180.0
    newVector = vector.copy()
    newVector[:, :, 0] = vector[:,:, 0] * np.cos(thetaRadians) - vector[:,:, 1] * np.sin(thetaRadians)
    newVector[:, :, 1] = vector[:,:, 0] * np.sin(thetaRadians) + vector[:,:, 1] * np.cos(thetaRadians)
    return newVector



def main():
    import argparse, os

    parser = argparse.ArgumentParser(description='Save CellType data in a numpy-array file.')
    parser.add_argument('-i', '--inputDir', dest='inputDir', default='./tmp/hdf5',
                       help='Input directory where HDF5 files are located [default:./tmp/hdf5]')
    parser.add_argument('-o', '--output', dest='outputFile', default='./ficsion_hdf5.npy',
                       help='Path to save the numpy file [default:./ficsion_hdf5.npy].')
    parser.add_argument('-t', '--cellType', dest='cellType', default='RBC',
                       help='CellType [default:RBC]')

    args = parser.parse_args()

    if not os.path.exists(args.inputDir):
        print 'ERROR: \n\t "', args.inputDir, '" does not exist.'
        sys.exit(1)
    numpy_array = readFicsionCellResults(args.inputDir, args.cellType, True)
    np.save(args.outputFile, numpy_array)
    return

if __name__ == '__main__':
    main()



