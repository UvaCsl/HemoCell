#!/usr/bin/env python 
# encoding: utf-8
#====================================================================
# Name        :  
# Version    : 
# Author      :  Lampros Mountrakis (L.Mountrakis@uva.nl)
# Date        : 
# Description : 
#
#
#====================================================================

import pylab as pl
import numpy as np
import sys
from glob import glob


def readStretchingData(fname):
    d = np.loadtxt(fname, delimiter=";",usecols = (0,1,2,3))
    t, f, DAX, DTX = [d[:,i] for i in range(4)]
    f *= 1e12
    DAX *= 1e6
    DTX *= 1e6
    return t, f, DAX, DTX


def readValidationData(dataDir="./validation/"):
    DA    = np.loadtxt(dataDir + "Stretching-FedosovCaswell2010b.DA.dat")
    DT    = np.loadtxt(dataDir + "Stretching-FedosovCaswell2010b.DT.dat")
    DAstd = np.loadtxt(dataDir + "Stretching-FedosovCaswell2010b.DA-STD.dat")
    DTstd = np.loadtxt(dataDir + "Stretching-FedosovCaswell2010b.DT-STD.dat")
    DAstd = DAstd[:,1] - DA[:,1]
    DTstd = DTstd[:,1] - DT[:,1]
    return DA, DT, DAstd, DTstd

def plotValidationData(dataDir="./validation/"):
    DA, DT, DAstd, DTstd = readValidationData(dataDir)
    pl.errorbar(DA[:,0], DA[:,1], yerr=DAstd, fmt='ro', label="DA")
    pl.errorbar(DT[:,0], DT[:,1], yerr=DTstd, fmt='bo', label="DT")
    pl.xlabel("Force [pN]")
    pl.ylabel("Deformation [um]")
    pass


def plotMultipleStretchingData(filenames, var='kBend', val='1000', var2='kWLC'):
    fchars = map(lambda z: dict(map(lambda x:x.split('-'), z.split('_'))[:-1]), filenames)
    clrs = 'rgbcmyk'
    for ifn, fn in enumerate(filenames):
        t, f, DAX, DTX = readStretchingData(fn)
        kVar2 = fchars[ifn][var2]
        pl.plot(f, DAX, clrs[ifn], label=var2 +" = %s"%(kVar2))
        pl.plot(f, DTX, clrs[ifn])
    pl.title(var + '-' + val)
    pl.axis([0, 200, 0, 25])
    pl.legend(loc="best")
    pl.savefig(var + '-' + val + '.png')
    pass



if 0:
    kWLC = {'0.1':0,'0.25':1,'0.5':2,'1.0':3}
    kBend = {'75':0,'150':1,'200':2,'1000':3}
    var='kBend'; val='1000'; var2='kWLC'
    filenames = glob('stretch00/results/rbcModel-1*' + var + '-' + val + '_*dat')
    plotValidationData(dataDir="./validation/")
    plotMultipleStretchingData(filenames, var, val, var2)
    pass


# pl.show()
if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except:
        print "Error: No or wrong filename as input."
    plotValidationData(dataDir="./validation/")
    t, f, DAX, DTX = readStretchingData(fname)
    pl.plot(f, DAX, label="DA sim")
    pl.plot(f, DTX, label="DT sim")
    pl.axis([0, 200, 0, 25])
    pl.legend(loc="best")
    pl.savefig(fname + '.png')
    print "Saving plot as ", fname + '.png'




