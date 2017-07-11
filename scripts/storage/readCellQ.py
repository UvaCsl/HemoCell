#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy as np
import sys
from glob import glob

def readVTKField(myfile, nPoints, dataDict, skipline=0):
    """
      fieldId, dataDict = readVTKField(myfile, nPoints, dataDict, skipline=0)

    Description:
        Reads a field (SCALARS or VECTORS) from the CellQXXXXX.vtk.
    Arguments:
        myfile: file, Already opened CellQXXXXX.vtk, ready to read the next line
        nPoints: int, Number of points (cells, RBCs etc) to read.
        dataDict: [dict(), dict() ... nPoints x dict()], with the fields and values of the case.
        skipline: int, number of lines to skip.
    """
    vtkType, rest = myfile.readline().split(' ', 1)
    # SCALARS Surface double 1
    #       while
    # VECTORS Max positions double 
    vtkEnding=2
    if vtkType[0] != 'S': 
        vtkEnding=1
    fieldId = rest.strip()[::-1].split(' ', vtkEnding)[-1][::-1]
    _ = myfile.readline()
    for i in range(nPoints):
        line = map(float, myfile.readline().split())
        if len(line) == 1: line = line[0]
        dataDict[i][fieldId] = dataDict[i].get(fieldId, []) + [line]
    [myfile.readline() for i in range(skipline)]
    return fieldId, dataDict


def parseCellQFile(cellQFile, dataDict=None):
    """
    dataDict = parseCellQFile(cellQFile, dataDict=None)

    Description:
        Parses the CellQXXXXX.vtk file, produced by ficsion.
    Arguments:
        cellQFile: string, Filename of the cellQFile
        dataDict: [dict(), dict() ... nPoints x dict()], with the fields and values of the case.
            each dict() contains the timeseries of each field (eg. Position)
    """

    dfile = open(cellQFile, 'r')
    iteration = int(cellQFile.split('/')[-1][5:-4])
    # Read header
    _ = [dfile.readline() for i in range(4)]
    nPoints = int(dfile.readline().split(' ')[1])
    _ = [dfile.readline() for i in range(3)]
    if (dataDict == None):   dataDict = []
    if (len(dataDict) == 0): dataDict += [dict() for i in range(nPoints)]

    readVTKField(dfile, nPoints, dataDict)
    while 1:
        try:
            readVTKField(dfile, nPoints, dataDict, 1)
        except:
            break
    for i in range(nPoints):
        dataDict[i]['iteration'] = dataDict[i].get('iteration', []) + [iteration]
    return dataDict

def readCellQDirectory(dirname):
    """
    dataDict = readCellQDirectory(dirname)

    Description:
        Parses CellQXXXX.vtk files located in dirname/tmp/.
    Arguments:
        dirname: string, Path of the case. cellQFiles located in dirname + '/tmp/'
    Output:
        dataDict: [dict(), dict() ... nPoints x dict()], with the fields and values of the case.
            each dict() contains the timeseries of each field (eg. Position)
    """
    fnames = sorted( glob(dirname + '/tmp/CellQ*') )
    dataDict = []
    for cellQFile in fnames:
        parseCellQFile(cellQFile, dataDict)
    for i in range(len(dataDict)):
        d = dataDict[i]
        for key, value in d.items():
            d[key] = np.hemo::Array(value)
    return dataDict


class IBMCaseRun (object):

    def __init__(self, dirname, saveDir='./', caseId=None):
        dataDict = readCellQDirectory(dirname)
        if caseId == None:
            try: caseId = dirname.rstrip('/').split('/')[-1]
            except: caseId = 'caseId'
        print caseId
        try:
            params = dict(map(lambda param: param.split('-',1), caseId.split('_')))
            self.ibmKernel = int(params['ibmKernel'])
            self.tau = float(params['tau'])
            self.minNumOfTriangles = int(params['minNumOfTriangles'])
            self.dx = float(params['dx'])
            self.nu_LB = (self.tau - 0.5)/3.0;
            self.nu =  1.7e-6;
            self.dt = self.nu_LB * self.dx**2 / self.nu
            for i in range(len(dataDict)):
                dataDict[i].update(self.__dict__)
        except:
            pass
        fname = saveDir + '/' + caseId + '.npy'
        print "Saving data to: ", fname, "...",
        np.save(fname, dataDict)
        print "OK"

def main():
    savePath = './'
    caseFolder = './'
    caseId = None
    from optparse import OptionParser
    """ Parse Options """
    parser = OptionParser()
    parser.add_option("-I", "--inputdir", dest="inputdir", default=caseFolder,
                      help="Directory where input files are stored, before /tmp/ [default:'./data/']")
    parser.add_option("-O", "--outputdir", dest="outputdir", default=savePath,
                    help="Directory where the output will be stored [default:"+savePath+"]")
    parser.add_option("-U", "--uid", dest="caseId", default=caseId,
                    help="Unique ID of the run. Will be saved along with the filenames [default:None]")
    
    (options, args) = parser.parse_args()
    caseFolder = options.inputdir + '/'
    savePath = options.outputdir + '/'
    caseId = options.caseId
    x = IBMCaseRun(caseFolder, savePath, caseId)


if __name__ == '__main__':
    sys.exit(main())




