#!/usr/bin/env python3
import sys
import os
from glob import glob
import errno
import numpy as np
import h5py as h5


class XMLIndentation(object):
    """
    Elementary class for keeping track of the indentation in an XML file (and not only).

    xmlInt = XMLIndentation(indent=2)
    
    Parameters:
    -----------
        indent: int, number of spaces for each indent.

    Examples:
    -----------
        stringToWrite = xmlInt.increase()  + "<parameters>"
        stringToWrite = xmlInt.current() + "<flowType> 0 </flowType>"
        stringToWrite = xmlInt.current() + "<Re> 0 </Re>"
        stringToWrite = xmlInt.decrease() + "</parameters>"

    """
    def __init__(self, indent=2):
        super(XMLIndentation, self).__init__()
        self.indentSize = indent
        self.currentIndent = 0
        self.current = lambda : " "*(self.indentSize * self.currentIndent)
        # Shortcuts
        self.inc = lambda : self.increase()
        self.dec = lambda : self.decrease()
        self.cur = lambda : self.current()
        pass

    def increase(self):
        ret = self.current()
        self.currentIndent += 1
        return ret

    def decrease(self):
        self.currentIndent -= 1
        return self.current()



def createH5TopologyAndGeometryFluid(xmlInt=XMLIndentation()):
    """
        String to be used for Topology and geometry tags of the XMF file related to the Fluid field.
    Output:
    ---------
        Returns string with substitutable field.
        Fields are:
            subDomainNx, subDomainNy, subDomainNz
            relativePositionX, relativePositionY, relativePositionZ
            dx (not used)
    """
    h5TopologyAndGeometry  = xmlInt.cur() + '<Topology TopologyType="3DCoRectMesh" NumberOfElements="%(subDomainNx)d %(subDomainNy)d %(subDomainNz)d"/>\n'
    h5TopologyAndGeometry += xmlInt.inc() + '<Geometry GeometryType="Origin_DxDyDz">\n'
    h5TopologyAndGeometry += xmlInt.inc() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
    h5TopologyAndGeometry += xmlInt.cur() + '%(relativePositionX)G %(relativePositionY)G %(relativePositionZ)G\n'
    h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
    h5TopologyAndGeometry += xmlInt.inc() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
# Use for SI units #    h5TopologyAndGeometry += xmlInt.inc() + '%(dx)f %(dx)f %(dx)f\n'
    h5TopologyAndGeometry += xmlInt.cur() + '%(dX)G %(dY)G %(dZ)G\n'
    h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
    h5TopologyAndGeometry += xmlInt.dec() + '</Geometry>\n'
    return h5TopologyAndGeometry


def createH5AttibuteFluid(xmlInt=XMLIndentation()):
    """
        String to be used for the Attribute field of the XMF file related to the Fluid field.
    Output:
    ---------
        Returns string with substitutable field.
        Fields are:
            attributeName, AttributeType
            subDomainNx, subDomainNy, subDomainNz
    """
    h5Attribute  = xmlInt.inc() + '<Attribute Name="%(attributeName)s" AttributeType="%(AttributeType)s" Center="Cell">\n'
    h5Attribute += xmlInt.inc() + '<DataItem Dimensions="%(subDomainNx)d %(subDomainNy)d %(subDomainNz)d %(rankString)s" Format="HDF">\n'
    h5Attribute += xmlInt.cur() + '%(pathToHDF5)s:/%(attributeName)s\n'
    h5Attribute += xmlInt.dec() + '</DataItem>\n'
    h5Attribute += xmlInt.dec() + '</Attribute>\n'
    return h5Attribute


def updateDictForXDMFStrings(h5File, h5dict):
    """
        Updates the h5dict in order to be used with the above strings.
        Works in combination with:
            readH5FileToDictionary(h5fname, close=True)
                and
            iteratePossibleDataSetsDict(h5dict)

    Usage:
    --------
        h5dict = updateDictForXDMFStrings(h5File, h5dict)

    Input
    ---------
        h5File: HDF5 file
        h5dict: Dictionary with the meta-information of the HDF5 file (no data contained)

    Output:
    ---------
        h5dict: Dictionary with the meta-information of the HDF5 file (no data contained)
    """
    h5dict['subDomainNx'] = h5dict["subdomainSize"][0]
    h5dict['subDomainNy'] = h5dict["subdomainSize"][1]
    h5dict['subDomainNz'] = h5dict["subdomainSize"][2]
    h5dict['dX'] = h5dict["dxdydz"][0]
    h5dict['dY'] = h5dict["dxdydz"][1]
    h5dict['dZ'] = h5dict["dxdydz"][2]
    h5dict['relativePositionX'] = h5dict["relativePosition"][0]
    h5dict['relativePositionY'] = h5dict["relativePosition"][1]
    h5dict['relativePositionZ'] = h5dict["relativePosition"][2]
    def getAttributeTypeAndRankStringFromObjectsShape(obj_shape):
        dimObjShape = len(obj_shape)
        if dimObjShape==3:
            attributeType, rankString = "Scalar", ""
        elif dimObjShape==4 and obj_shape[-1] == 1:
            attributeType, rankString = "Scalar", "1"
        elif dimObjShape==4 and obj_shape[-1] == 3:
            attributeType, rankString = "Vector", "3"
        elif dimObjShape==4 and obj_shape[-1] == 6:
            attributeType, rankString = "Tensor6", "6"
        elif dimObjShape==4 and obj_shape[-1] == 9:
            attributeType, rankString = "Tensor", "9"
        else:
            attributeType, rankString = "Matrix", " ".join(obj_shape[3:])
        return attributeType, rankString 
    def dictFromObject(obj):
        d = {
        "attributeName": obj[0], 
        "objectSubdomainSize": obj[1],
        "AttributeType": getAttributeTypeAndRankStringFromObjectsShape(obj[1].shape)[0],
        "rankString":getAttributeTypeAndRankStringFromObjectsShape(obj[1].shape)[1],
        }
        return d
    h5dict['DataSets'] =  list(map(dictFromObject,  list(h5File.items())))
    return h5dict


def readH5FileToDictionary(h5fname, close=True):
    """
    Reads the meta-information of an HDF5 file to a dictionary

    Usage:
    --------
        h5dict = readH5FileToDictionary(h5File, close=True)

    Input
    ---------
        h5File: HDF5 file
        close: boolean [optional], to close the HDF5 file or not. If not, the file is returned.

    Output:
    ---------
        h5dict: Dictionary with the meta-information of the HDF5 file (no data contained)
        h5File: [if close == False]. The HDF5 file
    """

    h5dict = {}
    h5File = h5.File(h5fname, 'r')
    h5dict['pathToHDF5'] = h5fname.replace('//','/')
    for key, val in list(h5File.attrs.items()):
        h5dict[key] = val[0] if len(val) == 1 else val
    updateDictForXDMFStrings(h5File, h5dict)
    if close:
        h5File.close()
        ret = h5dict
    else:
        ret = h5dict, h5File
    return ret      


def iteratePossibleDataSetsDict(h5dict):
    """
        Generator updating each time the h5dict with the values of the dictionary
            of h5dict['DataSets']. 
        Used to iterate through meta-information of Datasets 
            (velocity, density etc to Vector, Scalar etc)
    """
    for datasetDict in h5dict['DataSets']:
        h5dict.update(datasetDict)
        yield h5dict



class HDF5toXDMF_Fluid(object):
    """
        Class converting constructing an XMF file from HDF5 information (h5dict values)

            Makes use of the functions 
                    createH5TopologyAndGeometryFluid and createH5AttibuteFluid

    Examples
    ---------
        # Just from one file:
            xmfFile = HDF5toXDMF_Fluid('./test.xmf')
            xmfFile.writeSubDomain(h5dict, gridName="Plasma")
            xmfFile.close()
        # Just from one file:
            xmfFile = HDF5toXDMF_Fluid('./test.xmf')
            xmfFile.openCollection("Domain")
            xmfFile.writeSubDomain(h5dictProc0)
            xmfFile.writeSubDomain(h5dictProc1)
            xmfFile.closeCollection()
            xmfFile.close()
    """
    def __init__(self, fname=None):
        super(HDF5toXDMF_Fluid, self).__init__()
        self.xmlInt = XMLIndentation()
        if fname != None:
            self.open(fname)
        pass

    def open(self, fname):
        self.xdmfFile  = open(fname, 'w')
        stringToWrite  = self.xmlInt.cur() + '<?xml version="1.0" ?>\n'
        stringToWrite += self.xmlInt.cur() + '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
        stringToWrite += self.xmlInt.inc() + '<Xdmf Version="2.0">\n'
        stringToWrite += self.xmlInt.inc() + '<Domain>\n'
        self.xdmfFile.write(stringToWrite)
        return -1

    def close(self):
        stringToWrite  = self.xmlInt.dec() + '</Domain>\n'
        stringToWrite += self.xmlInt.dec() + '</Xdmf>\n'
        self.xdmfFile.write(stringToWrite)
        self.xdmfFile.close()
        return -1

    def openGrid(self, gridName="Domain", gridType="Uniform"):
        stringToWrite = self.xmlInt.inc() + '<Grid Name="%s" GridType="%s">\n'%(gridName, gridType)
        self.xdmfFile.write(stringToWrite)
        return -1

    def closeGrid(self):
        stringToWrite =  self.xmlInt.dec() + '</Grid>\n'
        self.xdmfFile.write(stringToWrite)
        return -1

    def writeSubDomain(self, h5dict, gridName="Subdomain "):
        pId = (" " + str(h5dict['processorId'])) if 'processorId' in h5dict else ""
        self.openGrid(gridName + pId)
        stringToWrite = createH5TopologyAndGeometryFluid(self.xmlInt)%(h5dict)
        for datasetDict in iteratePossibleDataSetsDict(h5dict):
            stringToWrite += createH5AttibuteFluid(self.xmlInt)%(datasetDict)
        self.xdmfFile.write(stringToWrite)
        self.closeGrid()
        return -1

    def openCollection(self, gridName="Domain"):
        self.openGrid(gridName=gridName, gridType="Collection")
        return -1

    def closeCollection(self):
        self.closeGrid()
        return -1



def createXDMF(fnameString, processorStrings, iterDir):
    """
    Reads file fitting the patters and saves an XMF file to the directory of the HDF5s.
    
        fnameString = './Fluid_00000000_%s.h5'
        processorStrings = ['p000', 'p001', 'p002']
        createXDMF(fnameString, processorStrings)

        Important: reads as many files as number of processors in the first file
    """
    fnameToSave =  ( (fnameString%(""))[:-6] + '.xmf' )
    fnameToSave = fnameToSave.replace('/hdf5/','/').replace(".p.",".").replace(".pgz.",".").replace("/" + iterDir + "/", "")
    if os.path.isfile(fnameToSave):
        return fnameToSave + " (existed)"
    xdmfFile = HDF5toXDMF_Fluid(fnameToSave)
    xdmfFile.openCollection("Domain")
    processorStrings = sorted(processorStrings)
    p = processorStrings[0]
    if not os.path.isfile(fnameString%(p,)) : fnameString = fnameString.replace('.p.','.pgz.')
    h5dict = readH5FileToDictionary(fnameString%(p,))
    #xdmfFile.writeSubDomain(h5dict)
    numberOfProcessors = len(processorStrings)
    for p in processorStrings[0:numberOfProcessors]:
        h5dict = readH5FileToDictionary(fnameString%(p,))
        xdmfFile.writeSubDomain(h5dict, "Subdomain " + p)
    xdmfFile.closeCollection()
    xdmfFile.close()
    print("Created file:", fnameToSave)


if __name__ == '__main__':
    try:
        if os.environ["ABSPATH"] == "1":
            dirname = os.path.abspath('./hdf5/') + '/'
        else:
            dirname = './hdf5/'
    except:
        dirname = './hdf5/'
    
    directories = sorted(os.listdir(dirname))
    identifier = 'Fluid'
    
    for iterDir in directories:
        fluidH5files = sorted( glob(dirname + '/' + iterDir + '/' + identifier + '*p*.h5') )
        fluidIDs = [x[:-3] for x in fluidH5files]
        iterationStrings, processorStrings  = list(zip(*[[f.split('.')[-3], f.split('.')[-1]] for f in fluidIDs]))
        iterationStrings, processorStrings = [sorted(set(l)) for l in (iterationStrings, processorStrings)]
        for iterString in iterationStrings:
            fnameString = dirname + '/' + iterDir + '/' + identifier + "." + iterString + ".p.%s.h5"
            createXDMF(fnameString, processorStrings, iterDir)



