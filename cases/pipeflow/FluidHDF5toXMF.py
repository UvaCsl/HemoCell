#!/usr/bin/env python
import sys
import os
from glob import glob
import errno
import numpy as np
import h5py as h5
from multiprocessing import Pool


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
  h5TopologyAndGeometry = ""
  h5TopologyAndGeometry += xmlInt.inc() + '<Geometry GeometryType="ORIGIN_DXDYDZ">\n'
  h5TopologyAndGeometry += xmlInt.inc() + '<DataItem DataType="Float" Dimensions="3" Format="XML">\n'
  h5TopologyAndGeometry += xmlInt.cur() + '%(rPX)d %(rPY)d %(rPZ)d\n'
  h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
  h5TopologyAndGeometry += xmlInt.inc() + '<DataItem DataType="Float" Dimensions="3" Format="XML">\n'
  h5TopologyAndGeometry += xmlInt.cur() + '1 1 1\n'
  h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
  h5TopologyAndGeometry += xmlInt.dec() + '</Geometry>\n'
  h5TopologyAndGeometry += xmlInt.cur() + '<Topology TopologyType="3DCoRectMesh" Dimensions="%(sdSX)d %(sdSY)d %(sdSZ)d" />\n'
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
  h5Attribute += xmlInt.inc() + '<DataItem Dimensions="%(sdSX)d %(sdSY)d %(sdSZ)d %(rankString)s" Format="HDF">\n'
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
  h5dict['sdSX'] = h5dict["subdomainSize"][0]
  h5dict['sdSY'] = h5dict["subdomainSize"][1]
  h5dict['sdSZ'] = h5dict["subdomainSize"][2]
  h5dict['rPX'] = h5dict["relativePosition"][0]
  h5dict['rPY'] = h5dict["relativePosition"][1]
  h5dict['rPZ'] = h5dict["relativePosition"][2]
  def getAttributeTypeAndRankStringFromObjectsShape(obj_shape):
    dimObjShape = len(obj_shape)
    if dimObjShape==1:
      attributeType, rankString = "Scalar", ""
    elif dimObjShape==2 and obj_shape[-1] == 3:
      attributeType, rankString = "Vector", "3"
    else:
      attributeType, rankString = "Vector", " ".join(str(obj_shape[3]))
    return attributeType, rankString 
  def dictFromObject(obj):
    d = {
    "attributeName": obj[0], 
    "objectSubdomainSize": obj[1],
    "AttributeType": getAttributeTypeAndRankStringFromObjectsShape(obj[1].shape)[0],
    "rankString":getAttributeTypeAndRankStringFromObjectsShape(obj[1].shape)[1],
    }
    return d
  h5dict['DataSets'] =  map(dictFromObject,  h5File.items())
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
  h5dict['pathToHDF5'] = h5fname
  for key, val in h5File.attrs.iteritems():
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
      if not datasetDict['attributeName'] == "Triangles":
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



def createXDMF(tuble):
  """
  Reads file fitting the patters and saves an XMF file to the directory of the HDF5s.
  
    fnameString = './Fluid_00000000_%s.h5'
    processorStrings = ['p000', 'p001', 'p002']
    createXDMF(fnameString, processorStrings)

    Important: reads as many files as number of processors in the first file
  """
  fnameString, processorStrings = tuble
  fnameToSave =  ( (fnameString%(""))[:-6] + '.xmf' )
  fnameToSave = fnameToSave.replace('/hdf5/','/').replace(".p.",".").replace(".pgz.",".")
  if os.path.isfile(fnameToSave):
            return fnameToSave + " (existed)"
  xdmfFile = HDF5toXDMF_Fluid(fnameToSave)
  xdmfFile.openCollection("Domain")
  processorStrings = sorted(processorStrings)
  p = processorStrings[0]
  if not os.path.isfile(fnameString%(p,)) : fnameString = fnameString.replace('.p.','.pgz.')
  h5dict = readH5FileToDictionary(fnameString%(p,))
  numberOfProcessors = len(processorStrings)
  for p in processorStrings[0:numberOfProcessors]:
    h5dict = readH5FileToDictionary(fnameString%(p,))
    xdmfFile.writeSubDomain(h5dict, "Subdomain " + p)
  xdmfFile.closeCollection()
  xdmfFile.close()
  print "Created file:", fnameToSave

if __name__ == '__main__':
  dirname = './hdf5/'
  identifier = "Fluid"
  fluidH5files = sorted( glob(dirname + identifier + '*.h5') )

  fluidIDs = map(lambda x: x[:-3], fluidH5files)
  iterationStrings, processorStrings  = zip(*map(lambda f: [f.split('.')[-3], f.split('.')[-1]], fluidIDs))
  iterationStrings, processorStrings = map(lambda l: sorted(set(l)), (iterationStrings, processorStrings))
  fnameStrings = []
  for iterString in iterationStrings:
    fnameStrings.append(dirname + identifier + "." + iterString + ".p.%s.h5")
    pool = Pool(16)
    pool.map(createXDMF,[(fnameString,processorStrings) for fnameString in fnameStrings])
    pool.close()


