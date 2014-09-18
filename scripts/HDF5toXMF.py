#!/usr/bin/env python
import sys
import os
from glob import glob
import errno
import numpy as np
import pylab as pl
import h5py as h5


class XMLIndentation(object):
	"""docstring for XMLIndentation"""
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
		self.currentIndent += 1
		return self.current()

	def decrease(self):
		self.currentIndent -= 1
		return self.current()



def updateDictForXDMFStrings(h5File, h5dict):
	h5dict['subDomainNx'] = h5dict["subdomainSize"][0]
	h5dict['subDomainNy'] = h5dict["subdomainSize"][1]
	h5dict['subDomainNz'] = h5dict["subdomainSize"][2]
	h5dict['relativePositionX'] = h5dict["relativePosition"][0]
	h5dict['relativePositionY'] = h5dict["relativePosition"][1]
	h5dict['relativePositionZ'] = h5dict["relativePosition"][2]
	def getAttributeTypeAndRankStringFromObjectsShape(obj_shape):
		dimObjShape = len(obj_shape)
		if dimObjShape==3:
			attributeType, rankString = "Scalar", ""
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
		"attributeName": obj.name[1:], 
		"objectSubdomainSize": obj.shape,
		"AttributeType": getAttributeTypeAndRankStringFromObjectsShape(obj.shape)[0],
		"rankString":getAttributeTypeAndRankStringFromObjectsShape(obj.shape)[1],
		}
		return d
	h5dict['DataSets'] =  map(dictFromObject,  h5File.listobjects())
	return h5dict

def iteratePossibleDataSetsDict(h5dict):
	for datasetDict in h5dict['DataSets']:
		h5dict.update(datasetDict)
		yield h5dict


def createH5TopologyAndGeometryFluid(xmlInt=XMLIndentation()):
	h5TopologyAndGeometry  = xmlInt.cur() + '<Topology TopologyType="3DCoRectMesh" NumberOfElements="%(subDomainNx)d %(subDomainNy)d %(subDomainNz)d"/>\n'
	h5TopologyAndGeometry += xmlInt.cur() + '<Geometry GeometryType="Origin_DxDyDz">\n'
	h5TopologyAndGeometry += xmlInt.inc() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
	h5TopologyAndGeometry += xmlInt.inc() + '%(relativePositionX)d %(relativePositionY)d %(relativePositionZ)d\n'
	h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
	h5TopologyAndGeometry += xmlInt.cur() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
#	h5TopologyAndGeometry += xmlInt.inc() + '%(dx)f %(dx)f %(dx)f\n'
	h5TopologyAndGeometry += xmlInt.inc() + '1 1 1'
	h5TopologyAndGeometry += xmlInt.dec() + '</DataItem>\n'
	h5TopologyAndGeometry += xmlInt.dec() + '</Geometry>\n'
	return h5TopologyAndGeometry


def createH5AttibuteFluid(xmlInt=XMLIndentation()):
	h5Attribute  = xmlInt.cur() + '<Attribute Name="%(attributeName)s" AttributeType="%(AttributeType)s" Center="Cell">\n'
	h5Attribute += xmlInt.inc() + '<DataItem Dimensions="%(subDomainNx)d %(subDomainNy)d %(subDomainNz)d %(rankString)s" Format="HDF">\n'
	h5Attribute += xmlInt.inc() + '%(pathToHDF5)s:/%(attributeName)s\n'
	h5Attribute += xmlInt.dec() + '</DataItem>\n'
	h5Attribute += xmlInt.dec() + '</Attribute>\n'
	return h5Attribute



class HDF5toXDMF_Fluid(object):
	"""docstring for HDF5toXDMF_Fluid"""
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
		stringToWrite += self.xmlInt.cur() + '<Xdmf Version="2.0">\n'
		stringToWrite += self.xmlInt.inc() + '<Domain>\n'
		self.xmlInt.inc()
		self.xdmfFile.write(stringToWrite)
		return -1

	def close(self):
		stringToWrite  = self.xmlInt.dec() + '</Domain>\n'
		stringToWrite += self.xmlInt.dec() + '</Xdmf>\n'
		self.xdmfFile.write(stringToWrite)
		self.xdmfFile.close()
		return -1

	def openGrid(self, gridName="Domain", gridType="Uniform"):
		stringToWrite = self.xmlInt.cur() + '<Grid Name="%s" GridType="%s">\n'%(gridName, gridType)
		self.xmlInt.inc()
		self.xdmfFile.write(stringToWrite)
		return -1

	def closeGrid(self):
		stringToWrite =  self.xmlInt.dec() + '</Grid>\n'
		self.xdmfFile.write(stringToWrite)
		return -1

	def writeSubDomain(self, h5dict):
		self.openGrid("Subdomain %d"%(h5dict['processorId'], ))
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


def readH5FileToDictionary(h5fname, close=True):
	h5dict = {}
	h5File = h5.File(h5fname, 'r')
	h5dict['pathToHDF5'] = h5fname
	for key, val in h5File.attrs.iteritems():
		h5dict[key] = val[0] if len(val) == 1 else val
	updateDictForXDMFStrings(h5File, h5dict)
	if close:
		h5File.close()
	return h5dict


def createXDMF(fnameString, processorStrings):
	fnameToSave =  ( (fnameString%("Collection"))[:-3] + '.xmf' ).split('/')[-1]
	xdmfFile = HDF5toXDMF_Fluid(fnameToSave)
	xdmfFile.openCollection("Domain")
	p = processorStrings[0]
	h5dict = readH5FileToDictionary(fnameString%(p,))
	xdmfFile.writeSubDomain(h5dict)
	numberOfProcessors = h5dict['numberOfProcessors']
	for p in processorStrings[1:numberOfProcessors]:
		h5dict = readH5FileToDictionary(fnameString%(p,))
		xdmfFile.writeSubDomain(h5dict)
	xdmfFile.closeCollection()
	xdmfFile.close()
	return fnameToSave




if __name__ == '__main__':
	dirname = './'
	identifier = 'Fluid'
	fluidH5files = sorted( glob(dirname + identifier + '*_*.h5') )
	print fluidH5files
	fluidIDs = map(lambda x: x[:-3], fluidH5files)
	iterationStrings, processorStrings  = zip(*map(lambda f: f.split('_')[-2:], fluidIDs))
	iterationStrings, processorStrings = map(lambda l: sorted(set(l)), (iterationStrings, processorStrings))
	print iterationStrings, processorStrings

	for iterString in iterationStrings:
		fnameString = dirname + identifier + "_" + iterString + "_%s.h5"
		print fnameString
		print createXDMF(fnameString, processorStrings)
		


