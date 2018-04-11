#!/usr/bin/env python3
import h5py as h5
import os
from glob import glob

class XMLIndentation(object):
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

class createXDMF(object):

  def readH5File(self,h5fname):
    h5dict = {}
    h5File = h5.File(h5fname, 'r')
    h5dict['pathToHDF5'] = h5fname.replace('//','/')
    for key, val in list(h5File.attrs.items()):
      h5dict[key] = val[0] if len(val) == 1 else val
    h5dict["datasets"] = []
    for ds, val in h5File.items():
      dimObjShape = len(val.shape)
      if dimObjShape==3:
        attributeType, rankString = "Scalar", ""
      elif dimObjShape==4 and val.shape[-1] == 1:
        attributeType, rankString = "Scalar", "1"
      elif dimObjShape==4 and val.shape[-1] == 3:
        attributeType, rankString = "Vector", "3"
      elif dimObjShape==4 and val.shape[-1] == 6:
        attributeType, rankString = "Tensor6", "6"
      elif dimObjShape==4 and val.shape[-1] == 9:
        attributeType, rankString = "Tensor", "9"
      else:
        attributeType, rankString = "Matrix", " ".join(val.shape[3:])
      d = {
      "attributeName": ds, 
      "AttributeType": attributeType ,
      "rankString":rankString
      }
      h5dict["datasets"].append(d)
      
      
    h5File.close()
    return h5dict
  def writeHeader(self):
    self.output += self.xmlInt.cur() + '<?xml version="1.0" ?>\n'
    #self.output += self.xmlInt.cur() + '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
    self.output += self.xmlInt.inc() + '<Xdmf Version="3.0">\n'
    self.output += self.xmlInt.inc() + '<Domain>\n'
  def writeFooter(self):
    self.output += self.xmlInt.dec() + '</Domain>\n'
    self.output += self.xmlInt.dec() + '</Xdmf>\n'
  def write_cached_subgrid(self,p):
    h5dict = self.readH5File(self.fnameString%(p,))
    if (not h5dict["processorId"] in self.cached_subgrid_output) :
      self.cached_subgrid_output[h5dict["processorId"]] = ""
    output = ""
    output += self.xmlInt.inc() + '<Grid Name="Subdomain %s %s" GridType="Uniform">\n'%(p,h5dict["processorId"])
    output += self.xmlInt.cur() + '<Topology TopologyType="3DCoRectMesh" NumberOfElements="%d %d %d"/>\n' % \
                                  tuple(h5dict["subdomainSize"])
    output += self.xmlInt.inc() + '<Geometry GeometryType="Origin_DxDyDz">\n'
    output += self.xmlInt.inc() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
    output += self.xmlInt.cur() + '%G %G %G\n'%tuple(h5dict["relativePosition"])
    output += self.xmlInt.dec() + '</DataItem>\n'
    output += self.xmlInt.inc() + '<DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n'
    output += self.xmlInt.cur() + '%G %G %G\n'%tuple(h5dict["dxdydz"])
    output += self.xmlInt.dec() + '</DataItem>\n'
    output += self.xmlInt.dec() + '</Geometry>\n'

    for ds in h5dict["datasets"]:
      output += self.xmlInt.inc() + '<Attribute Name="%s" AttributeType="%s" Center="Cell">\n'%(ds["attributeName"],ds["AttributeType"])
      output += self.xmlInt.inc() + '<DataItem Dimensions="%d %d %d %s" Format="HDF">\n'% \
           (h5dict["subdomainSize"][0],h5dict["subdomainSize"][1],h5dict["subdomainSize"][2], ds["rankString"])
      output += self.xmlInt.cur() + '%s:/%s\n'%(h5dict["pathToHDF5"],ds["attributeName"])
      output += self.xmlInt.dec() + '</DataItem>\n'
      output += self.xmlInt.dec() + '</Attribute>\n'

    output += self.xmlInt.dec() + '</Grid>\n'
    self.cached_subgrid_output[h5dict["processorId"]] += output

  def __init__(self, fnameString, processorStrings, iterDir):
    """
    Reads file fitting the patters and saves an XMF file to the directory of the HDF5s.
    
        fnameString = './Fluid_00000000_%s.h5'
        processorStrings = ['p000', 'p001', 'p002']
        createXDMF(fnameString, processorStrings)

        Important: reads as many files as number of processors in the first file
    """
    fnameToSave =  ( (fnameString%(""))[:-6] + '.xmf' )
    fnameToSave = fnameToSave.replace('/hdf5/','/').replace(".p.",".").replace("/" + iterDir + "/", "")
    processorStrings = sorted(processorStrings)
    self.fnameString = fnameString
    if os.path.isfile(fnameToSave):
      print("%s (existed)" % (fnameToSave))
    self.xdmf_file = open(fnameToSave, "w")
    self.xmlInt = XMLIndentation()
    self.output = ""
    self.cached_subgrid_output = {} #neccessary to sort subgrids per mpi proc

    self.writeHeader()
    for p in processorStrings:
      self.write_cached_subgrid(p)
    for key,value in self.cached_subgrid_output.items():
      self.output += self.xmlInt.inc() + '<Grid Name="CPU %s" GridType="Collection">\n'%(key)
      self.output += value
      self.output += self.xmlInt.dec() + '</Grid>\n'

    self.writeFooter()
    self.xdmf_file.write(self.output)
    self.xdmf_file.close()
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
