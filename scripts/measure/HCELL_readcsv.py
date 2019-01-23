import csv
import numpy as np
import pdb

class CELL(object):
	"""__init__() functions as the class constructor"""
	def __init__(self,
	position = None,
	area = None,
	volume = None,
	atomicblock =None,
	cid =None
	):
		self.position = position
		self.area = area
		self.volume = volume
		self.atomicblock = atomicblock
		self.cid = cid

def open_csv_files(r=True,p=True,ct3=True,begin=0,end=0,timestep=100000,datapath=".",nprocs=1):

	RBC = []
	PLT = []
	CT3 = []
	for t in range(begin,end,timestep):

		timepath = datapath+str(t).zfill(12)
		#timepath = datapath+str(t)
		print("reading in ... ")
		print(timepath)
		rposition,rarea,rvolume,ratomicblock,rcid = ([] for i in range(5))
		pposition,parea,pvolume,patomicblock,pcid = ([] for i in range(5))
		cposition,carea,cvolume,catomicblock,ccid = ([] for i in range(5))

		for n in range(nprocs):
			if r:
				with open(timepath+'/CellInfo_RBC_HO.'+str(n)+'.csv', 'r') as csvfile:
					next(csvfile,None)
					tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
					for row in tmpfiledata:
						rposition.extend([[row[0],row[1],row[2]]])
						rarea.extend([row[3]])
						rvolume.extend([row[4]])
						ratomicblock.extend([row[5]])
						rcid.extend([row[6]])
			if p:
				with open(timepath+'/CellInfo_PLT.'+str(n)+'.csv', 'r') as csvfile:
					next(csvfile,None)
					tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
					for row in tmpfiledata:
						pposition.extend([[row[0],row[1],row[2]]])
						parea.extend([row[3]])
						pvolume.extend([row[4]])
						patomicblock.extend([row[5]])
						pcid.extend([row[6]])
			if ct3:
				with open(timepath+'/CellInfo_RBC_stiff.'+str(n)+'.csv', 'r') as csvfile:
					next(csvfile,None)
					tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
					for row in tmpfiledata:
						cposition.extend([[row[0],row[1],row[2]]])
						carea.extend([row[3]])
						cvolume.extend([row[4]])
						catomicblock.extend([row[5]])
						ccid.extend([row[6]])
		if r:
			RBC.append(CELL(rposition,rarea,rvolume,ratomicblock,rcid))
		if p:
			PLT.append(CELL(pposition,parea,pvolume,patomicblock,pcid))
		if ct3:
			CT3.append(CELL(cposition,carea,cvolume,catomicblock,ccid))



	return(np.array(RBC),np.array(PLT),np.array(CT3))