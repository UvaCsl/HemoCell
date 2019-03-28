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
	cid =None,
	bcid =None,
	velocity =None
	):
		self.position = position
		self.area = area
		self.volume = volume
		self.atomicblock = atomicblock
		self.cid = cid
		self.bcid = bcid
		self.velocity = velocity


def open_csv_files(r=True,p=True,ct3=True,begin=0,end=0,timestep=100000,datapath=".",rbcname="RBC_HO",pltname="PLT",ct3name="RBC_stiff"):

	RBC = []
	PLT = []
	CT3 = []
	print("reading in ... ")
	print(datapath)
	for t in range(begin,end,timestep):

		#timepath = datapath+str(t).zfill(12)
		#timepath = datapath+str(t)
		#print(str(t).zfill(12))
		rposition,rarea,rvolume,ratomicblock,rcid,rbcid,rvelocity = ([] for i in range(7))
		pposition,parea,pvolume,patomicblock,pcid,pbcid,pvelocity = ([] for i in range(7))
		cposition,carea,cvolume,catomicblock,ccid,cbcid,cvelocity = ([] for i in range(7))
		if r:
			with open(datapath+rbcname+'.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
			#with open(timepath+'/CellInfo_RBC_HO.'+str(n)+'.csv', 'r') as csvfile:
				next(csvfile,None)
				tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
				for row in tmpfiledata:
					rposition.extend([[row[0],row[1],row[2]]])
					rarea.extend([row[3]])
					rvolume.extend([row[4]])
					ratomicblock.extend([row[5]])
					rcid.extend([row[6]])
					rbcid.extend([row[7]])
					rvelocity.extend([[row[8],row[9],row[10]]])

		if p:
			with open(datapath+pltname+'.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
			#with open(timepath+'/CellInfo_PLT.'+str(n)+'.csv', 'r') as csvfile:
				next(csvfile,None)
				tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
				for row in tmpfiledata:
					pposition.extend([[row[0],row[1],row[2]]])
					parea.extend([row[3]])
					pvolume.extend([row[4]])
					patomicblock.extend([row[5]])
					pcid.extend([row[6]])
					pbcid.extend([row[7]])
					pvelocity.extend([[row[8],row[9],row[10]]])
		if ct3:
			with open(datapath+ct3name+'.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
			#with open(timepath+'/CellInfo_RBC_stiff.'+str(n)+'.csv', 'r') as csvfile:
				next(csvfile,None)
				tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
				for row in tmpfiledata:
					cposition.extend([[row[0],row[1],row[2]]])
					carea.extend([row[3]])
					cvolume.extend([row[4]])
					catomicblock.extend([row[5]])
					ccid.extend([row[6]])
					cbcid.extend([row[7]])
					cvelocity.extend([[row[8],row[9],row[10]]])
		if r:
			RBC.append(CELL(rposition,rarea,rvolume,ratomicblock,rcid,rbcid,rvelocity))
		if p:
			PLT.append(CELL(pposition,parea,pvolume,patomicblock,pcid,pbcid,pvelocity))
		if ct3:
			CT3.append(CELL(cposition,carea,cvolume,catomicblock,ccid,cbcid,cvelocity))



	return(np.array(RBC),np.array(PLT),np.array(CT3))

def open_csv_files_OLD(r=True,p=True,ct3=True,begin=0,end=0,timestep=100000,datapath=".",nprocs=1,rbcname="RBC_HO",pltname="PLT",ct3name="RBC_stiff"):

	RBC = []
	PLT = []
	CT3 = []
	print("reading in ... ")
	print(datapath)
	for t in range(begin,end,timestep):

		timepath = datapath+str(t).zfill(12)
		#timepath = datapath+str(t)
		#print(str(t).zfill(12))
		rposition,rarea,rvolume,ratomicblock,rcid,rbcid,rvelocity = ([] for i in range(7))
		pposition,parea,pvolume,patomicblock,pcid,pbcid,pvelocity = ([] for i in range(7))
		cposition,carea,cvolume,catomicblock,ccid,cbcid,cvelocity = ([] for i in range(7))
		for n in range(0,nprocs):
			if r:
				#with open(datapath+'RBC_HO.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
				with open(timepath+'/CellInfo_'+rbcname+'.'+str(n)+'.csv', 'r') as csvfile:
					next(csvfile,None)
					tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
					for row in tmpfiledata:
						rposition.extend([[row[0],row[1],row[2]]])
						rarea.extend([row[3]])
						rvolume.extend([row[4]])
						ratomicblock.extend([row[5]])
						rcid.extend([row[6]])
			if p:
				#with open(datapath+'PLT.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
				with open(timepath+'/CellInfo_'+pltname+'.'+str(n)+'.csv', 'r') as csvfile:
					next(csvfile,None)
					tmpfiledata = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
					for row in tmpfiledata:
						pposition.extend([[row[0],row[1],row[2]]])
						parea.extend([row[3]])
						pvolume.extend([row[4]])
						patomicblock.extend([row[5]])
						pcid.extend([row[6]])
			if ct3:
				#with open(datapath+'RBC_stiff.'+str(t).zfill(12)+'.csv', 'r') as csvfile:
				with open(timepath+'/CellInfo_'+ct3name+'.'+str(n)+'.csv', 'r') as csvfile:
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