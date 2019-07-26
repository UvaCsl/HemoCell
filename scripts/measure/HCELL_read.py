#HCELL_read.py 
#Written by Ben Czaja 26.09.2018
#This code is parallelized to read in mulitple Nprocs at the same time
#This code is NOT parallelized on time
#If the number of run cores ==1 then the program will read in serial automatically
import numpy as np
import h5py
import csv
from multiprocessing import Pool
from functools import partial
import errno
import os
#===============================================================================
#===============================================================================

#Fluid class of attributes from the hdf5 files
class FLUID(object):
	def __init__(self, dt=None,dx=None,dxdydz=None,iteration=None,numberofcells=None,processorId=None,relativeposition=None,subdomainsize=None,position=None,velocity=None):
		self.dt = dt 
		self.dx = dx
		self.dxdydz = dxdydz
		self.iteration = iteration
		self.numberofcells = numberofcells
		self.processorId = processorId
		self.relativeposition = relativeposition
		self.subdomainsize = subdomainsize
		self.position = position
		self.velocity = velocity
#Cell class of attributes from the hdf5 files (Can be RBCs, WBCs, Platlets from HemoCell)
class LSP_CELL(object):
	def __init__(self,
	 dt=None,
	 dx=None,
	 iteration=None,
	 numberofparticles=None,
	 numberofprocessors=None,
	 numberoftriangles=None,
	 processorId=None,
	 Farea=None,
	 Fbend=None,
	 Flink=None,
	 Ftotal=None,
	 Fvisc=None,
	 Fvol=None,
	 triangles=None,
	 position=None,
	 cid =None
	 ):

		self.dt = dt
		self.dx = dx
		self.iteration = iteration
		self.numberofparticles =numberofparticles
		self.numberofprocessors = numberofprocessors
		self.numberoftriangles = numberoftriangles
		self.processorId= processorId
		self.Farea = Farea
		self.Fbend = Fbend
		self.Flink = Flink
		self.Ftotal = Ftotal
		self.Fvisc = Fvisc        
		self.Fvol = Fvol
		self.triangles = triangles
		self.position = position
		self.cid = cid

class CSVCELL(object):
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




#The read function
def core_read_processor(f,r,p,ct3,fluidoutputs,celloutputs,t,datapath,rbcname,pltname,ct3name,n):
	fdt,fdx,fdxdydz,fiter,fNpart,fprocid,frelpos,fsubdomsize,fpos,fvel = ([] for i in range(10))
	rdt,rdx,riter,rNpart,rNproc,rNtri,rprocid,rFarea,rFbend,rFlink,rFtotal,rFvisc,rFvol,rtriangles,rpos,rcid = ([] for i in range(16))
	pdt,pdx,piter,pNpart,pNproc,pNtri,pprocid,pFarea,pFbend,pFlink,pFtotal,pFvisc,pFvol,ptriangles,ppos,pcid = ([] for i in range(16))
	ct3dt,ct3dx,ct3iter,ct3Npart,ct3Nproc,ct3Ntri,ct3procid,ct3Farea,ct3Fbend,ct3Flink,ct3Ftotal,ct3Fvisc,ct3Fvol,ct3triangles,ct3pos,ct3cid = ([] for i in range(16))

	time_path = datapath+str(t).zfill(12)
	#time_path = datapath+str(t)
	if n == 0:
		print("reading in files ...")
		print(time_path)
	#Read in each processer filename
	#Read all fluid files
	if f:
		#open file
		#Read HemoCell v2.0 and v1.0 output. Try to catch files that do not exist. 
		try:
			fluid_file = h5py.File(time_path+"/Fluid." + str(t).zfill(12) +".p." +str(n) +".h5","r")	
		except:
			try:
				fluid_file = h5py.File(datapath+str(t)+"/Fluid." + str(t) +".p." +str(n) +".h5","r")
			except (OSError, IOError):
				raise

		print(fluid_file)

		#get data
		if "dt" in fluidoutputs:
			fdt.extend(fluid_file["dt"])
		if "dx" in fluidoutputs:
			fdx.extend(fluid_file["dx"])				
		if "dxdydz" in fluidoutputs:
			fdxdydz.extend(fluid_file["dxdydz"])
		if "iteration" in fluidoutputs:
			fiter.extend(fluid_file["iteration"])
		if "numberofcells" in fluidoutputs:
			fNpart.extend(fluid_file["numberOfcells"])	
		if "processorId" in fluidoutputs:
			fprocid.extend(fluid_file["processorId"])
		if "relativeposition" in fluidoutputs:
			frelpos.extend(fluid_file["relativePosition"])
		if "subdomainsize" in fluidoutputs:	
			fsubdomsize.extend(fluid_file["subdomainSize"])
		if "Velocity" in fluidoutputs:
			pos =[]
			vel = []
			tempvel = np.array(fluid_file["Velocity"])
			relpos = fluid_file.attrs.get('relativePosition')
	
			#THESE INDICIIES ARE REVERSED BECASUE OF PARAVIEW
			xblocks = np.shape(tempvel)[2]
			yblocks = np.shape(tempvel)[1]
			zblocks = np.shape(tempvel)[0]

			for xpos in range(xblocks):
				for ypos in range(yblocks):
					for zpos in range(zblocks):
						vel.append(np.array([tempvel[zpos][ypos][xpos][0],tempvel[zpos][ypos][xpos][1],tempvel[zpos][ypos][xpos][2]]))
						pos.append(np.array([xpos+relpos[2],ypos+relpos[1],zpos+relpos[0]]))
						#print(pos)

			fpos.extend(pos)
			fvel.extend(vel)

		#close file
		fluid_file.close()

	if r:		
		#open file
		#Read HemoCell v2.0 and v1.0 output. Try to catch files that do not exist. 
		try:
			RBC_file = h5py.File(time_path+"/"+rbcname+"." + str(t).zfill(12) +".p." +str(n) +".h5","r")	
		except:
			try:
				RBC_file = h5py.File(datapath+str(t)+"/"+rbcname+"." + str(t) +".p." +str(n) +".h5","r")
			except (OSError, IOError):
				raise

		print(RBC_file)
		#get data
		if "dt" in celloutputs:
			rdt.extend(RBC_file["dt"])
		if "dx" in celloutputs:
			rdx.extend(RBC_file["dx"])
		if "iteration" in celloutputs:			
			riter.extend(RBC_file["iteration"])
		if "numberOfParticles" in celloutputs:		
			rNpart.extend(RBC_file["numberOfParticles"])	
		if "numberOfProcessors" in celloutputs:		
			rNproc.extend(RBC_file["numberOfProcessors"])
		if "numberOfTriangles" in celloutputs:
			rNtri.extend(RBC_file["numberOfTriangles"])
		if "processorId" in celloutputs:
			rprocid.extend(RBC_file["processorId"])
		if "Area force" in celloutputs:
			rFarea.extend(RBC_file["Area force"])
		if "Bending force" in celloutputs:
			rFbend.extend(RBC_file["Bending force"])
		if "Link force" in celloutputs:
			rFlink.extend(RBC_file["Link force"])
		if "Total force" in celloutputs:
			rFtotal.extend(RBC_file["Total force"])
		if "Viscous force" in celloutputs:
			rFvisc.extend(RBC_file["Viscous force"])
		if "Volume force" in celloutputs:
			rFvol.extend(RBC_file["Volume force"])
		if "Triangles" in celloutputs:
			rtriangles.extend(RBC_file["Triangles"])
		if "Position" in celloutputs:
			rpos.extend(RBC_file["Position"])
		if "Cell Id" in celloutputs:
			rcid.extend(RBC_file["Cell Id"])

		#close file
		RBC_file.close()
	if p:		
		#open file
		#Read HemoCell v2.0 and v1.0 output. Try to catch files that do not exist. 
		try:
			platelet_file = h5py.File(time_path+"/"+pltname+"." + str(t).zfill(12) +".p." +str(n) +".h5","r")	
		except:
			try:
				platelet_file = h5py.File(datapath+str(t)+"/"+pltname+"." + str(t) +".p." +str(n) +".h5","r")
			except (OSError, IOError):
				raise

		print(platelet_file)
		#get data
		if "dt" in celloutputs:
			pdt.extend(platelet_file["dt"])
		if "dx" in celloutputs:
			pdx.extend(platelet_file["dx"])
		if "iteration" in celloutputs:			
			piter.extend(platelet_file["iteration"])
		if "numberOfParticles" in celloutputs:		
			pNpart.extend(platelet_file["numberOfParticles"])	
		if "numberOfProcessors" in celloutputs:		
			pNproc.extend(platelet_file["numberOfProcessors"])
		if "numberOfTriangles" in celloutputs:
			pNtri.extend(platelet_file["numberOfTriangles"])
		if "processorId" in celloutputs:
			pprocid.extend(platelet_file["processorId"])
		if "Area force" in celloutputs:
			pFarea.extend(platelet_file["Area force"])
		if "Bending force" in celloutputs:
			pFbend.extend(platelet_file["Bending force"])
		if "Link force" in celloutputs:
			pFlink.extend(platelet_file["Link force"])
		if "Total force" in celloutputs:
			pFtotal.extend(platelet_file["Total force"])
		if "Viscous force" in celloutputs:
			pFvisc.extend(platelet_file["Viscous force"])
		if "Volume force" in celloutputs:
			pFvol.extend(platelet_file["Volume force"])
		if "Triangles" in celloutputs:
			ptriangles.extend(platelet_file["Triangles"])
		if "Position" in celloutputs:
			ppos.extend(platelet_file["Position"])
		if "Cell Id" in celloutputs:
			pcid.extend(platelet_file["Cell Id"])

		#close file
		platelet_file.close()


	if ct3:		
		#open file
		#Read HemoCell v2.0 and v1.0 output. Try to catch files that do not exist. 		
		try:
			ct3_file = h5py.File(time_path+"/"+ct3name+"." + str(t).zfill(12) +".p." +str(n) +".h5","r")	
		except:
			try:
				ct3_file = h5py.File(datapath+str(t)+"/"+ct3name+"."+ str(t) +".p." +str(n) +".h5","r")
			except (OSError, IOError):
				raise

		print(ct3_file)
		#get data
		if "dt" in celloutputs:
			ct3dt.extend(ct3_file["dt"])
		if "dx" in celloutputs:
			ct3dx.extend(ct3_file["dx"])
		if "iteration" in celloutputs:			
			ct3iter.extend(ct3_file["iteration"])
		if "numberOfParticles" in celloutputs:		
			ct3Npart.extend(ct3_file["numberOfParticles"])	
		if "numberOfProcessors" in celloutputs:		
			ct3Nproc.extend(ct3_file["numberOfProcessors"])
		if "numberOfTriangles" in celloutputs:
			ct3Ntri.extend(ct3_file["numberOfTriangles"])
		if "processorId" in celloutputs:
			ct3procid.extend(ct3_file["processorId"])
		if "Area force" in celloutputs:
			ct3Farea.extend(ct3_file["Area force"])
		if "Bending force" in celloutputs:
			ct3Fbend.extend(ct3_file["Bending force"])
		if "Link force" in celloutputs:
			ct3Flink.extend(ct3_file["Link force"])
		if "Total force" in celloutputs:
			ct3Ftotal.extend(ct3_file["Total force"])
		if "Viscous force" in celloutputs:
			ct3Fvisc.extend(ct3_file["Viscous force"])
		if "Volume force" in celloutputs:
			ct3Fvol.extend(ct3_file["Volume force"])
		if "Triangles" in celloutputs:
			ct3triangles.extend(ct3_file["Triangles"])
		if "Position" in celloutputs:
			ct3pos.extend(ct3_file["Position"])
		if "Cell Id" in celloutputs:
			ct3cid.extend(ct3_file["Cell Id"])

		#close file
		ct3_file.close()

	return(fdt,fdx,fdxdydz,fiter,fNpart,fprocid,frelpos,fsubdomsize,fpos,fvel,
		rdt,rdx,riter,rNpart,rNproc,rNtri,rprocid,rFarea,rFbend,rFlink,rFtotal,rFvisc,rFvol,rtriangles,rpos,rcid,
		pdt,pdx,piter,pNpart,pNproc,pNtri,pprocid,pFarea,pFbend,pFlink,pFtotal,pFvisc,pFvol,ptriangles,ppos,pcid,
		ct3dt,ct3dx,ct3iter,ct3Npart,ct3Nproc,ct3Ntri,ct3procid,ct3Farea,ct3Fbend,ct3Flink,ct3Ftotal,ct3Fvisc,ct3Fvol,ct3triangles,ct3pos,ct3cid)



#The wrapper around the read function 
#This will read in a time range and all the run processors
def open_hdf5_files(f=True,r=True,p=True,ct3=True,fluidoutputs=[""],celloutputs=[""],read_procs=1,begin=0,end=0,timestep=100000,datapath=".",rbcname="RBC_HO",pltname="PLT",ct3name="RBC_stiff",nprocs=1):

	fluid =[]
	rbc =[]
	platelet =[]
	celltype3 =[]

	for t in range(begin,end,timestep):
		fdt,fdx,fdxdydz,fiter,fNpart,fprocid,frelpos,fsubdomsize,fpos,fvel = ([] for i in range(10))
		rdt,rdx,riter,rNpart,rNproc,rNtri,rprocid,rFarea,rFbend,rFlink,rFtotal,rFvisc,rFvol,rtriangles,rpos,rcid = ([] for i in range(16))
		pdt,pdx,piter,pNpart,pNproc,pNtri,pprocid,pFarea,pFbend,pFlink,pFtotal,pFvisc,pFvol,ptriangles,ppos,pcid = ([] for i in range(16))
		ct3dt,ct3dx,ct3iter,ct3Npart,ct3Nproc,ct3Ntri,ct3procid,ct3Farea,ct3Fbend,ct3Flink,ct3Ftotal,ct3Fvisc,ct3Fvol,ct3triangles,ct3pos,ct3cid = ([] for i in range(16))

		#Check if the  run procs are actually >1
		if read_procs >1:
			#If they are then split them on the Nread_procs
			pool = Pool(read_procs)
			iter_procs = (i for i in range(nprocs))
			#Prepare the function for pool.map
			read_function = partial(core_read_processor,f,r,p,ct3,fluidoutputs,celloutputs,t,datapath,rbcname,pltname,ct3name)
			#Read in parallel
			result = np.array(pool.map(read_function,iter_procs))
			pool.close()
			pool.join()
	
			#Crappy rejoining stuff
			for i in range(nprocs):
	
				fdt.extend(result[i,0])
				fdx.extend(result[i,1])
				fdxdydz.extend(result[i,2])
				fiter.extend(result[i,3])
				fNpart.extend(result[i,4])
				fprocid.extend(result[i,5])
				frelpos.extend(result[i,6])
				fsubdomsize.extend(result[i,7])
				fpos.extend(result[i,8])
				fvel.extend(result[i,9])
	
				rdt.extend(result[i,10])
				rdx.extend(result[i,11])
				riter.extend(result[i,12])
				rNpart.extend(result[i,13])
				rNproc.extend(result[i,14])
				rNtri.extend(result[i,15])
				rprocid.extend(result[i,16])
				rFarea.extend(result[i,17])
				rFbend.extend(result[i,18])
				rFlink.extend(result[i,19])
				rFtotal.extend(result[i,20])
				rFvisc.extend(result[i,21])
				rFvol.extend(result[i,22])
				rtriangles.extend(result[i,23])
				rpos.extend(result[i,24])
				rcid.extend(result[i,25])
	
				pdt.extend(result[i,26])
				pdx.extend(result[i,27])
				piter.extend(result[i,28])
				pNpart.extend(result[i,29])
				pNproc.extend(result[i,30])
				pNtri.extend(result[i,31])
				pprocid.extend(result[i,32])
				pFarea.extend(result[i,33])
				pFbend.extend(result[i,34])
				pFlink.extend(result[i,35])
				pFtotal.extend(result[i,36])
				pFvisc.extend(result[i,37])
				pFvol.extend(result[i,38])
				ptriangles.extend(result[i,39])
				ppos.extend(result[i,40])
				pcid.extend(result[i,41])

				ct3dt.extend(result[i,42])
				ct3dx.extend(result[i,43])
				ct3iter.extend(result[i,44])
				ct3Npart.extend(result[i,45])
				ct3Nproc.extend(result[i,46])
				ct3Ntri.extend(result[i,47])
				ct3procid.extend(result[i,48])
				ct3Farea.extend(result[i,49])
				ct3Fbend.extend(result[i,50])
				ct3Flink.extend(result[i,51])
				ct3Ftotal.extend(result[i,52])
				ct3Fvisc.extend(result[i,53])
				ct3Fvol.extend(result[i,54])
				ct3triangles.extend(result[i,54])
				ct3pos.extend(result[i,56])						
				ct3cid.extend(result[i,57])	

		elif nprocs ==1:
			#Read in serial if only 1 run processor
			result = core_read_processor(f,r,p,ct3,fluidoutputs,celloutputs,t,datapath,ct3name,0)

			fdt         = result[0]
			fdx         = result[1]
			fdxdydz     = result[2]
			fiter       = result[3]
			fNpart      = result[4]
			fprocid     = result[5]
			frelpos     = result[6]
			fsubdomsize = result[7]
			fpos        = result[8]
			fvel        = result[9]

			rdt         = result[10]
			rdx         = result[11]
			riter       = result[12]
			rNpart      = result[13]
			rNproc      = result[14]
			rNtri       = result[15]
			rprocid     = result[16]
			rFarea      = result[17]
			rFbend      = result[18]
			rFlink      = result[19]
			rFtotal     = result[20]
			rFvisc      = result[21]
			rFvol       = result[22]
			rtriangles  = result[23]
			rpos        = result[24]
			rcid        = result[25]

			pdt         = result[26]
			pdx         = result[27]
			piter       = result[28]
			pNpart      = result[29]
			pNproc      = result[30]
			pNtri       = result[31]
			pprocid     = result[32]
			pFarea      = result[33]
			pFbend      = result[34]
			pFlink      = result[35]
			pFtotal     = result[36]
			pFvisc      = result[37]
			pFvol       = result[38]
			ptriangles  = result[39]
			ppos        = result[40]
			pcid        = result[41]

			ct3dt         = result[42]
			ct3dx         = result[43]
			ct3iter       = result[44]
			ct3Npart      = result[45]
			ct3Nproc      = result[46]
			ct3Ntri       = result[47]
			ct3procid     = result[48]
			ct3Farea      = result[49]
			ct3Fbend      = result[50]
			ct3Flink      = result[51]
			ct3Ftotal     = result[52]
			ct3Fvisc      = result[53]
			ct3Fvol       = result[54]
			ct3triangles  = result[55]
			ct3pos        = result[56]
			ct3cid        = result[57]

		if f:
			fluid.append(FLUID(fdt,fdx,fdxdydz,fiter,fNpart,fprocid,frelpos,fsubdomsize,fpos,fvel))
		if r:
			rbc.append(LSP_CELL(rdt,rdx,riter,rNpart,rNproc,rNtri,rprocid,rFarea,rFbend,rFlink,rFtotal,rFvisc,rFvol,rtriangles,rpos,rcid))
		if p:
			platelet.append(LSP_CELL(pdt,pdx,piter,pNpart,pNproc,pNtri,pprocid,pFarea,pFbend,pFlink,pFtotal,pFvisc,pFvol,ptriangles,ppos,pcid))
		if ct3:
			celltype3.append(LSP_CELL(ct3dt,ct3dx,ct3iter,ct3Npart,ct3Nproc,ct3Ntri,ct3procid,ct3Farea,ct3Fbend,ct3Flink,ct3Ftotal,ct3Fvisc,ct3Fvol,ct3triangles,ct3pos,ct3cid))


	return np.array(fluid),np.array(rbc),np.array(platelet),np.array(celltype3)



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
			RBC.append(CSVCELL(rposition,rarea,rvolume,ratomicblock,rcid,rbcid,rvelocity))
		if p:
			PLT.append(CSVCELL(pposition,parea,pvolume,patomicblock,pcid,pbcid,pvelocity))
		if ct3:
			CT3.append(CSVCELL(cposition,carea,cvolume,catomicblock,ccid,cbcid,cvelocity))



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
			RBC.append(CSVCELL(rposition,rarea,rvolume,ratomicblock,rcid))
		if p:
			PLT.append(CSVCELL(pposition,parea,pvolume,patomicblock,pcid))
		if ct3:
			CT3.append(CSVCELL(cposition,carea,cvolume,catomicblock,ccid))



	return(np.array(RBC),np.array(PLT),np.array(CT3))