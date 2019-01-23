import numpy as np


def rectangle_volume_fraction(cell,CELLVOL,envelope=60.):
    '''Calculated X,Y,Z cell volume fraction in a rectangular chamber

    Keyword arguments:
    cell -- cell (class), can be indexed over time
    envelope -- envelope (float) region of the simulation domain
    CELLVOL -- volume (float) of the cell in micrometers^3

    Returns:
    X,Y,Z (numpy arrays) cell distribution histograms normalized to cell volume fraction.
    '''

    steps = 1  #step size in micrometers
    binsx = np.arange(-envelope,np.max(np.array(cell[0].position)[:,0])+envelope,steps)
    binsy = np.arange(-envelope,np.max(np.array(cell[0].position)[:,1])+envelope,steps)
    binsz = np.arange(-1.5,np.max(np.array(cell[0].position)[:,2])+steps+10,steps)
    
    Xsimlength = np.max(binsx) - np.min(binsx)
    Ysimlength = np.max(binsy) - np.min(binsy)
    Zsimlength = np.max(binsz) - np.min(binsz)

    Nlsps_per_cell = cell_statistics(cell[0])[0]
    
    cell_histx = []
    cell_histy = []
    cell_histz = []
    
    for i in range(len(cell)):
        #X slice
        stripVOL = steps*Ysimlength*Zsimlength
        Nlsps_in_strip = np.histogram(np.unique(np.array(cell[i].position)[:,0]),bins=binsx)[0]
        cell_histx.append(Nlsps_in_strip*CELLVOL/Nlsps_per_cell/stripVOL)
        
        #Y slice
        stripVOL = steps*Zsimlength*Xsimlength
        Nlsps_in_strip = np.histogram(np.unique(np.array(cell[i].position)[:,1]),bins=binsy)[0]
        cell_histy.append(Nlsps_in_strip*CELLVOL/Nlsps_per_cell/stripVOL)
        
        #Z slice
        stripVOL = steps*Ysimlength*Xsimlength
        Nlsps_in_strip = np.histogram(np.unique(np.array(cell[i].position)[:,2]),bins=binsz)[0]
        cell_histz.append(Nlsps_in_strip*CELLVOL/Nlsps_per_cell/stripVOL)
        
    cell_histx = np.mean(cell_histx,axis=0)
    cell_histy = np.mean(cell_histy,axis=0)
    cell_histz = np.mean(cell_histz,axis=0)

    return(cell_histx,cell_histy,cell_histz)






def pipeflow_MSD_cell_centers(cell,X,Y,Z,dx=1,refdir=0):
    R = np.sqrt(Z**2 +Y**2)
    L = X
    R0 = np.sqrt((Z*0.5)**2 +(Y*0.5)**2)
    R -= R0

    print("N of Cells deleted = ",len(cell[0].cid) - len(cell[-1].cid))
    print("% of Cells deleted = ",(len(cell[0].cid) - len(cell[-1].cid))/len(cell[0].cid)*100.)

    #Get final info to make sure that you calculate MSD for cells that are not deleted
    pos = np.array(cell[-1].position)
    cid_final = np.array(cell[-1].cid)
    r_final = np.sqrt((pos[:,1] - Y*.5)**2 + (pos[:,2] - Z*0.5)**2)
    tmp_final = np.array([r_final,cid_final])
    #sorted on increasing CID
    tmp_final = tmp_final.T[tmp_final.T[:,1].argsort()].T


    radial_and_cids_overtime = []
    for t in range(len(cell)):

        #get positions
        pos = np.array(cell[t].position)
        #Convert to r and make the center of tube the origin
        r = np.sqrt((pos[:,1] - Y*.5)**2 + (pos[:,2] - Z*0.5)**2)
        #create array of positions (r) and cids
        cid = np.array(cell[t].cid)
        tmp = np.array([r,cid])

        #now sort them on increasing CID
        tmp = tmp.T[tmp.T[:,1].argsort()].T

        tmp_rs_cids = []
        for i in tmp_final[1]:
            if len(np.where(tmp[1] == i)[0]) >0:
                idx = np.where(tmp[1] == i)[0][0]
                tmp_rs_cids.append(np.array([tmp[0][idx],tmp[1][idx]]))
        radial_and_cids_overtime.append(np.array(tmp_rs_cids))

    radial_and_cids_overtime = np.array(radial_and_cids_overtime)
    

    #now calculate the mean squared displacements
    rposses = []
    for i in range(len(radial_and_cids_overtime)):
        rposses.append(radial_and_cids_overtime[i][:,0])
    rposses = np.array(rposses)

    #This is how you get it for the entire ensemble
    #11.23.2018
    diff = np.diff(rposses.T)
    diff_sq = diff**2
    MSD = np.mean(diff_sq,axis=0)

    #BEN you need to figure out how to get this for different locations within the tube
    #11.23.2018
    return(MSD)




def pipeflow_radial_volumefraction(cell,X,Y,Z,dx,refdir=0):

    lsp_per_cell,N_cells = cell_statistics(cell[0])
    
    if lsp_per_cell >=640 and lsp_per_cell <=650:
        cell_vol = 90
    elif lsp_per_cell > 60 and lsp_per_cell < 70:
        cell_vol = 11

    R = np.sqrt(Z**2 +Y**2)
    L = X
    R0 = np.sqrt((Z*0.5)**2 +(Y*0.5)**2)
    R -= R0
    
    rbins = []
    for r in np.arange(R):
        print (r)
        tmpbin =[]
        for i in range(len(cell)):
            pos = np.array(cell[i].position)*dx*1e6
            
            xcell = pos[:,0]
            ycell = pos[:,1] - 0.5*Y
            zcell = pos[:,2] - 0.5*Z
            rcell = np.sqrt(ycell**2+zcell**2)
            #This is for RBCs being saved 3 times
            rcell = np.unique(rcell)
            rmask = (rcell <= r+1) & (rcell > r)
            if len(rcell[rmask]) == 0:
                tmpbin.append(0)
            else:
                tmpbin.append(len(rcell[rmask]))
        slice_area = np.pi*(r+1)**2 - np.pi*r**2
        # N_lsps per slice area
        #rbins.append(np.mean(tmpbin)*90./slice_area)
        rbins.append(np.mean(tmpbin)/lsp_per_cell*cell_vol/slice_area/L)
    
    R_hematocrit = np.array(rbins)
    
    return(R_hematocrit)



def pipe_volumefraction(cell,R,L):
    
    lsp_per_cell,N_cells = cell_statistics(cell)
    
    if lsp_per_cell >=640 and lsp_per_cell <=650:
        cell_vol = 90
    elif lsp_per_cell > 60 and lsp_per_cell < 70:
        cell_vol = 11

    pipe_volume = np.pi*R*R*L
    
    return(N_cells*cell_vol/pipe_volume)
    
def cell_statistics(cell):

    pos = np.array(cell.position)
    cids = np.array(cell.cid)
    cid = np.unique(cell.cid)
    mask = np.where(cid[0] == cids)[0]
    
    #This is to make sure there are no copies with the same CID
    lsp_per_cell = len(np.unique(pos[mask,0]))
    N_cells = len(np.unique(pos[:,0]))/lsp_per_cell
    
    
    return(lsp_per_cell,N_cells)


def pipeflow_radial_shearrate(fluid,X,Y,Z,dx,dt,refdir=0):
    
    R = np.sqrt(Y**2+Z**2)
    R0 = np.sqrt((Y*0.5)**2+(Z*0.5)**2)
    R -= R0
    L = X
    
    vr = []
    for r in np.arange(50):
        tmpvbin = []
        for i in range(len(fluid)):
            pos = np.array(fluid[i].position)*dx*1e6
            vel = np.array(fluid[i].velocity)
            rx = np.array(pos[:,0])
            ry = np.array(pos[:,1]) - 0.5*Y
            rz = np.array(pos[:,2]) - 0.5*Z
            vx = np.array(vel[:,0])
            vy = np.array(vel[:,1])
            vz = np.array(vel[:,2])
            
            r_fluid  = np.sqrt(rz**2 + ry**2)
            vr_fluid = np.sqrt(vz**2 + vy**2)
            
            r_mask = (r_fluid > r) & (r_fluid <= r+1)
            
            #pdb.set_trace()
            #print(np.mean(vr_fluid[r_mask]),np.max(vr_fluid[r_mask]),np.min(vr_fluid[r_mask]))
            #tmpvbin.append(np.sqrt(vz[r_mask]**2 + vy[r_mask]**2))
            tmpvbin.append(np.mean(vx[r_mask]))
            #print(np.mean(vx[r_mask]))
        vr.append(tmpvbin)
            
            
    #pdb.set_trace()
    vr = np.mean(vr,axis=1)*dx/dt
    SR = np.append(vr,0)
    #1e6 converting from micrometer to meter
    SR = np.gradient(SR*1e6)
            
    return(vr,SR)
    

#vr,SR = pipeflow_radial_shearrate(fluids,X=200,Y=100,Z=100,dx=dx,dt=dt)
    

def rectangle_velocity_profile(fluid,dx,dt,directvel,directpos,steppos=1.0):
    '''Get velocity in a direction from a rectangle

    Keyword arguments:
    fluid -- fluid (class), can be indexed over time
    dx --  lattice spatial resolution (float) in SI
    dt -- lattice temporal resolution (float) in SI
    directvel -- velocity direction (int) X=0,Y=1,Z=2
    directpos -- position axis (int) X=0,Y=1,Z=2
    steppos -- spatial resolution (float) to measure velocity

    Returns:
    Velocity (numpy array  dtype=float)
    '''
    #get the time average
    avgvel = np.array(fluid[0].velocity)
    for i in range(1,len(fluid)):
        avgvel += np.array(fluid[i].velocity)
    avgvel /= float(len(fluid))

    vel_mean = []
    pos = np.array(fluid[0].position)

    rectangle_side = np.arange(np.min(pos[:,directpos]),np.max(pos[:,directpos]))
    for i in rectangle_side:
        mask = (pos[:,directpos] > i-steppos) & (pos[:,directpos] <= i)
        vel_mean.append(np.mean(avgvel[mask,directvel]))
    

    
    return(np.array(vel_mean)*dx/dt)


def rectangle_velocity_profile_oldandslow(fluid,dx,dt,direction="zx",stepx=1,stepy=1,stepz=1):
    direct = 0 #Direction of the final veloctiy you want
    directone = 1 #The direction which you first take plane slices
    directtwo = 2 #The baseaxis you want compare velocity against
    stepone = 1 #The stepwidth you take in direction one
    steptwo = 1 #The stepwidth you take in direction two

    if direction == "zx":
        direct = 0
        directone = 1
        directtwo = 2
        stepone = stepy
        steptwo = stepz

    if direction == "yx":
        direct = 0
        directone = 2
        directtwo = 1
        stepone = stepz
        steptwo = stepy

    if direction == "xy":
        direct = 1
        directone = 2
        directtwo = 0
        stepone = stepz
        steptwo = stepx
        
    if direction == "zy":
        direct = 1
        directone = 0
        directtwo = 2
        stepone = stepx
        steptwo = stepz
        
    if direction == "xz":
        direct = 2
        directone = 1
        directtwo = 0
        stepone = stepy
        steptwo = stepx
        
    if direction == "yz":
        direct = 2
        directone = 0
        directtwo = 1
        stepone = stepx
        steptwo = stepy
        
        
    avgvel = np.array(fluid[0].velocity)
    for i in range(1,len(fluid)):
        avgvel += np.array(fluid[i].velocity)
    avgvel /= float(len(fluid))


    vel_mean = []
    pos = np.array(fluid[0].position)
    
    for onedx in np.arange(np.min(pos[:,directone]),np.max(pos[:,directone])+stepone,stepone):
        #mask to get data in x-z plane
        mask_plane = (pos[:,directone] < onedx +stepone) & ((pos[:,directone] >= onedx))
        vel_plane = []
        for twodx in np.arange(np.min(pos[:,directtwo]),np.max(pos[:,directtwo])+steptwo,steptwo):
            mask_base_axis = (pos[mask_plane,directtwo] < twodx + steptwo) & (pos[mask_plane,directtwo] >= twodx)
            vel_plane.append(np.mean(avgvel[mask_plane][mask_base_axis,direct]))

        vel_mean.append(np.array(vel_plane))
        
    vel_mean = np.mean(vel_mean,axis=0)
    
    return(vel_mean*dx/dt)
    

def Hd_from_Ht(Ht,diameter):
    a = 1.0 + 1.7*np.exp(-0.35*diameter) - 0.6*np.exp(-0.01*diameter)
    
    Hd1 = (np.sqrt(a**2 - 4.0*a*Ht + 4*Ht) + a)/(2.0*(a-1.0))

    Hd2 = (a - np.sqrt(a**2 - 4.0*a*Ht + 4*Ht))/(2.0*(a-1.0))
    
    return(Hd1,Hd2)

def Ht_from_Hd(Hd,diameter):
    
    a = 1.0 + 1.7*np.exp(-0.35*diameter) - 0.6*np.exp(-0.01*diameter)
    Ht = Hd**2 + Hd*(1.0-Hd)*a
    
    return(Ht)
    


    
