import h5py
import glob
files1 = glob.glob("case1/tmp/hdf5/RBC.*.h5")
files1.sort() 
files20 = glob.glob("case2/tmp/hdf5/RBC.*.000.h5")
files20.sort()
files21 = glob.glob("case2/tmp/hdf5/RBC.*.001.h5")
files21.sort()

for i in range(0,len(files1)):
    f1 = h5py.File(files1[i], "r")
    f20= h5py.File(files20[i], "r")
    f21= h5py.File(files21[i], "r")
    v1 = {}
    v2 = {}
    for j in range(0,len(f1['VertexId'][:])):
        v1[f1['VertexId'][j]] = f1['pbcPosition'][i]
    try:
        for j in range(0,len(f20['VertexId'][:])):
            v2[f20['VertexId'][j]] = f20['pbcPosition'][i]
    except:
        pass
    try:
        for j in range(0,len(f21['VertexId'][:])):
            v2[f21['VertexId'][j]] = f21['pbcPosition'][i]
    except:
        pass

    for v1k, v1v in v1.iteritems():
        if not v2[v1k][0] == v1v[0]:
            print 'Difference for vertex,value1,value2',v1k,v1v,v2[v1k]
            continue
        if not v2[v1k][1] == v1v[1]:
            print 'Difference for vertex,value1,value2',v1k,v1v,v2[v1k]
            continue
        if not v2[v1k][2] == v1v[2]:
            print 'Difference for vertex,value1,value2',v1k,v1v,v2[v1k]
            continue

