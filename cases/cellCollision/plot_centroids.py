import csv
import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import os

  # Set the cell type names to look for
cell_types = ["RBC_HO", "PLT"]
  # Set the last iteration where we had all cells
maxIter = 420000
  # Shear rate
shear = 200
  # dt
dt = 1e-7


def getCoord(csv_file):
    ifile = open(csv_file, "rb")
    reader = csv.reader(ifile)
    reader.next() #header
    coords = map(float, reader.next())
    ifile.close()
    
    return coords


if __name__ == "__main__":
    print "For now, only single result file / save is supported. If you have multiple csv/save, concatenate them first."

    import sys

    folder = "tmp"

    if len(sys.argv) > 1:
        folder = sys.argv[1]
    
    files = {}
    for i in cell_types:
        files[i]=[]
       
        
    for root, dirnames, filenames in sorted(os.walk(folder+'/csv')):
        for filename in fnmatch.filter(filenames, '*.csv'):
            for i in cell_types:
                if i in filename:
                    print filename
                    iteration = int(os.path.basename(os.path.normpath(root)))
                    if iteration < maxIter:
                        files[i].append( [iteration*dt] + getCoord(os.path.join(root, filename)) )


    for i in files:
        files[i] = np.array(files[i])
        dy = files[i][-1, 2] - files[i][0, 2]
        plt.plot(files[i][:, 0], files[i][:, 2], label=i + "  (dy: " + str(dy)+")")  

    plt.xlabel("time [s]")
    plt.ylabel("y position [um]")
    plt.title("Shear rate of " + str(shear) + " s^-1")
    plt.legend(loc=3)            
    plt.savefig("collision.png", dpi=300)
    #plt.show()  
