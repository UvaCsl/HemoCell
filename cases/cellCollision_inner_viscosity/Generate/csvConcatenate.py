import csv
import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import os
# import sys

  # Set the cell type names to look for
cellTypes = ["RBC_HO", "PLT"]
  # Set the last iteration where we had all cells
maxIter = 700000
  # Shear rate
shear = 200
  # dt
dt = 1e-7

folder = "tmp" # Folder where the results are stored


def concatenate():
    """
    For ease of working with the output csv, we output all the seperate csv files
    to a single csv per particle type which also contains the timesteps
    """

    # First create the relevant data structure containing strings for timesteps
    files = {}
    for i in cellTypes:
        files[i]={}

    for root, dirnames, filenames in sorted(os.walk(folder+'/csv')):
        # print type(root), type(filenames), type(dirnames)
        for filename in fnmatch.filter(filenames, '*.csv'):
            for i in cellTypes:
                if i in filename:
                    timestep = root.split('/')
                    if int(timestep[-1]) > maxIter:
                        continue
                    files[i][timestep[-1]] = ''

    # Go over the files now and actually add to our dicts
    for root, dirnames, filenames in sorted(os.walk(folder+'/csv')):
        for filename in fnmatch.filter(filenames, '*.csv'):
            for i in cellTypes:
                if i in filename:
                    timestep = root.split('/')[-1]
                    if int(timestep) > maxIter:
                        continue
                    csvFile = root + '/' + filename

                    ifile = open(csvFile, "rb")
                    reader = csv.reader(ifile)
                    reader.next() #header
                    try:
                        temp = ','.join(reader.next())
                        timeReal = '{0:.6f},'.format(int(timestep)*dt)
                        files[i][timestep] += (timeReal)
                        files[i][timestep] += temp

                    except:
                        pass
                    # coords = map(float, reader.next())
                    ifile.close()

    for i in cellTypes:
        with open(i+'_data.csv', 'wb') as csvfile:
            csvWriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)

            # csvWriter.writerow()
            for key in sorted(files[i].keys(), key= lambda ke: int(ke)):
                csvWriter.writerow(files[i][key].split(','))
            
if __name__ == "__main__":
    concatenate()



