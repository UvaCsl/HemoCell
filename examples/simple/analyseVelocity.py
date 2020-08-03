import numpy as np
import h5py
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import sys

class Fluid:
    def __init__(self, path):
        self.path = path
        files = [f for f in listdir(path) if (isfile(join(path, f)))]
        self.fileNames = [f for f in files if 'Fluid' in f]
        self.files = [h5py.File((self.path + f), 'r') for f in self.fileNames]
        self.zdim, self.ydim, self.xdim = self.files[0].attrs.get('subdomainSize')
        self.velocity = [np.array(f['Velocity']) for f in self.files]
        # print(list(self.files[0].items()))

    def getVelocity(self):
        velocities = []
        positions = []
        for i, velocity in enumerate(self.velocity):
            fluidFile = self.files[i]
            relz, rely, relx = fluidFile.attrs.get('relativePosition')
            xWidth = np.shape(velocity)[2]
            yWidth = np.shape(velocity)[1]
            zWidth = np.shape(velocity)[0]

            for x in range(xWidth):
                for y in range(yWidth):
                    for z in range(zWidth):
                        velocities.append(np.array([velocity[z][y][x][0], velocity[z][y][x][1], velocity[z][y][x][2]]))
                        positions.append(np.array([x + relx, y + rely, z + relz]))

        return velocities, positions

def createVelocityGraph(velocities, positions, xdim, ydim, zdim, velocityAxis = 0, positionAxis = 0):
    xpos = []
    ypos = []
    pos = [0, 1, 2]
    pos.remove(positionAxis)

    for i in range(len(positions)):
        if positions[i][pos[0]] == 10.5 and positions[i][pos[1]] == 10.5:
            xpos.append(positions[i][positionAxis])
            ypos.append(velocities[i][velocityAxis])

    plt.figure()
    plt.plot(xpos, ypos)
    plt.show()

if __name__ == '__main__':
    path = './tmp/hdf5/000000000500/'
    fluid = Fluid(path)
    velocities, positions = fluid.getVelocity()
    createVelocityGraph(velocities, positions, fluid.xdim, fluid.ydim, fluid.zdim, 2, 2)
    # fnameString = './tmp/hdf5/000000000500/Fluid.000000000500.p.0.h5'
    # file = h5File(fnameString)
    # file.keyToDict('Velocity')
    # print(file.dict)
