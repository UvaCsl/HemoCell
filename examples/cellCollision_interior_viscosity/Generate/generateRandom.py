
import random as rand

# Quick script to go over all the cellCollision folders and change the initial orientations of the platelet
# And rbc's

d = 2.0  # distance between y-coordinates of particles
shear= 4000.0
visco=1

center=12.0
folderRoot = '$HOME/cellCollider/'+'gam{0:.0f}_d{1:.0f}_visco{2}'.format(shear, d, visco)+'/'

tn = d/2
iters = 500000

# For setting up the LISA script
nodes = 2
cores = 16  # Amount of simulations one wants to run per node

for i in range(cores*nodes):
    folderPath = 'cellCollision' + str(i+1)
    print folderPath


    with open(folderPath+'/RBC.pos', 'wb') as f:
        f.write('1\n')

        f.write('5.5 {0:.1f} 7.0 {1} {2} {3}'.format(center+tn, rand.randint(0,359), rand.randint(0,359), rand.randint(0,359)))


    with open(folderPath+'/PLT.pos', 'wb') as f:
        f.write('1\n')
        f.write('21.5 {0:.1f} 7.0 {1} {2} {3}'.format(center-tn, rand.randint(0,359), rand.randint(0,359), rand.randint(0,359)))


    with open(folderPath+'/config.xml', 'wb') as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<hemocell>\n')
        f.write('<parameters>\n')
        f.write('\t <warmup> 0 </warmup> <!-- Number of LBM iterations to prepare fluid field. -->\n')
        f.write('\t <enableInteriorViscosity> 1 </enableInteriorViscosity>\n')
        f.write('</parameters>\n')
        f.write('<ibm>\n')

        f.write('\t <radius> 3.91e-6 </radius> <!-- Radius of the particle in [m] (dx) [3.3e-6, 3.91e-6, XX and 4.284 for shapes [0,1,2,3] respectively -->\n')
        f.write('\t<minNumOfTriangles> 600 </minNumOfTriangles> <!--Minimun numbers of triangles per cell. Not always exact.-->\n')
        f.write('</ibm>\n')

        f.write('\n')
        f.write('<domain>\n')
        f.write('\t<shearrate> {0:.1f} </shearrate>   <!--Shear rate for the fluid domain. [s^-1] [25]. -->\n'.format(shear))
        f.write('\t <rhoP> 1025 </rhoP>   <!--Density of the surrounding fluid, Physical units [kg/m^3]-->\n')
        f.write('\t <nuP> 1.1e-6 </nuP>   <!-- Kinematic viscosity of the surrounding fluid, physical units [m^2/s]-->\n')
        f.write('\t <dx> 0.5e-6 </dx> <!--Physical length of 1 Lattice Unit -->\n')
        f.write('\t <dt> 1e-7 </dt> <!-- Time step for the LBM system. A negative value will set Tau=1 and calc. the corresponding time-step. -->\n')
        f.write('\t <timeStepSize> 1 </timeStepSize> <!-- Update particle material model after how many fluid time steps. [Integer] -->\n')
        f.write('\t <particleEnvelope>20</particleEnvelope>\n')
        f.write('\t <kBT>4.100531391e-21</kBT> <!-- in SI, m2 kg s-2 (or J) for T=300 -->\n')
        f.write('</domain>\n')
        f.write('\n')

        f.write('<sim>\n')
        f.write('\t <tmax> {0} </tmax> <!-- total number of iterations -->\n'.format(iters))
        f.write('\t <tmeas> 500 </tmeas> <!-- interval after which data is written -->\n') 
        f.write('\t <tcheckpoint>1000000</tcheckpoint>\n')
        f.write('\t <interiorViscosity>' + str(visco) + '</interiorViscosity>\n')
	f.write('\t <interorViscosityEntireGrid>200</interiorViscosityEntireGrid>\n')
        f.write('</sim>\n')

        f.write('</hemocell>\n')

## Generate a LISA jobscript

with open('gam{0:.0f}_d{1:.0f}_visco{2}'.format(shear, d, visco), 'wb') as f:
    f.write('#PBS -S /bin/bash\n')
    f.write('#\n')
    f.write('#PBS -lnodes={0}:cores{1}\n'.format(nodes, cores))
    f.write('#PBS -lwalltime=10:00:00\n')

    f.write('module load hdf5/intel/1.8.16-parallel\n')
    f.write('module load openmpi/gnu\n')
    f.write('module list\n')


    for i in range(cores*nodes):
        folderPath = 'cellCollision' + str(i+1)
	# Root of the folder with all cellCollisionFolders
        f.write('(\n')
        f.write('cd '+folderRoot+folderPath+'\n')
        f.write('echo "Start the program in folder {0}"\n'.format(i+1))
        f.write('./cellCollision config.xml\n')
        f.write(') &\n')
        f.write('\n')

    f.write('wait')
