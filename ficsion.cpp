/* This file is part of the Palabos library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "ficsionInit.hh"
#include "cellsInit.hh"
#include "immersedCells3D.hh"
#include "immersedCellsFunctional3D.hh"
#include "immersedCellsReductions.hh"


using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

plint borderWidth     = 1;  // Because Guo acts in a one-cell layer.
// Requirement: margin>=borderWidth.
plint extraLayer      = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint extendedEnvelopeWidth = 2;  // Because Guo needs 2-cell neighbor access.
const plint particleEnvelopeWidth = 2;

void readFicsionXML(XMLreader document,T & shellDensity, T & k_rest, T & k_stretch, T & k_shear, T & k_bend,
        T & u, T & Re, plint & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & shape, T & radius,
        plint & tmax, plint & tmeas, plint & npar )
    {
    document["cell"]["shellDensity"].read(shellDensity);
    document["cell"]["k_rest"].read(k_rest);
    document["cell"]["k_stretch"].read(k_stretch);
    document["cell"]["k_shear"].read(k_shear);
    document["cell"]["k_bend"].read(k_bend);
    document["parameters"]["u"].read(u);
    document["parameters"]["Re"].read(Re);
    document["parameters"]["N"].read(N);
    document["parameters"]["lx"].read(lx);
    document["parameters"]["ly"].read(ly);
    document["parameters"]["lz"].read(lz);
    document["ibm"]["forceToFluid"].read(forceToFluid);
    document["ibm"]["shape"].read(shape);
    document["ibm"]["radius"].read(radius);
    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    document["sim"]["npar"].read(npar);
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().setStlFilesHaveLowerBound(true);
    global::IOpolicy().setLowerBoundForStlFiles(-1.);

    plint N;
    plint forceToFluid, shape;
    plint tmax, tmeas, npar;
    T u, Re;
    T lx, ly, lz;
    T shellDensity, k_rest, k_stretch, k_shear, k_bend;
    T radius;
    T dt; dt = 0;


    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "reading.." <<std::endl;
    readFicsionXML(document, shellDensity, k_rest, k_stretch, k_shear, k_bend,
            u, Re, N, lx, ly, lz,  forceToFluid, shape, radius, tmax, tmeas, npar);

    IncomprFlowParam<T> parameters(
            u, // u
            Re, // Re
            N,   // N
            lx,        // lx
            ly,        // ly
            lz         // lz
    );

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();
    T tau = parameters.getTau();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz,
                                                            extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;
    iniLattice(lattice, parameters, *boundaryCondition);
    MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
	MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );
    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > immersedParticles (
            particleManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    Box3D inlet(0, 3, 0, ny-1, 0, nz-1);
    Box3D outlet(nx-2, nx-1, 0, ny-1, 0, nz-1);


    std::vector<Array<T,3> > centers;
    std::vector<T> radii;
    positionCells(shape, radius, npar, parameters, centers, radii);

//  === Create Mesh, particles and CellModel ===
    plint numOfCellsPerInlet = radii.size(); // number used for the generation of Cells at inlet
    std::vector<plint> cellIds;
    plint numPartsPerCell = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    TriangleBoundary3D<T> Cells = createCompleteMesh(centers, radii, cellIds, numPartsPerCell, parameters, shape, 200);

    //  plint nProcessors = 1 ;
    //  plint nProcessors = MPI::COMM_WORLD.Get_size() ;
    //  plint LUs=(1+nx)*(1+ny)*(1+nz);
    // plint nTriangles = Cells.getMesh().getNumTriangles();


    pcout << "Mesh Created" << std::endl;
	generateCells(immersedParticles, immersedParticles.getBoundingBox(), cellIds, Cells, numPartsPerCell, numOfCellsPerInlet, slice);

    std::vector<plint> numParts(cellIds.size()); // Count number of particles per Cell
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
            numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), cellIds[iA]);
            pcout << "Cell: " << iA << ", Particles: " << numParts[iA] << std::endl;
	}
    plint totParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox()); //Total number of particles

    std::vector<T> cellsVolume; // Count Volume per Cell
    std::vector<T> cellsSurface; // Count Volume per Cell
    std::vector<T> cellsMeanEdgeDistance;
    std::vector<T> cellsMeanAngle, cellsMeanTriangleArea;
    countCellVolume(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsVolume);
    countCellMeanTriangleArea(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanTriangleArea);
    countCellMeanAngle(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanAngle);
    countCellMeanEdgeDistance(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanEdgeDistance);
    for (pluint iA = 0; iA < cellsVolume.size(); ++iA) {
            pcout << "Cell: " << cellIds[iA]
                  << ", Volume: " << cellsVolume[iA]
                  << ", Mean Triangle Surface: " << cellsMeanTriangleArea[iA]
                  << ", Mean Edge Distance : " << cellsMeanEdgeDistance[iA]
                  << ", Mean Angle: " << cellsMeanAngle[iA]
                  << std::endl;
    }

    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;
    particleArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);

    CellModel3D<T> cellModel(shellDensity, k_rest, k_stretch, k_shear, k_bend, k_shear, k_shear);

    pcout << std::endl << "Starting simulation" << std::endl;
    global::timer("sim").start();
    applyProcessingFunctional ( // copy fluid velocity on particles
        new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(),
        immersedParticles.getBoundingBox(), particleLatticeArg);
    pcout << "Timer; iteration; LU; Cells; Vertices; Triangles; Processors; dt" << std::endl;
    /* ********************* Main Loop ***************************************** * */
    for (plint i=0; i<tmax; ++i) {

         applyProcessingFunctional ( // compute force applied on the particles by springs
             new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
                 Cells, cellModel.clone() ), // used because pushSelect is not used
             immersedParticles.getBoundingBox(), particleArg );
        if (forceToFluid != 0) { // Force from the Cell dynamics to the Fluid
            setExternalVector( lattice, lattice.getBoundingBox(),
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));
            applyProcessingFunctional ( // compute force applied on the fluid by the particles
                    new ForceToFluid3D<T,DESCRIPTOR> (),
                    immersedParticles.getBoundingBox(), particleLatticeArg );
        }
        lattice.collideAndStream();

        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(),
            immersedParticles.getBoundingBox(), particleLatticeArg);

        applyProcessingFunctional ( // advance particles in time according to a velocity, acceleration, ...
            new AdvanceParticlesFunctional3D<T,DESCRIPTOR>,
            immersedParticles.getBoundingBox(), particleArg );

        deleteCell(immersedParticles, outlet, numParts, cellIds, centers, radii );
        Cells.pushSelect(0,1);
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        Cells.popSelect();

        if (i%tmeas==0) { // Output (or screen information every tmeas
            dt = global::timer("sim").stop();
            plint totParticlesNow = 0;
            totParticlesNow = countParticles(immersedParticles, immersedParticles.getBoundingBox());
            pcout << i << " totParticles = " << totParticles << " ";
            // PLB_ASSERT(totParticles == totParticlesNow); //Assert if some particles are outside of the domain

            cellsVolume.clear(); cellsMeanTriangleArea.clear(); cellsMeanEdgeDistance.clear();cellsMeanAngle.clear();
            countCellVolume(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsVolume);
            countCellMeanTriangleArea(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanTriangleArea);
            countCellMeanAngle(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanAngle);
            countCellMeanEdgeDistance(Cells, immersedParticles, immersedParticles.getBoundingBox(), cellIds, cellsMeanEdgeDistance);
            for (pluint iA = 0; iA < cellsVolume.size(); ++iA) {
                    pcout << "Cell: " << cellIds[iA]
                          << ", Volume: " << cellsVolume[iA]
                          << ", Mean Triangle Surface: " << cellsMeanTriangleArea[iA]
                          << ", Mean Edge Distance : " << cellsMeanEdgeDistance[iA]
                          << ", Mean Angle: " << cellsMeanAngle[iA]
                          << std::endl;
            }
            Cells.pushSelect(0,1); // Print the coordinates of one particles
            pcout << "x: " << Cells.getMesh().getVertex(0)[0] << " y: " << Cells.getMesh().getVertex(0)[1] <<
                     " z: " << Cells.getMesh().getVertex(0)[2] <<std::endl;
            Cells.popSelect();

//            pcout << "Timer (w/o Output); " << i <<"; " << LUs << "; " << radii.size() << "; " << itotParticles << "; " << nTriangles << "; " <<nProcessors << "; " <<  dt << ";" << std::endl;
//            pcout << "Write Particle VTK. " << std::endl; ;

            std::vector<std::string> force_scalarNames;
            force_scalarNames.push_back("pressure");
            force_scalarNames.push_back("wss");
            std::vector<std::string> velocity_scalarNames;
            std::vector<std::string> force_vectorNames;
            force_vectorNames.push_back("force");
            std::vector<std::string> velocity_vectorNames;
            velocity_vectorNames.push_back("velocity");
            Cells.pushSelect(0,1);
            Cells.getMesh().writeAsciiSTL(global::directories().getOutputDir()+createFileName("Mesh",i,6)+".stl");
            Cells.popSelect();

            // serialize the particle information to write them.
            // a correspondance between the mesh and the particles is made. (Needs rescale)
            bool dynamicMesh = true;
            plint tag = -1; // Take all triangles.
            writeImmersedSurfaceVTK (
                Cells,
                *getParticlePosAndVelocity(immersedParticles),
                velocity_scalarNames, velocity_vectorNames,
                global::directories().getOutputDir()+createFileName("RBC",i,6)+".vtk", dynamicMesh, tag );
            writeVTK(lattice, parameters, i);
            // === Checkpoint ===
            //    parallelIO::save(immersedParticles, "immersedParticles.dat", true);
            //    parallelIO::save(Cells, "Cells.dat", true);
            //    parallelIO::load("immersedParticles.dat", immersedParticles, true);
            //    parallelIO::load("Cells.dat", Cells, true);
            // ==================
            global::timer("sim").restart();
        }
    }
//    MPI_Finalize();
}
