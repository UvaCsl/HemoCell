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

void readFicsionXML(XMLreader document,T & shellDensity,
        T & k_shear, T & k_bend, T & k_stretch, T & k_WLC, T & k_elastic, T & k_volume, T & k_surface,
        T & rho_p, T & u, plint & flowType, T & shearRate, T & Re, T & Re_p, T & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & shape, T & radius,
        plint & tmax, plint & tmeas, plint & npar)
    {
    T nu_p, tau, dx;
    T dt, nu_lb;
    T nx, ny, nz;

    document["cell"]["shellDensity"].read(shellDensity);
    document["cell"]["k_WLC"].read(k_WLC);
    document["cell"]["k_elastic"].read(k_elastic);
    document["cell"]["k_shear"].read(k_shear);
    document["cell"]["k_bend"].read(k_bend);
    document["cell"]["k_stretch"].read(k_stretch);
    document["cell"]["k_volume"].read(k_volume);
    document["cell"]["k_surface"].read(k_surface);
    document["parameters"]["Re"].read(Re);
    document["parameters"]["shearRate"].read(shearRate);
    document["parameters"]["flowType"].read(flowType);
    document["parameters"]["rho_p"].read(rho_p);
    document["parameters"]["nu_p"].read(nu_p);
    document["parameters"]["tau"].read(tau);
    document["parameters"]["dx"].read(dx);
    document["parameters"]["nx"].read(nx);
    document["parameters"]["ny"].read(ny);
    document["parameters"]["nz"].read(nz);
    document["ibm"]["forceToFluid"].read(forceToFluid);
    document["ibm"]["shape"].read(shape);
    document["ibm"]["radius"].read(radius);
    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    document["sim"]["npar"].read(npar);

    nu_lb = (tau-0.5)/3.0;
    dt = (nu_lb/nu_p)*dx*dx;
    u = dt*1.0/dx;
    Re_p = 1.0/nu_p;
    N = int(1.0/dx);
    lx = nx * dx;
    ly = ny * dx;
    lz = nz * dx;
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::directories().setLogOutDir("./tmp/");
    global::IOpolicy().setStlFilesHaveLowerBound(true);
    global::IOpolicy().setLowerBoundForStlFiles(-1.);

    std::string logFileName = global::directories().getLogOutDir() + "plbCells.log";
    plb_ofstream logFile(logFileName.c_str());
    logFile << "# Iteration; Volume, Surface, Mean Triangle Surface, Mean Edge Distance, Mean Angle [^o], Mean Edge Distance [LU], Max Edge Distance [LU]" << std::endl;

    plint forceToFluid, shape;
    plint tmax, tmeas, npar;
    T dtIteration = 0;
    T shellDensity, k_shear, k_bend, k_stretch, k_WLC, k_elastic,  k_volume, k_surface, eta_m=0;
    T u, Re, Re_p, N, lx, ly, lz;
    T rho_p;
    T radius;
    plint flowType;
    T shearRate, shearRate_p;

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "reading.." <<std::endl;
    readFicsionXML(document, shellDensity, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface,
            rho_p, u, flowType, shearRate_p, Re, Re_p, N, lx, ly, lz,  forceToFluid, shape, radius, tmax, tmeas, npar);
    IncomprFlowParam<T> parameters(
            u, // u
            Re_p, // Inverse viscosity (1/nu_p)
            N,   // N
            lx,        // lx
            ly,        // ly
            lz         // lz
    );

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();
    T tau = parameters.getTau();
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    T dm = rho_p * (dx*dx*dx);
    kBT = kBT_p / ( dx*dx * dm /(dt*dt) );
    dNewton = (dm*dx/(dt*dt)) ;
    shearRate = shearRate_p * dt;

    writeFicsionLogFile(parameters, "log", Re, shearRate, flowType);
    pcout << "kT = " << kBT << std::endl;

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
    if (flowType == 1) {
        iniLatticeSquareCouette(lattice, parameters, *boundaryCondition, shearRate);
    } else {
        iniLatticeSquarePoiseuille(lattice, parameters, *boundaryCondition, Re);
    }
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
    positionCells(shape, radius, npar, parameters, centers, radii, flowType);

//  === Create Mesh, particles and CellModel ===
    plint numOfCellsPerInlet = radii.size(); // number used for the generation of Cells at inlet
    std::vector<plint> cellIds;
    plint numPartsPerCell = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    TriangleBoundary3D<T> Cells = createCompleteMesh(centers, radii, cellIds, numPartsPerCell, parameters, shape, 500);
    pcout << "Mesh Created" << std::endl;
    generateCells(immersedParticles, immersedParticles.getBoundingBox(), cellIds, Cells, numPartsPerCell, numOfCellsPerInlet, slice);

    //  plint nProcessors = 1 ;
    //  plint nProcessors = MPI::COMM_WORLD.Get_size() ;
    //  plint LUs=(1+nx)*(1+ny)*(1+nz);
    // plint nTriangles = Cells.getMesh().getNumTriangles();

    std::vector<plint> numParts(cellIds.size()); // Count number of particles per Cell
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
            numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), cellIds[iA]);
            pcout << "Cell: " << iA << ", Particles: " << numParts[iA] << std::endl;
    }
    plint totParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox()); //Total number of particles

    std::vector<T> cellsVolume, cellsSurface;
    std::vector<T> cellsMeanEdgeDistance, cellsMaxEdgeDistance, cellsMeanAngle, cellsMeanTriangleArea;
    T eqArea, eqLength, eqAngle, eqVolume, eqSurface;
    calculateCellMeasures(Cells, immersedParticles, cellIds, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                        cellsMaxEdgeDistance, cellsMeanAngle);
    eqArea = cellsMeanTriangleArea[0];     eqLength = cellsMeanEdgeDistance[0];
    eqAngle = cellsMeanAngle[0];     eqVolume = cellsVolume[0];
    eqSurface = cellsSurface[0];
    printCellMeasures(0, Cells, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                           cellsMaxEdgeDistance, cellsMeanAngle, eqVolume, eqSurface, eqArea, eqLength) ;

    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;
    particleArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);

    T persistenceLength = 7.5e-9 * sqrt(23865.0/(numParts[0]-2)) / dx;
    T maxLength = 2.2*eqLength;
    pcout << "persistenceLength = " << persistenceLength
            << " maxLength = " << maxLength << std::endl;
    PLB_PRECONDITION( maxLength < 1.0 );
//    eqVolume = pow(eqSurface,1.5)/(6*pi);
    k_WLC *= 1;
    k_elastic *= 1;
    k_bend *= kBT;
    k_volume /= dNewton/(dx*dx);
    k_surface /= dNewton/dx;
    /* == */
    k_shear /= dNewton/dx;
    k_stretch /= dNewton;
    CellModel3D<T> cellModel(shellDensity, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, \
                                                eqArea, eqLength, eqAngle, eqVolume, eqSurface, eta_m,
                                                maxLength, persistenceLength);
    pcout << std::endl << "Starting simulation" << std::endl;
    global::timer("sim").start();
    pcout << "Timer; iteration; LU; Cells; Vertices; Triangles; Processors; dt" << std::endl;
    /* ********************* Main Loop ***************************************** * */
    T eqVolumeInitial, eqVolumeFinal, ifinal=10000.;
    eqVolumeInitial = eqVolume;
    eqVolumeFinal = eqVolume; // 0.65 * eqVolumeInitial ;
//    eqVolumeFinal = 0.65 * eqVolumeInitial ;
    applyProcessingFunctional ( // copy fluid velocity on particles
        new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(),
        immersedParticles.getBoundingBox(), particleLatticeArg);


    for (plint i=0; i<tmax; ++i) {
        if (i<=ifinal) {
            eqVolume = eqVolumeInitial + (i*1.0)/ifinal * (eqVolumeFinal - eqVolumeInitial) ;
            cellModel.setEquilibriumVolume(eqVolume);
        }
         applyProcessingFunctional ( // compute force applied on the particles by springs
             new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
                 Cells, cellModel.clone(), cellsVolume, cellsSurface),
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
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);

        if (i%tmeas==0) { // Output (or screen information every tmeas
            dtIteration = global::timer("sim").stop();
            plint totParticlesNow = 0;
            totParticlesNow = countParticles(immersedParticles, immersedParticles.getBoundingBox());
            pcout << i << " totParticles = " << totParticles << std::endl;
            // PLB_ASSERT(totParticles == totParticlesNow); //Assert if some particles are outside of the domain
            calculateCellMeasures(Cells, immersedParticles, cellIds, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                                cellsMaxEdgeDistance, cellsMeanAngle);
            printCellMeasures(i, Cells, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                                   cellsMaxEdgeDistance, cellsMeanAngle, eqVolume, eqSurface, eqArea, eqLength) ;
            writeCellLog(i, logFile,
                         cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                         cellsMaxEdgeDistance, cellsMeanAngle,
                         eqVolume, eqSurface, eqArea, eqLength) ;

            std::vector<std::string> force_scalarNames;
            force_scalarNames.push_back("pressure");
            force_scalarNames.push_back("wss");
            std::vector<std::string> velocity_scalarNames;
            std::vector<std::string> force_vectorNames;
            force_vectorNames.push_back("force");
            std::vector<std::string> velocity_vectorNames;
            velocity_vectorNames.push_back("velocity");
            writeMeshAsciiSTL(Cells, global::directories().getOutputDir()+createFileName("Mesh",i,6)+".stl");
//            Cells.getMesh().writeAsciiSTL(global::directories().getOutputDir()+createFileName("Mesh",i,6)+".stl");
            // serialize the particle information to write them.
            // a correspondance between the mesh and the particles is made. (Needs rescale)
//            bool dynamicMesh = true;
//            plint tag = -1; // Take all triangles.
//            writeImmersedSurfaceVTK (
//                Cells,
//                *getParticlePosAndVelocity(immersedParticles),
//                velocity_scalarNames, velocity_vectorNames,
//                global::directories().getOutputDir()+createFileName("RBC",i,6)+".vtk", dynamicMesh, tag );
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
