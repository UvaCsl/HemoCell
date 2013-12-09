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
#include <stdio.h>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "ficsionInit.h"
#include "cellQuantities3D.h"
#include "immersedCells3D.h"
#include "immersedCellsReductions.h"
#include "meshToParticleField3D.h"
#include <string>
#include <map>

#include "shapeMemoryModel3D.h"
#include "cellStretching3D.h"
#include "cellModel3D.hh"

#include "cellInShearFlow3D.h"
#include "meshMetrics.h"
#include "cellForceChecking3D.h"


using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

void readFicsionXML(XMLreader documentXML,std::string & caseId, plint & rbcModel, T & shellDensity, T & k_rest,
        T & k_shear, T & k_bend, T & k_stretch, T & k_WLC, T & eqLengthRatio, T & k_rep, T & k_elastic, T & k_volume, T & k_surface, T & eta_m,
        T & rho_p, T & u, plint & flowType, T & Re, T & shearRate, T & stretchForce, std::vector<T> & eulerAngles, T & Re_p, T & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & ibmKernel, plint & shape, std::string & cellPath, T & radius, T & deflationRatio, pluint & relaxationTime,
        plint & minNumOfTriangles, pluint & tmax, plint & tmeas, plint & npar, plint & goAndStop)
    {
    T nu_p, tau, dx;
    T dt, nu_lb;
    T nx, ny, nz;
    XMLreaderProxy document = documentXML["ficsion"];
    document["caseId"].read(caseId);
    document["cellModel"]["rbcModel"].read(rbcModel);
    document["cellModel"]["shellDensity"].read(shellDensity);
    document["cellModel"]["kWLC"].read(k_WLC);
    document["cellModel"]["eqLengthRatio"].read(eqLengthRatio);
    document["cellModel"]["kRep"].read(k_rep);
    document["cellModel"]["kElastic"].read(k_elastic);
    document["cellModel"]["kBend"].read(k_bend);
    document["cellModel"]["kVolume"].read(k_volume);
    document["cellModel"]["kSurface"].read(k_surface);
    document["cellModel"]["etaM"].read(eta_m);
    document["cellModel"]["kRest"].read(k_rest);
    document["cellModel"]["kShear"].read(k_shear);
    document["cellModel"]["kStretch"].read(k_stretch);
    document["parameters"]["flowType"].read(flowType);
    document["parameters"]["Re"].read(Re);
    document["parameters"]["shearRate"].read(shearRate);
    document["parameters"]["stretchForce"].read(stretchForce); // In picoNewton
    stretchForce *= 1e-12;
    document["parameters"]["eulerAngles"].read(eulerAngles);
    if (eulerAngles.size() != 3) {
        eulerAngles.resize(3, 0.0);
    }
    eulerAngles[0] *= pi/180.;
    eulerAngles[1] *= pi/180.;
    eulerAngles[2] *= pi/180.;
    document["parameters"]["deflationRatio"].read(deflationRatio);
    document["parameters"]["relaxationTime"].read(relaxationTime);
    document["ibm"]["forceToFluid"].read(forceToFluid);
    document["ibm"]["ibmKernel"].read(ibmKernel);
    document["ibm"]["shape"].read(shape);
    if (2 == shape) {
        document["ibm"]["cellPath"].read(cellPath);
    }
    document["ibm"]["radius"].read(radius);
    document["ibm"]["minNumOfTriangles"].read(minNumOfTriangles);
    document["domain"]["rhoP"].read(rho_p);
    document["domain"]["nuP"].read(nu_p);
    document["domain"]["tau"].read(tau);
    document["domain"]["dx"].read(dx);
    document["domain"]["nx"].read(nx);
    document["domain"]["ny"].read(ny);
    document["domain"]["nz"].read(nz);
    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    document["sim"]["npar"].read(npar);

    radius = radius*1.0/dx; // Transform from [m] to [LU]
    nu_lb = (tau-0.5)/3.0;
    dt = (nu_lb/nu_p)*dx*dx;
    u = dt*1.0/dx;
    Re_p = 1.0/nu_p;
    N = int(1.0/dx);
    lx = nx * dx;
    ly = ny * dx;
    lz = nz * dx;
    goAndStop = 0;
    if (flowType == 7) { // Fischer2004 setup, Go-and-stop
        goAndStop = 1;
        flowType = 6;
    }
    if ( (flowType == 4) or (flowType == 5) or (flowType == 8) ) { // Cell Stretching Analysis
        shearRate = 0;
    } else if (flowType == 6) { // Tumbling Tank Treading Measurements
        if (shearRate > 0) {
            tmax  = 100.0/(shearRate*dt); // Measurements for 100 Strain rates
            tmeas = 0.02/(shearRate*dt); // 50 Measurements per Strain rates
        }
    }
    if (goAndStop == 1) {
        tmax += 20.0/dt; // Add 20 seconds to the total time.
    }
}



int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp2/");
    global::directories().setLogOutDir("./tmp2/");
    global::directories().setInputDir("./tmp2/");
    global::IOpolicy().setStlFilesHaveLowerBound(true);
    global::IOpolicy().setLowerBoundForStlFiles(-1.);
//    testInPlane(); PLB_ASSERT(false);

    global::timer("simulation").start();
    std::string performanceLogFileName = global::directories().getLogOutDir() + "performance.log";
    plb_ofstream performanceLogFile(performanceLogFileName.c_str());
    std::string meshQualityFileName = global::directories().getLogOutDir() + "plbMeshQuality.log";
    plb_ofstream meshQualityFile(meshQualityFileName.c_str());

    plint forceToFluid, shape, cellNumTriangles, ibmKernel;
    plint rbcModel;
    std::string caseId;
    std::string cellPath;
    pluint tmax;
    plint tmeas, npar;
    T dtIteration = 0;
    T shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_rep, k_elastic,  k_volume, k_surface, eta_m;
    T eqLengthRatio;
    T u, Re, Re_p, N, lx, ly, lz;
    T rho_p;
    T radius;
    T deflationRatio;
    pluint relaxationTime;
    plint flowType;
    T shearRate, shearRate_p;
    Array<T,3> stretchForce(0,0,0);
    T stretchForceScalar, stretchForce_p;
    std::vector<T> eulerAngles;
    plint goAndStop;

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "reading.." <<std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
            k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface, eta_m,
            rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz,  forceToFluid, ibmKernel, shape, cellPath, radius, deflationRatio, relaxationTime,
            cellNumTriangles, tmax, tmeas, npar, goAndStop);
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
    dNewton = (dm*dx/(dt*dt)) ;
//    kBT = kBT_p / ( dNewton * dx );
    kBT = kBT_p / ( dm * dx*dx/(dt*dt) );
    shearRate = shearRate_p * dt;
    stretchForceScalar = stretchForce_p / dNewton;
    writeFicsionLogFile(parameters, "log", Re, shearRate, flowType);
    pcout << "dx = " << dx << ", " <<
             "dt = " << dt << ", " <<
             "dm = " << dt << ", " <<
             "kT = " << kBT <<
             std::endl;

    /* LATTICE INTIALIZATIONS */
    plint extendedEnvelopeWidth = 2;  // Because Guo needs 2-cell neighbor access.
    plint particleEnvelopeWidth = 2;
    if ((flowType == 0) or (flowType == 1)) {
        extendedEnvelopeWidth = 4;
        particleEnvelopeWidth = 4;
    }

    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;
    if (flowType == 0) { iniLatticeSquarePoiseuille(lattice, parameters, *boundaryCondition, Re); }
    else { iniLatticeSquareCouette(lattice, parameters, *boundaryCondition, shearRate); }
    MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
	MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );
    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > immersedParticles (
            particleManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );
    lattice.periodicity().toggleAll(true);
    immersedParticles.periodicity().toggle(0, true);
    Box3D inlet(0, 3, 0, ny-1, 0, nz-1);
    Box3D outlet(nx-2, nx-1, 0, ny-1, 0, nz-1);
    /* Measure Cell Variables */
    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;
    particleArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);



//  === Create Mesh, particles and CellModel ===
    std::vector<Array<T,3> > centers;
    std::vector<T> radii;
    positionCells(shape, radius, npar, parameters, centers, radii, flowType);
    plint numOfCellsPerInlet = radii.size(); // number used for the generation of Cells at inlet
    std::vector<plint> cellIds;
    plint cellNumVertices = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    std::vector<T> eulerAnglesNONE(3);
    if (flowType == 9) {
        eulerAnglesNONE[0] = 0.;
        eulerAnglesNONE[1] = 0.;
        eulerAnglesNONE[2] = 0.;
        flowType = 1;
    } else {
        eulerAnglesNONE[0] = eulerAngles[0];
        eulerAnglesNONE[1] = eulerAngles[1];
        eulerAnglesNONE[2] = eulerAngles[2];
    }
    TriangleBoundary3D<T> DeCells = createCompleteMesh(centers, radii, eulerAngles, cellIds, parameters, shape, cellPath, cellNumTriangles, cellNumVertices);
    cellIds.clear();
    TriangleBoundary3D<T> Cells = createCompleteMesh(centers, radii, eulerAnglesNONE, cellIds, parameters, shape, cellPath, cellNumTriangles, cellNumVertices);
    pcout << "Mesh Created" << std::endl;
    generateCells(immersedParticles, immersedParticles.getBoundingBox(), cellIds, Cells, cellNumVertices, numOfCellsPerInlet, slice);
    MeshMetrics<T> meshmetric(Cells);
    meshmetric.setResultFile(meshQualityFile);
    meshmetric.write(meshQualityFile);
    writeMeshAsciiSTL(Cells, global::directories().getOutputDir()+createFileName("Mesh",0,8)+".stl");

    performanceLogFile << "# Nx*Ny*Nz; " << nx * ny * nz << std::endl;
    performanceLogFile << "# Nparticles; " << Cells.getMesh().getNumVertices() << std::endl;
    performanceLogFile << "# Nx; Ny; Nz; " << nx << "; " << ny << "; "<< nz << std::endl;
    performanceLogFile << "# Ncells; " << centers.size() << std::endl;

    std::vector<plint> numParts(cellIds.size()); // Count number of particles per Cell
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
            numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), cellIds[iA]);
            pcout << "Cell: " << iA << ", Particles: " << numParts[iA] << std::endl;
    }
    CellFieldQuantityHolder<T,DESCRIPTOR> cqh(npar, numParts[0]);
    // plint totParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox()); //Total number of particles

    /* Find particles in the domain (+envelopes) and store them in tagToParticle3D */
    std::map<plint, Particle3D<T,DESCRIPTOR>*> tagToParticle3D;
    applyProcessingFunctional ( // advance particles in time according to velocity
        new AdvanceParticlesEveryWhereFunctional3D<T,DESCRIPTOR>,
        immersedParticles.getBoundingBox(), particleArg );
    applyProcessingFunctional ( // update mesh position
        new CopyParticleToMeshVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
        immersedParticles.getBoundingBox(), particleArg);
    applyProcessingFunctional (
        new MapVertexToParticle3D<T,DESCRIPTOR> (
            Cells, tagToParticle3D),
        immersedParticles.getBoundingBox(), particleArg );

    std::string logFileName = global::directories().getLogOutDir() + "plbCells.log";
    CellQuantities3D<T,DESCRIPTOR, DenseParticleField3D> rbcQuantities(Cells, immersedParticles, cellIds, npar, tagToParticle3D, logFileName, dx, dt);
    rbcQuantities.calculateAll();


    /* INITIALIZE MODELS */
    k_WLC *= 1.0;     k_rep *= 1.0;     k_elastic *= 1.0;     k_bend *= 1.0;
    k_volume *= 1.0;     k_surface *= 1.0;     k_shear *= 1.0;
    eta_m /= dNewton*dt/dx;     k_stretch /= dNewton;    k_rest /= dNewton/dx;
    T eqArea, eqLength, eqAngle, eqVolume, eqSurface, eqTileSpan;
    eqArea = rbcQuantities.getCellsMeanTriangleArea()[0];     eqLength = rbcQuantities.getCellsMeanEdgeDistance()[0];
    eqAngle = rbcQuantities.getCellsMeanAngle()[0];     eqVolume = rbcQuantities.getCellsVolume()[0];
    eqSurface = rbcQuantities.getCellsSurface()[0]; eqTileSpan = rbcQuantities.getCellsMeanTileSpan()[0];
    T persistenceLengthFine = 7.5e-9  / dx;
    // T eqLengthRatio = 3.17; // According to Pivkin2008
    T maxLength = eqLengthRatio*eqLength;
    pcout << "eqLengthRatio:" << eqLengthRatio  << ", maxLength [LU]:" << maxLength << std::endl;
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    // PLB_PRECONDITION( maxLength < 2.0 );

    vector<T> eqAreaPerTriangle(cellNumTriangles);
    map<plint,T> eqLengthPerEdge, eqAnglePerEdge;
    RBCModel3D<T> *cellModel;
    getCellShapeQuantitiesFromMesh(DeCells, eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumVertices);
    if (rbcModel == 2) {
        rbcModel = 0;
        eqAngle = acos( (sqrt(3.)*(cellNumVertices-2.0) - 5*pi)/(sqrt(3.)*(cellNumVertices-2.0) - 3*pi) );
        map<plint,T>::reverse_iterator iter = eqAnglePerEdge.rbegin();
        for (iter = eqAnglePerEdge.rbegin(); iter != eqAnglePerEdge.rend(); ++iter) {
            eqAnglePerEdge[iter->first] = eqAngle;
        }
    }
    if (rbcModel == 0) {
       cellModel = new ShapeMemoryModel3D<T>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m, \
                       eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, eqVolume, eqSurface, eqTileSpan,
                       persistenceLengthFine, eqLengthRatio, cellNumTriangles, cellNumVertices);
    } else  { // if (rbcModel == 1) {
       cellModel = new CellModel3D<T>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m, \
                       eqArea, eqLength, eqAngle, eqVolume, eqSurface, eqTileSpan,
                       persistenceLengthFine, eqLengthRatio, cellNumTriangles, cellNumVertices);
    }
    pcout << "mu_0 = " << cellModel->getMembraneShearModulus()*dNewton/dx << std::endl;
    pcout << "K = " << cellModel->getMembraneElasticAreaCompressionModulus()*dNewton/dx << std::endl;
    pcout << "YoungsModulus = " << cellModel->getYoungsModulus()*dNewton/dx << std::endl;
    pcout << "Poisson ratio = " << cellModel->getPoissonRatio() << std::endl;


    pcout << std::endl << "Starting simulation" << std::endl;
    pcout << "Timer; iteration; LU; Cells; Vertices; Triangles; Processors; dt" << std::endl;

    /* Deflate if I say so */
    T eqVolumeInitial, eqVolumeFinal;
    eqVolumeInitial = eqVolume;
    eqVolumeFinal = deflationRatio * eqVolumeInitial ;
    applyProcessingFunctional (
       new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (Cells, cellModel->clone(), rbcQuantities.getCellsVolume(), rbcQuantities.getCellsSurface()),
       immersedParticles.getBoundingBox(), particleArg );


    plint timesToStretch = 40;
    CellStretching3D<T,DESCRIPTOR, DenseParticleField3D>  rbcStretch(Cells, immersedParticles, 0.05*numParts[0], flowType, dx, dt, dNewton,
                                                                                            &tagToParticle3D, stretchForceScalar, timesToStretch);
    SingleCellInShearFlow<T,DESCRIPTOR, DenseParticleField3D> shearFlow(Cells, immersedParticles, cellIds,
                                    rbcQuantities.getCellsCenter(), rbcQuantities.getCellsVolume(), 0.05*numParts[0], flowType, dx, dt, T(dNewton), &tagToParticle3D);
    shearFlow.writeHeader( (flowType==6) );
    rbcQuantities.print(0, eqVolume, eqSurface, eqArea, eqLength);
    rbcQuantities.write(0, eqVolumeFinal, eqSurface, eqArea, eqLength) ; // Write Log file for the cell Particles
    pcout << "== Inertia Tensor == "<< std::endl;
    pcout <<  shearFlow.getInertiaTensor()[0][0] << "\t" <<  shearFlow.getInertiaTensor()[0][1] << "\t" <<  shearFlow.getInertiaTensor()[0][2] << std::endl;
    pcout <<  shearFlow.getInertiaTensor()[0][3] << "\t" <<  shearFlow.getInertiaTensor()[0][4] << "\t" <<  shearFlow.getInertiaTensor()[0][5] << std::endl;
    pcout <<  shearFlow.getInertiaTensor()[0][6] << "\t" <<  shearFlow.getInertiaTensor()[0][7] << "\t" <<  shearFlow.getInertiaTensor()[0][8] << std::endl;
    pcout << "== Cell Diameters ==" << std::endl;
    pcout <<  shearFlow.getDiameters()[0][0] << "\t" <<  shearFlow.getDiameters()[0][1] << "\t" <<  shearFlow.getDiameters()[0][2] << std::endl;
    pcout << "== Cell Angles ==" << std::endl;
    pcout <<  shearFlow.getTumblingAngles()[0][0]*180/pi << "\t" <<  shearFlow.getTumblingAngles()[0][1]*180/pi << "\t" <<  shearFlow.getTumblingAngles()[0][2]*180/pi << std::endl;

    bool checkpointed=0;
    /* ********************* Load CheckPoint ***************************************** * */
    if (checkpointed) {
        pcout << "Loading checkpoint... " ;
        parallelIO::load("CHLattice", lattice, true);
        parallelIO::load("CHImmersedParticles", immersedParticles, true);
        rbcQuantities.calculateVolumeAndSurfaceAndCenters();
        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(ibmKernel),
            immersedParticles.getBoundingBox(), particleLatticeArg);
        applyProcessingFunctional ( // Compute Cell Model forces
           new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (Cells, cellModel->clone(), rbcQuantities.getCellsVolume(), rbcQuantities.getCellsSurface()),
           immersedParticles.getBoundingBox(), particleArg );
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToMeshVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        applyProcessingFunctional (
            new MapVertexToParticle3D<T,DESCRIPTOR> (
                Cells, tagToParticle3D),
            immersedParticles.getBoundingBox(), particleArg );
        rbcStretch.applyForce(0, cellModel->getDensity());

        pcout << "OK" << std::endl;
    }

    /* ********************* Main Loop ***************************************** * */
    dtIteration = global::timer("simulation").stop();
    pcout << "Time to initialize: " << dtIteration <<std::endl;
    performanceLogFile << "Init" << "; " << 0 << "; "<< dtIteration << std::endl;
    global::timer("mainLoop").start();
    /* =============================== Main Loop ===================================*/
    for (pluint i=0; i<tmax+1; ++i) {
        if (goAndStop != 0 && i == pluint(tmax - 20/dt)) { // If stopAndGo experiment, turn off the shear
            pcout << "[FICSION]: Switching off shear rate (Fischer2004)" << std::endl;
            shearRate = 0.0; tmeas = 0.1/dt;
            changeCouetteShearRate(lattice, parameters, *boundaryCondition, shearRate);
        }
        /* =============================== OUTPUT ===================================*/
        if (i%tmeas==0) {
            // Stop mainLoop timers and start output timers
            global::timer("mainLoop").stop();
            global::timer("output").restart();
            if (i%(10*tmeas)==0) {
                pcout << "Saving checkpoint... " ;
                // === Checkpoint ===
                std::string fromFileName = global::directories().getOutputDir() + "CHLattice.dat";
                std::string toFileName = global::directories().getOutputDir() + "CHLattice_old.dat";
                if (rename(fromFileName.c_str(), toFileName.c_str()) != 0) { pcout << " (error renaming file) "; }
                fromFileName = global::directories().getOutputDir() + "CHLattice.plb";
                toFileName = global::directories().getOutputDir() + "CHLattice_old.plb";
                if (rename(fromFileName.c_str(), toFileName.c_str()) != 0) { pcout << " (error renaming file) "; }
                fromFileName = global::directories().getOutputDir() + "CHImmersedParticles.dat";
                toFileName = global::directories().getOutputDir() + "CHImmersedParticles_old.dat";
                if (rename(fromFileName.c_str(), toFileName.c_str()) != 0) { pcout << " (error renaming file) "; }
                fromFileName = global::directories().getOutputDir() + "CHImmersedParticles.plb";
                toFileName = global::directories().getOutputDir() + "CHImmersedParticles_old.plb";
                if (rename(fromFileName.c_str(), toFileName.c_str()) != 0) { pcout << " (error renaming file) "; }
                parallelIO::save(lattice, "CHLattice", true);
                parallelIO::save(immersedParticles, "CHImmersedParticles", true);
//                writeMeshAsciiSTL(Cells, global::directories().getOutputDir()+"CH_Mesh.stl");
                pcout << "OK" << std::endl;

            }
            rbcQuantities.write(i, eqVolumeFinal, eqSurface, eqArea, eqLength) ;
            rbcStretch.write(i, rbcQuantities.getCellsMeanEdgeDistance()[0], rbcQuantities.getCellsMaxEdgeDistance()[0]); // SINGLE CELL STRETCHING
            if (flowType==6) { // SINGLE CELL IN SHEAR FLOW
                shearFlow.updateQuantities(i, cellIds, rbcQuantities.getCellsCenter(), rbcQuantities.getCellsVolume());
                shearFlow.write();
                pcout << "[SF] Diameters :" <<  shearFlow.getDiameters()[0][0] << "\t" <<  shearFlow.getDiameters()[0][1] << "\t" <<  shearFlow.getDiameters()[0][2] ;
                pcout << ", Angles :" <<  shearFlow.getTumblingAngles()[0][0]*180/pi << "\t" <<  shearFlow.getTumblingAngles()[0][1]*180/pi << "\t" <<  shearFlow.getTumblingAngles()[0][2]*180/pi << std::endl;
            }
            writeImmersedSurfaceVTK (
                Cells, immersedParticles,
                global::directories().getOutputDir()+createFileName("RBC",i,8)+".vtk");
            writeVTK(lattice, parameters, i);
            // WRITE PERFORMANCE OUTPUT
            dtIteration = global::timer("mainLoop").stop(); global::timer("mainLoop").reset();
            if (i>0) { performanceLogFile << "Iteration" << "; " << i-tmeas << "; "<< dtIteration*1.0/tmeas << std::endl; }
            pcout << "Iteration " << i << ", time "<< dtIteration*1.0/tmeas << std::endl;
            dtIteration = global::timer("LBM").stop(); global::timer("LBM").reset();
            if (i>0) { performanceLogFile << "LBM" << "; " << i << "; "<< dtIteration*1.0/tmeas  << std::endl; }
            dtIteration = global::timer("IBM").stop(); global::timer("IBM").reset();
            if (i>0) { performanceLogFile << "IBM" << "; " << i << "; "<< dtIteration*1.0/tmeas  << std::endl; }
            dtIteration = global::timer("Quantities").stop(); global::timer("Quantities").reset();
            if (i>0) { performanceLogFile << "Quantities" << "; " << i << "; "<< dtIteration*1.0/tmeas  << std::endl; }
            dtIteration = global::timer("Model").stop(); global::timer("Model").reset();
            if (i>0) { performanceLogFile << "Model" << "; " << i << "; "<< dtIteration*1.0/tmeas  << std::endl; }
            dtIteration = global::timer("output").stop(); global::timer("output").reset();
            performanceLogFile << "Output" << "; " << i << "; "<< dtIteration << std::endl;
            global::timer("mainLoop").start();
        }
        /* ====================== Update for Cell Stretching behaviour ===================================*/
        rbcStretch.hasConverged(i);


        /* =============================== MAIN ACTIONS ===================================*/
        // #0# Equilibration
        if (i<=relaxationTime) {
            eqVolume = eqVolumeInitial + (i*1.0)/relaxationTime * (eqVolumeFinal - eqVolumeInitial) ;
            cellModel->setEquilibriumVolume(eqVolume);
            if (flowType == 3) {
                rbcStretch.setStretchScalarForce( (i*1.0)/relaxationTime * stretchForceScalar );
            }
        }


        // #2# IBM Spreading
        global::timer("IBM").start();
        if (forceToFluid != 0) { // Force from the Cell dynamics to the Fluid
            setExternalVector( lattice, lattice.getBoundingBox(),
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));
            applyProcessingFunctional ( // compute force applied on the fluid by the particles
                    new ForceToFluid3D<T,DESCRIPTOR> (ibmKernel),
                    immersedParticles.getBoundingBox(), particleLatticeArg );
        }
        global::timer("IBM").stop();


        // #3# LBM
        global::timer("LBM").start();
        lattice.collideAndStream();
        dtIteration = global::timer("LBM").stop();

        global::timer("IBM").start();
        // #4# IBM Interpolation
        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(ibmKernel),
            immersedParticles.getBoundingBox(), particleLatticeArg);
        // #5# Position Update
        applyProcessingFunctional ( // advance particles in time according to velocity
            new AdvanceParticlesEveryWhereFunctional3D<T,DESCRIPTOR>,
            immersedParticles.getBoundingBox(), particleArg );
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToMeshVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        dtIteration = global::timer("IBM").stop();

//        applyProcessingFunctional (
//            new MeshToParticleField3D<T,DESCRIPTOR> (cqh),
//            immersedParticles.getBoundingBox(), particleArg );
//        std::map<plint, plint> cid2mcid = cqh.getCellIdToMeshCellId();
//        std::map<plint, plint>::iterator iter = cid2mcid.begin();
//    	for (; iter != cid2mcid.end();) {
//    		std::cout << "cid2mcid " << i << " " << cid2mcid.size() << " cid " << iter->first << ", mid " << iter->second << std::endl;
//    		iter++;
//    	}
        global::timer("Quantities").start();
        applyProcessingFunctional (
            new MapVertexToParticle3D<T,DESCRIPTOR> (
                Cells, tagToParticle3D),
            immersedParticles.getBoundingBox(), particleArg );
//        pcout << "tagToParticle3D: " << tagToParticle3D.size() << std::endl;
        if (flowType==6) { rbcQuantities.calculateVolumeAndSurfaceAndCenters(); }
        else { rbcQuantities.calculateVolumeAndSurface();  }
        dtIteration = global::timer("Quantities").stop();

        global::timer("Model").start();
        // #1# Membrane Model + Stretching
        applyProcessingFunctional (
                        new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (Cells, cellModel->clone(), rbcQuantities.getCellsVolume(), rbcQuantities.getCellsSurface()),
                        immersedParticles.getBoundingBox(), particleArg );
        rbcStretch.applyForce(i, cellModel->getDensity());
        dtIteration = global::timer("Model").stop();
    }
    pcout << "Simulation finished." << std::endl;
//    MPI_Finalize();
}
