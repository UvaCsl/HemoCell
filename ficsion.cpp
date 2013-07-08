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
#include <string>
#include <map>

#include "shapeMemoryModel3D.hh"
#include "shapeMemoryModelFunctional3D.hh"
#include "cellStretchingForces3D.hh"
#include "cellModel3D.hh"

#include "cellInShearFlow3D.h"


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

void readFicsionXML(XMLreader documentXML,std::string & caseId, plint & rbcModel, T & shellDensity, T & k_rest,
        T & k_shear, T & k_bend, T & k_stretch, T & k_WLC, T & eqLengthRatio, T & k_rep, T & k_elastic, T & k_volume, T & k_surface, T & eta_m,
        T & rho_p, T & u, plint & flowType, T & Re, T & shearRate, T & stretchForce, std::vector<T> & eulerAngles, T & Re_p, T & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & ibmKernel, plint & shape, std::string & cellPath, T & radius, T & deflationRatio, plint & relaxationTime,
        plint & minNumOfTriangles, pluint & tmax, plint & tmeas, plint & npar, plint & goAndStop, plint & fastStretchRelease)
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
    fastStretchRelease = 0;
    if (flowType == 8) { // Fischer2004 setup, Go-and-stop
        fastStretchRelease = 1;
        flowType = 5;
    }
    if (flowType == 7) { // Fischer2004 setup, Go-and-stop
        goAndStop = 1;
        flowType = 6;
    }
    if ( (flowType == 4) or (flowType == 5) ) { // Cell Stretching Analysis
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
    global::directories().setOutputDir("./tmp/");
    global::directories().setLogOutDir("./tmp/");
    global::IOpolicy().setStlFilesHaveLowerBound(true);
    global::IOpolicy().setLowerBoundForStlFiles(-1.);

    std::string logFileName = global::directories().getLogOutDir() + "plbCells.log";
    plb_ofstream logFile(logFileName.c_str());

    std::string stretchLogFileName = global::directories().getLogOutDir() + "stretchDeformation.log";
    std::string stretchResultFileName = global::directories().getLogOutDir() + "stretchResults.log";
    std::string stretchResealedFileName = global::directories().getLogOutDir() + "stretchRelease.log";
    std::string shearResultFileName = global::directories().getLogOutDir() + "shearResults.log";
    plb_ofstream stretchLogFile(stretchLogFileName.c_str());
    plb_ofstream stretchResultFile(stretchResultFileName.c_str());
    plb_ofstream stretchResealedFile(stretchResealedFileName.c_str());
    plb_ofstream shearResultFile(shearResultFileName.c_str());

    plint forceToFluid, shape, cellNumTriangles, ibmKernel;
    plint rbcModel;
    std::string caseId;
    std::string cellPath;
    pluint tmax;
    plint tmeas, npar;
    // T dtIteration = 0;
    T shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_rep, k_elastic,  k_volume, k_surface, eta_m;
    T eqLengthRatio;
    T u, Re, Re_p, N, lx, ly, lz;
    T rho_p;
    T radius;
    T deflationRatio;
    plint relaxationTime;
    plint flowType;
    T shearRate, shearRate_p;
    Array<T,3> stretchForce(0,0,0);
    T stretchForceScalar, stretchForce_p;
    std::vector<T> eulerAngles;
    plint goAndStop, fastStretchRelease;
    std::vector<plint> goAndStopVertices;
    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "reading.." <<std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
            k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface, eta_m,
            rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz,  forceToFluid, ibmKernel, shape, cellPath, radius, deflationRatio, relaxationTime,
            cellNumTriangles, tmax, tmeas, npar, goAndStop, fastStretchRelease);
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
    if (flowType == 0) {
        iniLatticeSquarePoiseuille(lattice, parameters, *boundaryCondition, Re);
    } else {
        iniLatticeSquareCouette(lattice, parameters, *boundaryCondition, shearRate);
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
    plint cellNumVertices = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    TriangleBoundary3D<T> Cells = createCompleteMesh(centers, radii, eulerAngles, cellIds, parameters, shape, cellPath, cellNumTriangles, cellNumVertices);
    pcout << "Mesh Created" << std::endl;
    generateCells(immersedParticles, immersedParticles.getBoundingBox(), cellIds, Cells, cellNumVertices, numOfCellsPerInlet, slice);

    std::vector<plint> numParts(cellIds.size()); // Count number of particles per Cell
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
            numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), cellIds[iA]);
            pcout << "Cell: " << iA << ", Particles: " << numParts[iA] << std::endl;
    }
    // plint totParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox()); //Total number of particles

    /* Measure Cell Variables */
    std::vector<T> cellsVolume, cellsSurface;
    std::vector<T> cellsMeanEdgeDistance, cellsMaxEdgeDistance, cellsMeanAngle, cellsMeanTriangleArea, cellsMeanTileSpan;
    std::vector< Array<T,3> > cellsCenter, cellsVelocity;
    T eqArea, eqLength, eqAngle, eqVolume, eqSurface, eqTileSpan;
    calculateCellMeasures(Cells, immersedParticles, cellIds, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                        cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity, cellsMeanTileSpan);
    eqArea = cellsMeanTriangleArea[0];     eqLength = cellsMeanEdgeDistance[0];
    eqAngle = cellsMeanAngle[0];     eqVolume = cellsVolume[0];
    eqSurface = cellsSurface[0];	eqTileSpan = cellsMeanTileSpan[0];
    printCellMeasures(0, Cells, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                           cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity, eqVolume, eqSurface, eqArea, eqLength,
                           dx, dt) ;

    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;
    particleArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);

    std::vector<plint> outerLeftTags, outerRightTags;
    std::vector<plint> outerFrontTags, outerBackTags;
    std::vector<T> stretchingDeformations;
    std::vector<std::vector<plint>*> lateralCellParticleTags;
    lateralCellParticleTags.push_back(&outerLeftTags);
    lateralCellParticleTags.push_back(&outerRightTags);
    lateralCellParticleTags.push_back(&outerFrontTags);
    lateralCellParticleTags.push_back(&outerBackTags);
    plint numParticlesPerSide = plint(0.02*numParts[0]);
    if ((flowType == 3) || (flowType == 4) || (flowType == 5) ) {
        PLB_PRECONDITION( npar == 1 && MPI::COMM_WORLD.Get_size() == 1 );
        applyProcessingFunctional (
            new FindTagsOfLateralCellParticles3D<T,DESCRIPTOR>(numParticlesPerSide, &outerLeftTags, &outerRightTags, 0),
            immersedParticles.getBoundingBox(), particleArg );
        applyProcessingFunctional ( 
            new FindTagsOfLateralCellParticles3D<T,DESCRIPTOR>(numParticlesPerSide, &outerFrontTags, &outerBackTags, 2),
            immersedParticles.getBoundingBox(), particleArg );
        stretchLogFile << setprecision(20) << "# Force [N] = " << stretchForce_p << std::endl << "# t [sec]; D_A [m]; D_T [m]; Mean Edge Distance [LU]; Max Edge Distance [LU]; " << std::endl;
    }
    for (pluint ipa = 0; ipa < outerLeftTags.size(); ++ipa) {
        pcout << ipa << " : " <<outerLeftTags[ipa] ;
        pcout << ", " <<outerRightTags[ipa] << std::endl;
    }


    T persistenceLengthFine = 7.5e-9  / dx;
    // T eqLengthRatio = 3.17; // According to Pivkin2008
    T maxLength = eqLengthRatio*eqLength;
    pcout << "eqLengthRatio:" << eqLengthRatio  << ", maxLength [LU]:" << maxLength << std::endl;
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    // PLB_PRECONDITION( maxLength < 2.0 );
    k_WLC *= 1.0;     k_rep *= 1.0;     k_elastic *= 1.0;     k_bend *= 1.0;
    k_volume *= 1.0;     k_surface *= 1.0;     k_shear *= 1.0;
    /* == */
    eta_m /= dNewton*dt/dx;
    k_stretch /= dNewton;
    k_rest /= dNewton/dx;
    vector<T> eqAreaPerTriangle(cellNumTriangles);
    map<plint,T> eqLengthPerEdge, eqAnglePerEdge;
    T eta_mPrevious = 0.;
    if (flowType == 5)  {
        eta_mPrevious = eta_m;
        eta_m = 0.0;
    }
    ShellModel3D<T> *cellModel;
    if (rbcModel == 0) {
       getCellShapeQuantitiesFromMesh(Cells, eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumVertices);
       cellModel = new ShapeMemoryModel3D<T>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m, \
                       eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, eqVolume, eqSurface, eqTileSpan,
                       persistenceLengthFine, eqLengthRatio, cellNumTriangles, cellNumVertices);
    } else  { // if (rbcModel == 1) {
       cellModel = new CellModel3D<T>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m, \
                       eqArea, eqLength, eqAngle, eqVolume, eqSurface, eqTileSpan,
                       persistenceLengthFine, eqLengthRatio, cellNumTriangles, cellNumVertices);
    }
    if (rbcModel == 0) {
        pcout << "mu_0 = " << dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->getMembraneShearModulus()*dNewton/dx << std::endl;
        pcout << "K = " << dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->getMembraneElasticAreaCompressionModulus()*dNewton/dx << std::endl;
        pcout << "YoungsModulus = " << dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->getYoungsModulus()*dNewton/dx << std::endl;
        pcout << "Poisson ratio = " << dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->getPoissonRatio() << std::endl;
    } else if (rbcModel == 1) {
        pcout << "mu_0 = " << dynamic_cast<CellModel3D<T>*>(cellModel)->getMembraneShearModulus()*dNewton/dx << std::endl;
        pcout << "K = " << dynamic_cast<CellModel3D<T>*>(cellModel)->getMembraneElasticAreaCompressionModulus()*dNewton/dx << std::endl;
        pcout << "YoungsModulus = " << dynamic_cast<CellModel3D<T>*>(cellModel)->getYoungsModulus()*dNewton/dx << std::endl;
        pcout << "Poisson ratio = " << dynamic_cast<CellModel3D<T>*>(cellModel)->getPoissonRatio() << std::endl;
    }


    if (goAndStop != 0) {
        goAndStopVertices.push_back(pluint(0));
        goAndStopVertices.push_back(pluint(cellNumVertices/2));
        goAndStopVertices.push_back(pluint(cellNumVertices-1));
    }

    pcout << std::endl << "Starting simulation" << std::endl;
    global::timer("sim").start();
    pcout << "Timer; iteration; LU; Cells; Vertices; Triangles; Processors; dt" << std::endl;

    /* Deflate if I say so */
    T eqVolumeInitial, eqVolumeFinal;
    eqVolumeInitial = eqVolume;
    eqVolumeFinal = deflationRatio * eqVolumeInitial ;
   if (rbcModel == 0) {
       applyProcessingFunctional (
           new ComputeShapeMemoryModelForce3D<T,DESCRIPTOR> (
               Cells, dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->clone(), cellsVolume, cellsSurface),
           immersedParticles.getBoundingBox(), particleArg );
   } else if (rbcModel == 1) {
       applyProcessingFunctional (
           new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
               Cells, dynamic_cast<CellModel3D<T>*>(cellModel)->clone(), cellsVolume, cellsSurface),
           immersedParticles.getBoundingBox(), particleArg );
   }


    writeCellLog(0, logFile,
                 cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                 cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity,
                 eqVolumeFinal, eqSurface, eqArea, eqLength) ; // Write Log file for the cell Particles

    // ======== Used for exploring cell stretching behaviour ========== //
    util::ValueTracer<T> convergeX((T)1,(T)100,1.0e-5), convergeY((T)1,(T)100,1.0e-5);
    T dStretchingForce = 5.0e-12/dNewton;
    plint timesToStretch = 40;
    if (fastStretchRelease) {
        dStretchingForce = 5.0e-12/dNewton;
        timesToStretch = 1;
    }
    plint statIter = 10;
    plint stretchReleased = 0;
    if ( (flowType == 4) or (flowType == 5) ) {
        stretchForce.resetToZero();
        stretchForceScalar = 0.0;
    }

    std::vector< std::vector<T> > cellInertia;
    std::vector< std::vector<T> > cellsEllipsoidFitSemiAxes, cellsEllipsoidFitAngles;
    std::vector<T> difference;
    computeEllipsoidFit (Cells, immersedParticles, cellIds, cellsCenter, cellsVolume,
                        cellsEllipsoidFitAngles, cellsEllipsoidFitSemiAxes, cellInertia, difference);

    T inertiaAnalytical =  2./5.*(eqVolume)*pow(radius,2);
    pcout << "INERTIA TENSOR (theoretical = " << inertiaAnalytical << ") " << inertiaAnalytical/cellInertia[0][0] << std::endl;
    for (plint i = 0; i < 3; ++i) {
            pcout <<  cellInertia[0][3*i] << "\t" <<  cellInertia[0][3*i +1] << "\t" <<  cellInertia[0][3*i + 2] << std::endl;
    }
    pcout << "a^2 + b^2  " << std::endl;
    pcout <<  cellsEllipsoidFitSemiAxes[0][0] << "\t" <<  cellsEllipsoidFitSemiAxes[0][1] << "\t" <<  cellsEllipsoidFitSemiAxes[0][2] << std::endl;
    pcout <<  cellsEllipsoidFitAngles[0][0]*180/pi << "\t" <<  cellsEllipsoidFitAngles[0][1]*180/pi << "\t" <<  cellsEllipsoidFitAngles[0][2]*180/pi << std::endl;

    SingleCellInShearFlow<T,DESCRIPTOR,DenseParticleField3D> shearFlow(Cells, immersedParticles, cellIds, cellsCenter, cellsVolume);
    shearFlow.set_dxdtdm(dx,dt,dm);
    if (flowType==6) {
        shearFlow.writeHeader(shearResultFile);
    }
    /* ********************* Main Loop ***************************************** * */
    for (pluint i=0; i<tmax+1; ++i) {
        if (goAndStop != 0 && i == pluint(tmax - 20/dt)) { // If stopAndGo experiment, turn off the shear
            shearRate = 0.05;
            changeCouetteShearRate(lattice, parameters, *boundaryCondition, shearRate);
            tmeas = 0.1/dt;
        }
        /* =============================== OUTPUT ===================================*/
        if (i%tmeas==0) {
            if (goAndStop != 0) {
                writeImmersedPointsVTK(Cells, goAndStopVertices, dx,
                    global::directories().getOutputDir()+createFileName("GoAndStop.",i,10)+".vtk");
            }
            // dtIteration = global::timer("sim").stop();
            // plint totParticlesNow = 0;
            // totParticlesNow = countParticles(immersedParticles, immersedParticles.getBoundingBox());
            // pcout << i << " totParticles = " << totParticles << std::endl;
            // PLB_ASSERT(totParticles == totParticlesNow); //Assert if some particles are outside of the domain
            calculateCellMeasures(Cells, immersedParticles, cellIds, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                                cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity, cellsMeanTileSpan);
            printCellMeasures(i, Cells, cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                                   cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity, eqVolumeFinal, eqSurface, eqArea, eqLength,
                                   dx, dt) ;
            writeCellLog(i, logFile,
                         cellsVolume, cellsSurface, cellsMeanTriangleArea, cellsMeanEdgeDistance,
                         cellsMaxEdgeDistance, cellsMeanAngle, cellsCenter, cellsVelocity,
                         eqVolumeFinal, eqSurface, eqArea, eqLength) ;
            if ((flowType == 3) || (flowType == 4) || (flowType == 5) ) {
                    applyProcessingFunctional (
                        new MeasureCellStretchDeformation3D<T,DESCRIPTOR>(lateralCellParticleTags, &stretchingDeformations),
                        immersedParticles.getBoundingBox(), particleArg );
                    pcout << "StrDeform [m] [";
                    for (pluint iDirection = 0; iDirection < stretchingDeformations.size(); ++iDirection) {
                        pcout << "(" << iDirection << ", " << stretchingDeformations[iDirection]*dx << "), " ;
                   } pcout << "]" << std::endl;
                   stretchLogFile << setprecision(20) << i*dt
                                  << "; " << stretchForceScalar*dNewton
                                  << "; " << stretchingDeformations[0]*dx
                                  << "; " << stretchingDeformations[1]*dx << ";  "
                                  << "; " << cellsMeanEdgeDistance[0] << ";  "
                                  << "; " << cellsMaxEdgeDistance[0] << ";  "
                                                                    << std::endl;
                   if (stretchReleased) {
                       stretchResealedFile << setprecision(20) << i*dt
                               << "; " << stretchForceScalar*dNewton
                               << "; " << stretchingDeformations[0]*dx
                               << "; " << stretchingDeformations[1]*dx << ";  "
                               << "; " << cellsMeanEdgeDistance[0] << ";  "
                               << "; " << cellsMaxEdgeDistance[0] << ";  "
                                                                       << std::endl;
                   }
            }
            std::vector<std::string> force_scalarNames;
            std::vector<std::string> velocity_scalarNames;
            std::vector<std::string> velocity_vectorNames;
            // writeMeshAsciiSTL(Cells, global::directories().getOutputDir()+createFileName("Mesh",i,8)+".stl");
            // serialize the particle information to write them.
            // a correspondance between the mesh and the particles is made. (Needs rescale)
            bool dynamicMesh = true;
            plint tag = -1; // Take all triangles.
            writeImmersedSurfaceVTK (
                Cells,
                immersedParticles,
                velocity_scalarNames, velocity_vectorNames,
                global::directories().getOutputDir()+createFileName("RBC",i,8)+".vtk", dynamicMesh, tag);
            writeVTK(lattice, parameters, i);
            // === Checkpoint ===
            //    parallelIO::save(immersedParticles, "immersedParticles.dat", true);
            //    parallelIO::save(Cells, "Cells.dat", true);
            //    parallelIO::load("immersedParticles.dat", immersedParticles, true);
            //    parallelIO::load("Cells.dat", Cells, true);
            // ==================
            global::timer("sim").restart();

            /* ====================== Update for single cell in shear flow ===================================*/
            if (flowType==6) {
                shearFlow.addData(i, Cells, immersedParticles, cellIds, cellsCenter, cellsVolume);
                shearFlow.write(shearResultFile);
            }

        }
        /* ====================== Update for Cell Stretching behaviour ===================================*/
        if (( (flowType == 4) or (flowType == 5) ) && (i % statIter == 0)) {
            applyProcessingFunctional (
                new MeasureCellStretchDeformation3D<T,DESCRIPTOR>(lateralCellParticleTags, &stretchingDeformations),
                immersedParticles.getBoundingBox(), particleArg );
            convergeX.takeValue(stretchingDeformations[0],false);
            convergeY.takeValue(stretchingDeformations[1],false);
            if (convergeX.hasConverged() && convergeY.hasConverged() ) {
                convergeX.resetValues();
                convergeY.resetValues();

                std::vector<std::string> force_scalarNames;
                std::vector<std::string> velocity_scalarNames;
                std::vector<std::string> velocity_vectorNames;
                if (not stretchReleased) {
                    writeImmersedSurfaceVTK (
                        Cells,
                        immersedParticles,
                        velocity_scalarNames, velocity_vectorNames,
                        global::directories().getOutputDir()+createFileName("stretch.", plint(stretchForceScalar*dNewton*1.0e12), 3)+".pN.vtk", true, -1);
                }
                else{
                    writeImmersedSurfaceVTK (
                        Cells,
                        immersedParticles,
                        velocity_scalarNames, velocity_vectorNames,
                        global::directories().getOutputDir()+createFileName("stretchR.", plint(stretchForceScalar*dNewton*1.0e12), 3)+".pN.vtk", true, -1);
                }
                stretchResultFile << setprecision(20) << i*dt
                        << "; " << stretchForceScalar*dNewton
                        << "; " << stretchingDeformations[0]*dx
                        << "; " << stretchingDeformations[1]*dx << ";  "
                        << "; " << cellsMeanEdgeDistance[0] << ";  "
                        << "; " << cellsMaxEdgeDistance[0] << ";  "
                                                                << std::endl;
                pcout << "# StretchForce: " << stretchForceScalar*dNewton*1.0e12 << "pN converged in iteration " << i << "." << std::endl;
                stretchForceScalar += dStretchingForce;
                stretchForce = Array<T,3>(stretchForceScalar,0,0);
                if (stretchForceScalar > timesToStretch*dStretchingForce ) {
                    if (not stretchReleased) {
                        stretchForce = Array<T,3>(0,0,0);
                        stretchReleased = 1;
                        stretchResealedFile << setprecision(20) << i*dt
                                << "; " << stretchForceScalar*dNewton
                                << "; " << stretchingDeformations[0]*dx
                                << "; " << stretchingDeformations[1]*dx << ";  "
                                << "; " << cellsMeanEdgeDistance[0] << ";  "
                                << "; " << cellsMaxEdgeDistance[0] << ";  "
                                << std::endl;
                        if (flowType == 5) {
                            if (rbcModel == 0) {
                                (dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel))->setMembraneShearViscosity(eta_mPrevious);
                            } else if (rbcModel == 1) {
                                (dynamic_cast<CellModel3D<T>*>(cellModel))->setMembraneShearViscosity(eta_mPrevious);
                            }
                        }
                    }
                    else {
                        stretchResealedFile << setprecision(20) << i*dt
                                << "; " << stretchForceScalar*dNewton
                                << "; " << stretchingDeformations[0]*dx
                                << "; " << stretchingDeformations[1]*dx << ";  "
                                << "; " << cellsMeanEdgeDistance[0] << ";  "
                                << "; " << cellsMaxEdgeDistance[0] << ";  "
                                << std::endl;
                        pcout << "Simulation is over.\n";
                        break;
                    }
                }
            }
        }
        /* =============================== MAIN ACTIONS ===================================*/
        // #0# Equilibration
        if (i<=relaxationTime) {
            eqVolume = eqVolumeInitial + (i*1.0)/relaxationTime * (eqVolumeFinal - eqVolumeInitial) ;
            if (rbcModel == 0) {
                (dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel))->setEquilibriumVolume(eqVolume);
            } else if (rbcModel == 1) {
                (dynamic_cast<CellModel3D<T>*>(cellModel))->setEquilibriumVolume(eqVolume);
            }
            if (flowType == 3) {
                    stretchForce = (i*1.0)/relaxationTime * Array<T,3>(stretchForceScalar,0,0);
            }
        }
        // #2# IBM Spreading
        if (forceToFluid != 0) { // Force from the Cell dynamics to the Fluid
            setExternalVector( lattice, lattice.getBoundingBox(),
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));
            applyProcessingFunctional ( // compute force applied on the fluid by the particles
                    new ForceToFluid3D<T,DESCRIPTOR> (ibmKernel),
                    immersedParticles.getBoundingBox(), particleLatticeArg );
        }
        // #3# LBM
        lattice.collideAndStream();
        // #4# IBM Interpolation
        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedCell3D<T,DESCRIPTOR>(ibmKernel),
            immersedParticles.getBoundingBox(), particleLatticeArg);
        // #5# Position Update
        applyProcessingFunctional ( // advance particles in time according to velocity
            new AdvanceParticlesFunctional3D<T,DESCRIPTOR>,
            immersedParticles.getBoundingBox(), particleArg );
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        // #1# Membrane Model
        if (rbcModel == 0) {
            applyProcessingFunctional (
                new ComputeShapeMemoryModelForce3D<T,DESCRIPTOR> (
                    Cells, dynamic_cast<ShapeMemoryModel3D<T>*>(cellModel)->clone(), cellsVolume, cellsSurface),
                immersedParticles.getBoundingBox(), particleArg );
        } else if (rbcModel == 1) {
            applyProcessingFunctional (
                new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
                    Cells, dynamic_cast<CellModel3D<T>*>(cellModel)->clone(), cellsVolume, cellsSurface),
                immersedParticles.getBoundingBox(), particleArg );
        }

        if ((flowType == 3) || (flowType == 4) || (flowType == 5)) {
            applyProcessingFunctional ( // compute force applied on the some particles by the stretching force
                    new ApplyStretchingForce3D<T,DESCRIPTOR>(outerLeftTags, outerRightTags, stretchForce, cellModel->getDensity()),
                    immersedParticles.getBoundingBox(), particleArg );
        }
    }
    pcout << "Simulation finished." << std::endl;
//    MPI_Finalize();
}
