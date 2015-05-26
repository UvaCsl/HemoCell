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

#include "ficsion.h"
#include "initRBCsAndPLTs.hh"
typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor


template<typename T>
class StenosisDomain3D : public DomainFunctional3D {
public:
    StenosisDomain3D(plint nx_, plint ny_, plint nz_, plint stenosisHeightLU_, plint stenosisLengthLU_)
        : nx(nx_),
          ny(ny_),
          nz(nz_),
          stenosisHeightLU(stenosisHeightLU_),
          stenosisLengthLU(stenosisLengthLU_)
    { }
    // The function-call operator is overridden to specify the location
    //   of bounce-back nodes.
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        if (0==iY) { return 1; };
        plint smoothingFactor=5;
        if (iX < stenosisLengthLU) {
            bool toReturn=0;
            T midpoint = ny/2.0;
            T halfHeight = stenosisHeightLU/2;
            T midpointDistance = sqrt( (iY-midpoint)*(iY-midpoint) );
            if (midpointDistance < halfHeight ) { return 0; };
            if (midpointDistance > halfHeight+smoothingFactor ) { return 1; };
            if ( (iX >= smoothingFactor) and (iX <= stenosisLengthLU-smoothingFactor) ) { return 1; }
            if ( (iX < smoothingFactor) and (midpointDistance >= halfHeight+smoothingFactor - iX ) ) { return 1; };
            if ( (iX > stenosisLengthLU-smoothingFactor) and (midpointDistance >= halfHeight+smoothingFactor + (iX-stenosisLengthLU) ) ) { return 1; };
            return 0;

            // Smooth
            if (midpoint)
                if (iY > (ny-stenosisHeightLU-2*smoothingFactor)/2 ) { toReturn=0; };
                if (iY < (ny+stenosisHeightLU+2*smoothingFactor)/2 ) { toReturn=0; };
                if (iX>smoothingFactor and iX < stenosisLengthLU-smoothingFactor) { return 1; }

                vector<plint> centersX(4), centersY(4);
                centersX[0] = smoothingFactor; centersY[0] = (ny-stenosisHeightLU-2*smoothingFactor)/2;
                centersX[1] = smoothingFactor; centersY[1] = (ny+stenosisHeightLU+2*smoothingFactor)/2;
                centersX[2] = stenosisLengthLU-smoothingFactor; centersY[2] = (ny-stenosisHeightLU-2*smoothingFactor)/2;
                centersX[3] = stenosisLengthLU-smoothingFactor; centersY[3] = (ny+stenosisHeightLU+2*smoothingFactor)/2;
                for (int ic = 0; ic < 4; ++ic) {
                    T cx = (centersX[ic]-iX); T cy = (centersY[ic] - iY);
                    if ( cx*cx + cy*cy < smoothingFactor*smoothingFactor ) { return 1; }
                }
        }
        return 0;
    }
    virtual StenosisDomain3D<T>* clone() const {
        return new StenosisDomain3D<T>(*this);
    }
private:
    plint nx, ny, nz;
    plint stenosisHeightLU, stenosisLengthLU;
};


template<typename T, template<class U> class Descriptor>
void iniLatticePoiseuilleStenosisWithBodyForce(MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T poiseuilleForce, T stenosisHeightLU, T stenosisLengthLU)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D left   = Box3D(0, nx-1, 0,    0,    0, nz-1);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);

    pcout << "(iniLatticePoiseuilleStenosisWithBodyForce) initializing lattice" << std::endl;
    defineDynamics( lattice, lattice.getBoundingBox(), new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    defineDynamics( lattice, lattice.getBoundingBox(),
                    new StenosisDomain3D<T>(nx, ny, nz, stenosisHeightLU, stenosisLengthLU),
                    new BounceBack<T,Descriptor> );

    T rhoInit=1.0; Array<T,3> uInit(0.0,0.0,0.0);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rhoInit, uInit);

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(poiseuilleForce, 0.0, 0.0));

    lattice.initialize();

}



void readFicsionXML(XMLreader & documentXML,std::string & caseId, plint & rbcModel, T & shellDensity, T & k_rest,
        T & k_shear, T & k_bend, T & k_stretch, T & k_WLC, T & eqLengthRatio, T & k_rep, T & k_elastic, T & k_volume, T & k_surface, T & eta_m,
        T & rho_p, T & u, plint & flowType, T & Re, T & shearRate, T & stretchForce, Array<T,3> & eulerAngles, T & Re_p, T & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & ibmKernel, plint & ibmScheme, plint & shape, std::string & cellPath, T & radius, T & deflationRatio, pluint & relaxationTime,
        plint & minNumOfTriangles, pluint & tmax, plint & tmeas, T & hct, plint & npar, plint & flowParam, bool & checkpointed, T & stenosisHeightLU, T & stenosisLengthLU)
    {
    T nu_p, tau, dx;
    T dt, nu_lb;
    std::string firstField = (*(documentXML.getChildren( documentXML.getFirstId() )[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    if (firstField=="ficsion") { checkpointed = 0; }
    else { checkpointed = 1; }

    XMLreaderProxy document = checkpointed?documentXML["Checkpoint"]["ficsion"]:documentXML["ficsion"];
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
    flowType = 0;
    document["parameters"]["Re"].read(Re);
    document["parameters"]["shearRate"].read(shearRate);
    document["parameters"]["stretchForce"].read(stretchForce); // In picoNewton
    stretchForce *= 1e-12;
    std::vector<T> ea;
    document["parameters"]["eulerAngles"].read(ea);
    if (ea.size() != 3) {
        ea.resize(3, 0.0);
    }
    eulerAngles = Array<T,3>(ea[0], ea[1], ea[2]);
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
    try {
        document["ibm"]["ibmScheme"].read(ibmScheme);
    } catch(const plb::PlbIOException & message) {
        ibmScheme=0;
    }
    document["domain"]["rhoP"].read(rho_p);
    document["domain"]["nuP"].read(nu_p);
    document["domain"]["tau"].read(tau);
    document["domain"]["dx"].read(dx);
    // Read lx, ly, lz --or nx, ny, nz
    try {
        document["domain"]["lx"].read(lx);
        document["domain"]["ly"].read(ly);
        document["domain"]["lz"].read(lz);
    } catch(const plb::PlbIOException & message) {
        T nx, ny, nz;
        document["domain"]["nx"].read(nx);
        document["domain"]["ny"].read(ny);
        document["domain"]["nz"].read(nz);
        lx = nx * dx;
        ly = ny * dx;
        lz = nz * dx;
    }

    T stenosisHeightFraction, stenosisLengthFraction;
    document["domain"]["stenosisHeightFraction"].read(stenosisHeightFraction);
    document["domain"]["stenosisLengthFraction"].read(stenosisLengthFraction);


    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    try {
        document["sim"]["npar"].read(npar);
        hct = 0;
    } catch(const plb::PlbIOException & message) {
        document["sim"]["hct"].read(hct);
        npar = 0;
    }
    hct /= 100.0;

    radius = radius*1.0/dx; // Transform from [m] to [LU]
    nu_lb = (tau-0.5)/3.0;
    dt = (nu_lb/nu_p)*dx*dx;
    u = dt*1.0/dx;
    Re_p = 1.0/nu_p;
    N = int(1.0/dx);
    flowParam = 0; // flowType/10;
//    flowType = flowType%10;
    if ( (flowType == 3) or (flowType == 4) or (flowType == 5) or (flowType == 8) ) { // Cell Stretching Analysis
        shearRate = 0;
        Re = 0;
    } else if (flowType == 6) { // Tumbling Tank Treading Measurements
        if (shearRate > 0) {
            tmax  = 100.0/(shearRate*dt); // Measurements for 100 Strain rates
            tmeas = 0.02/(shearRate*dt); // 50 Measurements per Strain rates
        }
    }
    if (flowParam == 7) {
        tmax += 20.0/dt; // Add 20 seconds to the total time for Fischer2004.
    }

    stenosisHeightLU = stenosisHeightFraction * ly / dx;
    stenosisLengthLU = stenosisLengthFraction * lx / dx;

}





int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::timer("ficsion_init").start();

    global::directories().setOutputDir("./tmp/");
    global::directories().setLogOutDir("./log/");
    global::directories().setInputDir("./");

    global::IOpolicy().activateParallelIO(true);
    global::IOpolicy().setStlFilesHaveLowerBound(false);
//    global::IOpolicy().setLowerBoundForStlFiles(-1.);
//    testInPlane(); PLB_ASSERT(false);
    std::string outputDir = global::directories().getOutputDir();
    std::string inputDir = global::directories().getInputDir();
    std::string logOutDir = global::directories().getLogOutDir();
    mkpath((outputDir + "/hdf5/").c_str(), 0777);
    mkpath(logOutDir.c_str(), 0777);

    plint forceToFluid, shape, cellNumTriangles, ibmKernel, ibmScheme;
    plint rbcModel;
    std::string caseId;
    std::string cellPath;
    pluint tmax;
    plint tmeas;
    T hct = 0;
    plint npar = 0;
//    T dtIteration = 0;
    T shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_rep, k_elastic,  k_volume, k_surface, eta_m;
    T eqLengthRatio;
    T u, Re, Re_p, N, lx, ly, lz;
    T poiseuilleForce=0;
    T rho_p;
    T radius;
    T deflationRatio;
    pluint relaxationTime;
    plint flowType;
    T shearRate, shearRate_p;
    Array<T,3> stretchForce(0,0,0);
    T stretchForceScalar;
    T stretchForce_p;
    Array<T,3> eulerAngles;
    plint flowParam;
    bool checkpointed=0;
    T stenosisHeightLU, stenosisLengthLU;
    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "(main) reading.." <<std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
            k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface, eta_m,
            rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz,  forceToFluid, ibmKernel, ibmScheme, shape, cellPath, radius, deflationRatio, relaxationTime,
            cellNumTriangles, tmax, tmeas, hct, npar, flowParam, checkpointed, stenosisHeightLU, stenosisLengthLU);
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
    shearRate = shearRate_p * dt;
    stretchForceScalar = stretchForce_p / dNewton;
    pcout << "(main) dx = " << dx << ", " <<
             "dt = " << dt << ", " <<
             "dm = " << dt << ", " <<
             "kT = " << kBT <<
             std::endl;

    /* ------------------ *
     * Initialize Lattice *
     * ------------------ */
    plint extendedEnvelopeWidth=4;
//    if (ibmKernel==2) {
//        extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with width 2.
//    }
//    else {
//        extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with width 2.
//    }
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    lattice.periodicity().toggleAll(true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    T nu_tmp = parameters.getLatticeNu();
    poiseuilleForce = 8 * (nu_tmp*nu_tmp) * Re / (stenosisHeightLU*stenosisHeightLU*stenosisHeightLU) ;
    T Umax = poiseuilleForce * (stenosisHeightLU-1) * (stenosisHeightLU-1) / (8*parameters.getLatticeNu());
    tmeas = plint(stenosisHeightLU*1.0/Umax)/40;
    iniLatticePoiseuilleStenosisWithBodyForce<T, DESCRIPTOR>(lattice, parameters, *boundaryCondition, poiseuilleForce, stenosisHeightLU-2, stenosisLengthLU+2);

    pcout << getMultiBlockInfo(lattice) << std::endl;

    /*
     * Initialize model *
     */

    /* CREATE MESH */
//    Array<T,3> eqShapeRotationEulerAngles = eulerAngles;
    if (flowParam == 9) { eulerAngles= Array<T,3>(0.,0.,0.); }
    // Radius in LU


    T persistenceLengthFine = 7.5e-9 ; // In meters
    k_rest= 0;

    std::vector<ConstitutiveModel<T, DESCRIPTOR>* > cellModels;
    std::vector<CellField3D<T, DESCRIPTOR>* > cellFields;

    std::vector<T> eqVolumes;
//    =======================  Create RBC
    TriangleBoundary3D<T> Cells = constructMeshElement(shape, radius, cellNumTriangles, dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = Cells.getMesh();
    MeshMetrics<T> meshmetric(meshElement);    meshmetric.write();
    eqVolumes.push_back(meshmetric.getVolume());
    plint numVerticesPerCell = meshElement.getNumVertices();
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    cellModels.push_back(new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
        persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement));
    cellFields.push_back(new CellField3D<T, DESCRIPTOR>(lattice, meshElement, hct, cellModels[0], ibmKernel, "RBC"));

//    =======================  Create Platelet
        T pltRadius = 1.15e-6/dx;
        T aspectRatio = 1.0 / (2*pltRadius);
        TriangleBoundary3D<T> PLTCells = constructMeshElement(6, pltRadius, 200, dx, cellPath, eulerAngles, aspectRatio);
        TriangularSurfaceMesh<T> pltMeshElement = PLTCells.getMesh();
        eqVolumes.push_back(MeshMetrics<T>(pltMeshElement).getVolume());
        cellModels.push_back(new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend*5, k_stretch, k_WLC*5.0, k_elastic, k_volume, k_surface, eta_m,
            persistenceLengthFine, eqLengthRatio, dx, dt, dm, pltMeshElement) );
        cellFields.push_back(new CellField3D<T, DESCRIPTOR>(lattice, pltMeshElement, 0.001, cellModels[cellModels.size()-1], ibmKernel, "PLT"));



    CellField3D<T, DESCRIPTOR> & RBCField = *cellFields[0];
    CellStretch<T, DESCRIPTOR> cellStretch(RBCField, stretchForceScalar, 0.1);

    FcnCheckpoint<T, DESCRIPTOR> checkpointer(document);
    plint initIter=0;
    checkpointer.load(document, lattice, cellFields, initIter);
    if (not checkpointer.wasCheckpointed()) {
        pcout << "(main) initializing"<< std::endl;
        std::vector<Array<T,3> > cellsOrigin;
        cellsOrigin.push_back( Array<T,3>(nx*0.5, ny*0.5, nz*0.5) );
        if (flowType==3) { RBCField.initialize(cellsOrigin); }
        else if (hct>0) {
//            RBCField.grow(0);
//           orderedPositionCellField3D(cellFields);
//            randomPositionCellFieldsForGrowth3D(cellFields);
            orderedPositionMultipleCellField3D(cellFields);
//            Dot3D latticeSize(nx, ny, nz);
//            positionOrderedRBCsAndPLTs3D(cellFields, latticeSize);
        }
        else { RBCField.initialize(cellsOrigin); }
        // IMPORTANT: Runs fluid for 1k iterations to kickstart the flow.
        iniLatticePoiseuilleStenosisWithBodyForce<T, DESCRIPTOR>(lattice, parameters, *boundaryCondition, poiseuilleForce, stenosisHeightLU, stenosisLengthLU);
         for (plint itrt=0; itrt<200; ++itrt) { lattice.collideAndStream(); }
        checkpointer.save(lattice, cellFields, initIter);
    }
    for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
        cellFields[iCell]->setParticleUpdateScheme(ibmScheme);
    }

//    if (rbcModel == 3) {
//        // Has a problem with checkpointing
//      (dynamic_cast<RestModel3D<T,DESCRIPTOR>*>(cellModel))->freezeVertices(RBCField);
//    }
    pcout << std::endl ;
    T domainVol = nx*ny*nz; // computeSum(*flagMatrix);

    std::auto_ptr<plb::MultiScalarField3D<T> > flagMatrix = computeDensity(lattice, lattice.getBoundingBox());
    domainVol =  computeSum(*flagMatrix);

    pcout << "(main) Number of fluid cells: " << domainVol << ", " << domainVol*1.0/(nx*ny*nz) << "%." << std::endl;
    pcout << "(main) Ma= " << Umax*sqrt(3.0) << std::endl;
    pcout << "(main) Re= " << Re << std::endl;
    for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
        plint nCells = cellFields[iCell]->getNumberOfCells_Global();
        pcout << "(main) Hematocrit [x100%]: " << nCells*eqVolumes[iCell]*100.0/(domainVol) << std::endl;
        pcout << "(main) nCells (global) = " << nCells << ", pid: " << global::mpi().getRank() ;
        pcout << ", Volume = " << eqVolumes[iCell] << std::endl;
    }
//    tmeas = 1;
    pcout << std::endl << "(main) Starting simulation i=" << initIter << ", tmeas=" << tmeas << std::endl;

    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
                                                        createBoundaryParticleField3D(lattice);

    /* Repulsive force */
    T k_int = 0.005, DeltaX=1.0, R=1.7, k=1.5;
    PowerLawForce<T> PLF(k_int, DeltaX, R, k);

    /*      Sync all quantities    */
    SyncRequirements everyCCR(allReductions);
    for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
        cellFields[iCell]->synchronizeCellQuantities(everyCCR);
    }
    /*            I/O              */
    global::timer("HDFOutput").start();
    writeHDF5(lattice, parameters, initIter);
    for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
        writeCellField3D_HDF5(*cellFields[iCell], dx, dt, initIter);
        writeCell3D_HDF5(*cellFields[iCell], dx, dt, initIter);
    }

    global::timer("HDFOutput").stop();

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, -1, numVerticesPerCell);
    /* --------------------------- */
    global::timer("mainLoop").start();
    global::profiler().turnOn();
//#ifdef PLB_DEBUG // Palabos has this bug. It's missing the "envelope-update" is the profiler.
    if (flowType==2) { global::profiler().turnOff(); }
//#endif
    for (pluint iter=initIter; iter<tmax+1; ++iter) {
        // #1# Membrane Model
//       RBCField.applyConstitutiveModel();
//       RBCField.applyCellCellForce(PLF, R);
        for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
           cellFields[iCell]->applyConstitutiveModel();
           applyWallCellForce<T, DESCRIPTOR>(PLF, R, *boundaryParticleField3D, *cellFields[iCell]);
        }

       if (flowType==3) { cellStretch.stretch(); }
        // #2# IBM Spreading
       cellFields[0]->setFluidExternalForce(poiseuilleForce);
       for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
           cellFields[iCell]->spreadForceIBM();
       }
        // #3# LBM
        if ((iter+1)%tmeas==0 && flowType==11) { lattice.toggleInternalStatistics(true); }
        global::timer("LBM").start();
        lattice.collideAndStream();
        global::timer("LBM").stop();
        for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
            // #4# IBM Interpolation
            cellFields[iCell]->interpolateVelocityIBM();
            // #5# Position Update
            cellFields[iCell]->advanceParticles();
        }

        // #6# Output
        if ((iter+1)%tmeas==0) {
            SyncRequirements everyCCR(allReductions);
            for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
                cellFields[iCell]->synchronizeCellQuantities(everyCCR);
            }
            global::timer("HDFOutput").start();
            writeHDF5(lattice, parameters, iter+1);
            for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
                writeCellField3D_HDF5(*cellFields[iCell], dx, dt, iter+1);
                writeCell3D_HDF5(*cellFields[iCell], dx, dt, iter+1);
            }
            global::timer("HDFOutput").stop();
            if ((iter+1)%(2*tmeas)==0) {
                global::timer("Checkpoint").start();
                checkpointer.save(lattice, cellFields, iter+1);
                global::timer("Checkpoint").stop();
            }
            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter+1);
            global::profiler().writeReport();
            pcout << "(main) Iteration:" << iter + 1 << "; time "<< dtIteration*1.0/tmeas ;
//            pcout << "; Volume (" << RBCField[0]->getVolume() << ")";
            pcout << std::endl;
        } else {
            for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
                cellFields[iCell]->synchronizeCellQuantities();
            }
        }
        if ((iter+1)%tmeas==0 && flowType==11) { lattice.toggleInternalStatistics(false); }
    }
    for (pluint iCell=0; iCell<cellFields.size(); ++iCell) {
        delete cellFields[iCell];
        delete cellModels[iCell];
    }
    simpleProfiler.writeIteration(tmax+1);
    global::profiler().writeReport();
    pcout << "Simulation finished." << std::endl;
}
