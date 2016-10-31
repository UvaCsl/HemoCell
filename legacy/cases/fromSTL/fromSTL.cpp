#include "ficsion.h"

typedef double T;
typedef Array<T, 3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

void readFicsionXML(XMLreader &documentXML, std::string &caseId, plint &rbcModel, T &shellDensity, T &k_rest,
                    T &k_shear, T &k_bend, T &k_stretch, T &k_WLC, T &eqLengthRatio, T &k_rep, T &k_elastic,
                    T &k_volume, T &k_surface, T &eta_m,
                    T &rho_p, T &u, plint &flowType, T &Re, T &shearRate, T &stretchForce, Array<T, 3> &eulerAngles,
                    T &Re_p, T &N, T &lx, T &ly, T &lz,
                    plint &forceToFluid, plint &ibmKernel, plint &ibmScheme, plint &shape, std::string &cellPath,
                    T &radius, T &deflationRatio, pluint &relaxationTime,
                    plint &minNumOfTriangles, pluint &tmax, plint &tmeas, T &hct, plint &npar, plint &flowParam,
                    bool &checkpointed) {
    T nu_p, tau, dx;
    T dt, nu_lb;
    std::string firstField = (*(documentXML.getChildren(
            documentXML.getFirstId())[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    if (firstField == "ficsion") { checkpointed = 0; }
    else { checkpointed = 1; }

    XMLreaderProxy document = checkpointed ? documentXML["Checkpoint"]["ficsion"] : documentXML["ficsion"];
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
    eulerAngles = Array<T, 3>(ea[0], ea[1], ea[2]);
    eulerAngles[0] *= pi / 180.;
    eulerAngles[1] *= pi / 180.;
    eulerAngles[2] *= pi / 180.;
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
    } catch (const plb::PlbIOException &message) {
        ibmScheme = 0;
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
    } catch (const plb::PlbIOException &message) {
        T nx, ny, nz;
        document["domain"]["nx"].read(nx);
        document["domain"]["ny"].read(ny);
        document["domain"]["nz"].read(nz);
        lx = nx * dx;
        ly = ny * dx;
        lz = nz * dx;
    }
    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    try {
        document["sim"]["npar"].read(npar);
        hct = 0;
    } catch (const plb::PlbIOException &message) {
        document["sim"]["hct"].read(hct);
        npar = 0;
    }
    hct /= 100.0;

    radius = radius * 1.0 / dx; // Transform from [m] to [LU]
    nu_lb = (tau - 0.5) / 3.0;
    dt = (nu_lb / nu_p) * dx * dx;
    u = dt * 1.0 / dx;
    Re_p = 1.0 / nu_p;
    N = int(1.0 / dx);
    flowParam = 0; // flowType/10;
//    flowType = flowType%10;
    if ((flowType == 3) or (flowType == 4) or (flowType == 5) or (flowType == 8)) { // Cell Stretching Analysis
        shearRate = 0;
        Re = 0;
    } else if (flowType == 6) { // Tumbling Tank Treading Measurements
        if (shearRate > 0) {
            tmax = 100.0 / (shearRate * dt); // Measurements for 100 Strain rates
            tmeas = 0.02 / (shearRate * dt); // 50 Measurements per Strain rates
        }
    }
    if (flowParam == 7) {
        tmax += 20.0 / dt; // Add 20 seconds to the total time for Fischer2004.
    }
}

//VoxelizedDomain3D(TriangleBoundary3D<T> const& boundary_,
//                  int flowType_, plint extraLayer_, plint borderWidth_,
//                  plint envelopeWidth_, plint blockSize_,
//                  plint gridLevel_=0, bool dynamicMesh_ = false);




template<typename Tp>
class CopyFromNeighbor : public BoxProcessingFunctional3D_S<Tp> {
public:
    CopyFromNeighbor(Array<plint, 3> offset_) : offset(offset_) { };

    virtual void process(Box3D domain, ScalarField3D<Tp> &field1);

    virtual CopyFromNeighbor<Tp> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

    virtual BlockDomain::DomainT appliesTo() const;

private:
    Array<plint, 3> offset;
};


/* *************** Data Functionals for scalar-fields **************** */

template<typename Tp>
void CopyFromNeighbor<Tp>::process(
        Box3D domain, ScalarField3D<Tp> &field1) {
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field1.get(iX, iY, iZ) = field1.get(iX + offset[0], iY + offset[1], iZ + offset[2]);
            }
        }
    }
}

template<typename Tp>
CopyFromNeighbor<Tp> *CopyFromNeighbor<Tp>::clone() const {
    return new CopyFromNeighbor<Tp>(*this);
}

template<typename Tp>
void CopyFromNeighbor<Tp>::getTypeOfModification(std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::allVariables;
}

template<typename Tp>
BlockDomain::DomainT CopyFromNeighbor<Tp>::appliesTo() const {
    return BlockDomain::bulk;
}


void getFlagMatrixFromSTL(std::string meshFileName, plint extendedEnvelopeWidth, plint nyRef,
                          VoxelizedDomain3D<T> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix) {
    plint yDirection = 1; // Reference direction for stl resolution
    plint extraLayer = 0;  // Make the bounding box larger; for visualization purposes
    //   only. For the simulation, it is OK to have extraLayer=0.
    plint blockSize = -1; // Zero means: no sparse representation.
    plint borderWidth = 1;  // Because the Guo boundary condition acts in a one-cell layer.
    // Requirement: margin>=borderWidth.
    plint margin = 1;  // Extra margin of allocated cells around the obstacle.

    TriangleSet<T> *triangleSet = new TriangleSet<T>(meshFileName, DBL);

    DEFscaledMesh<T> *defMesh =
            new DEFscaledMesh<T>(*triangleSet, nyRef, yDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    voxelizedDomain = new VoxelizedDomain3D<T>(
            boundary, voxelFlag::inside, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain->getVoxelMatrix()) << std::endl;


    flagMatrix = new MultiScalarField3D<int>((MultiBlock3D &) voxelizedDomain->getVoxelMatrix());
    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix->getBoundingBox(), 1);
    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix->getBoundingBox(), 1);


    Box3D domainBox = flagMatrix->getBoundingBox();
    plint nx = domainBox.getNx();
    plint ny = domainBox.getNy();
    plint nz = domainBox.getNz();

    Box3D domain(0, 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor<int>(Array<plint, 3>(1, 0, 0)), domain, *flagMatrix);

    domain = Box3D(nx - 2, nx - 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor<int>(Array<plint, 3>(-1, 0, 0)), domain, *flagMatrix);


    VtkImageOutput3D<T> vtkOut("test.vtk");
    vtkOut.writeData<float>(*copyConvert<int, double>(*extractSubDomain(*flagMatrix, flagMatrix->getBoundingBox())),
                            "flag", 1.);

}


int main(int argc, char *argv[]) {
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
    T shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_rep, k_elastic, k_volume, k_surface, eta_m;
    T eqLengthRatio;
    T u, Re, Re_p, N, lx, ly, lz;
    T poiseuilleForce = 0;
    T rho_p;
    T radius;
    T deflationRatio;
    pluint relaxationTime;
    plint flowType;
    T shearRate, shearRate_p;
    Array<T, 3> stretchForce(0, 0, 0);
    T stretchForceScalar;
    T stretchForce_p;
    Array<T, 3> eulerAngles;
    plint flowParam;
    bool checkpointed = 0;

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "(main) reading.." << std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
                   k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface,
                   eta_m,
                   rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz, forceToFluid,
                   ibmKernel, ibmScheme, shape, cellPath, radius, deflationRatio, relaxationTime,
                   cellNumTriangles, tmax, tmeas, hct, npar, flowParam, checkpointed);
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
    T dm = rho_p * (dx * dx * dx);
    dNewton = (dm * dx / (dt * dt));
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
    plint extendedEnvelopeWidth = 4;
//    if (ibmKernel==2) {
//        extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with width 2.
//    }
//    else {
//        extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with width 2.
//    }
    std::string meshFileName = "tube.stl";
    plb::MultiScalarField3D<int> *flagMatrix = 0;
    VoxelizedDomain3D<T> *voxelizedDomain = 0;
    getFlagMatrixFromSTL(meshFileName, extendedEnvelopeWidth, ny, voxelizedDomain, flagMatrix);
    Box3D domainBox = flagMatrix->getBoundingBox();
    nx = domainBox.getNx();
    ny = domainBox.getNy();
    nz = domainBox.getNz();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
//            defaultMultiBlockPolicy3D().getMultiBlockManagement(flagMatrix->getBoundingBox(), extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega()));

//    MultiBlockLattice3D<T,DESCRIPTOR> lattice = *generateMultiBlockLattice<T,DESCRIPTOR> (
//                voxelizedDomain->getVoxelMatrix(), extendedEnvelopeWidth, new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));

    defineDynamics(lattice, *flagMatrix, lattice.getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

    lattice.periodicity().toggleAll(true);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 1., Array<T, 3>(0., 0., 0.));

    T nu_tmp = parameters.getLatticeNu();
    poiseuilleForce = 8 * (nu_tmp * nu_tmp) * Re / (ny * ny * ny);
    if (tmeas == 0) {
        tmeas = plint(ny * ny * 1.0 / (4 * nu_tmp * Re)) / 40;
    }

    setExternalVector(lattice, lattice.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    lattice.initialize();

    /*
     * Initialize model *
     */

    /* CREATE MESH */
//    Array<T,3> eqShapeRotationEulerAngles = eulerAngles;
    if (flowParam == 9) { eulerAngles = Array<T, 3>(0., 0., 0.); }
    // Radius in LU


    T persistenceLengthFine = 7.5e-9; // In meters
    k_rest = 0;

    std::vector<ConstitutiveModel<T, DESCRIPTOR> *> cellModels;
    std::vector<CellField3D<T, DESCRIPTOR> *> cellFields;

    std::vector<T> eqVolumes;
//    =======================  Create RBC
    TriangleBoundary3D<T> RBCCells = constructMeshElement(shape, radius, cellNumTriangles, dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = RBCCells.getMesh();
    MeshMetrics<T> meshmetric(meshElement);
    meshmetric.write();
    eqVolumes.push_back(meshmetric.getVolume());
    plint numVerticesPerCell = meshElement.getNumVertices();
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    cellModels.push_back(
            new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic,
                                                  k_volume, k_surface, eta_m,
                                                  persistenceLengthFine, eqLengthRatio, dx, dt, dm, meshElement));
    cellFields.push_back(new CellField3D<T, DESCRIPTOR>(lattice, meshElement, hct, cellModels[0], ibmKernel, "RBC"));

//    =======================  Create Platelet
    T pltRadius = 1.15e-6 / dx;
    T aspectRatio = 1.0 / (2 * pltRadius);
    TriangleBoundary3D<T> PLTCells = constructMeshElement(6, pltRadius, 200, dx, cellPath, eulerAngles, aspectRatio);
    TriangularSurfaceMesh<T> pltMeshElement = PLTCells.getMesh();
    eqVolumes.push_back(MeshMetrics<T>(pltMeshElement).getVolume());
    cellModels.push_back(
            new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend * 5, k_stretch, k_WLC * 5.0,
                                                  k_elastic, k_volume, k_surface, eta_m,
                                                  persistenceLengthFine, eqLengthRatio, dx, dt, dm, pltMeshElement));
    cellFields.push_back(new CellField3D<T, DESCRIPTOR>(lattice, pltMeshElement, 0.0025 * hct,
                                                        cellModels[cellModels.size() - 1], ibmKernel, "PLT"));


    CellField3D<T, DESCRIPTOR> &RBCField = *cellFields[0];
    CellStretch<T, DESCRIPTOR> cellStretch(RBCField, stretchForceScalar, 0.1);

    FcnCheckpoint<T, DESCRIPTOR> checkpointer(document);
    plint initIter = 0;
    checkpointer.load(document, lattice, cellFields, initIter);
    if (not checkpointer.wasCheckpointed()) {
        pcout << "(main) initializing" << std::endl;
        std::vector<Array<T, 3> > cellsOrigin;
        cellsOrigin.push_back(Array<T, 3>(nx * 0.5, ny * 0.5, nz * 0.5));
        if (flowType == 3) { RBCField.initialize(cellsOrigin); }
        else if (hct > 0) {
//            RBCField.grow(0);
//           orderedPositionCellField3D(cellFields);
//            randomPositionCellFieldsForGrowth3D(cellFields);
            //orderedPositionMultipleCellField3D(cellFields);
            randomPositionMultipleCellField3D(cellFields, hct);
        }
        else { RBCField.initialize(cellsOrigin); }
        checkpointer.save(lattice, cellFields, initIter);
    }
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        cellFields[iCell]->setParticleUpdateScheme(ibmScheme);
    }

//    if (rbcModel == 3) {
//        // Has a problem with checkpointing
//      (dynamic_cast<RestModel3D<T,DESCRIPTOR>*>(cellModel))->freezeVertices(RBCField);
//    }
    pcout << std::endl;
    plint domainVol = computeSum(*flagMatrix);
    pcout << "Number of fluid cells: " << domainVol << std::endl;
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        plint nCells = cellFields[iCell]->getNumberOfCells_Global();
        pcout << "(main) Hematocrit [x100%]: " << nCells * eqVolumes[iCell] * 100.0 / (domainVol) << std::endl;
        pcout << "(main) nCells (global) = " << nCells << ", pid: " << global::mpi().getRank();
        pcout << ", Volume = " << eqVolumes[iCell] << std::endl;
    }

    /* Repulsive force */
    T k_int = 0.00025, DeltaX = 1.0, R = 0.75, k = 1.5;
    PowerLawForce<T> PLF(k_int, DeltaX, R, k);

    /*      Sync all quantities    */
    SyncRequirements everyCCR(allReductions);
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        cellFields[iCell]->synchronizeCellQuantities(everyCCR);
    }
    /*            I/O              */
    global::timer("HDFOutput").start();
    bool invertXZ_for_XDMF = true;
    writeHDF5(lattice, parameters, initIter, invertXZ_for_XDMF);
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        writeCellField3D_HDF5(*cellFields[iCell], dx, dt, initIter);
        writeCell3D_HDF5(*cellFields[iCell], dx, dt, initIter);
    }
    global::timer("HDFOutput").stop();


    pcout << std::endl << "(main) Starting simulation i=" << initIter << ", tmeas=" << tmeas << std::endl;
    MultiParticleField3D<LightParticleField3D<T, DESCRIPTOR> > *boundaryParticleField3D =
            createBoundaryParticleField3D(lattice);
    writeParticleField3D_HDF5(*boundaryParticleField3D, dx, dt, 0, "BoundaryParticles");

    for (plint itrt = 0; itrt < 200; ++itrt) { lattice.collideAndStream(); }

//    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
//                                                        createBoundaryParticleField3D(lattice);

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, -1, numVerticesPerCell);
    /* --------------------------- */
    global::timer("mainLoop").start();
    for (pluint iter = initIter; iter < tmax + 1; ++iter) {
        // #1# Membrane Model
//       RBCField.applyConstitutiveModel();
//       RBCField.applyCellCellForce(PLF, R);
        for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
            cellFields[iCell]->applyConstitutiveModel();
        }

        if (flowType == 3) { cellStretch.stretch(); }
        // #2# IBM Spreading
        cellFields[0]->setFluidExternalForce(poiseuilleForce);
        for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
            cellFields[iCell]->spreadForceIBM();
        }
        // #3# LBM
        if ((iter + 1) % tmeas == 0 && flowType == 11) { lattice.toggleInternalStatistics(true); }
        global::timer("LBM").start();
        lattice.collideAndStream();
        global::timer("LBM").stop();
        for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
            // #4# IBM Interpolation
            cellFields[iCell]->interpolateVelocityIBM();
            // #5# Position Update
            cellFields[iCell]->advanceParticles();
        }

        // #6# Output
        if ((iter + 1) % tmeas == 0) {
            SyncRequirements everyCCR(allReductions);
            for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                cellFields[iCell]->synchronizeCellQuantities(everyCCR);
            }
            global::timer("HDFOutput").start();
            bool invertXZ_for_XDMF = true;
            writeHDF5(lattice, parameters, iter + 1, invertXZ_for_XDMF);
            for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                writeCellField3D_HDF5(*cellFields[iCell], dx, dt, iter + 1);
                writeCell3D_HDF5(*cellFields[iCell], dx, dt, iter + 1);
            }
            global::timer("HDFOutput").stop();
            if ((iter + 1) % (2 * tmeas) == 0) {
                global::timer("Checkpoint").start();
                checkpointer.save(lattice, cellFields, iter + 1);
                global::timer("Checkpoint").stop();
            }
            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter + 1);
            pcout << "(main) Iteration:" << iter + 1 << "; time " << dtIteration * 1.0 / tmeas;
//            pcout << "; Volume (" << RBCField[0]->getVolume() << ")";
            pcout << std::endl;
        } else {
            for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                cellFields[iCell]->synchronizeCellQuantities();
            }
        }
        if ((iter + 1) % tmeas == 0 && flowType == 11) { lattice.toggleInternalStatistics(false); }
    }
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        delete cellFields[iCell];
    }
    simpleProfiler.writeIteration(tmax + 1);
    pcout << "Simulation finished." << std::endl;
}