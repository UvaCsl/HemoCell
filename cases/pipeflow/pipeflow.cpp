//-----------------
//Add Definition overrides here

//-----------------
#include "hemocell.h"
#include "rbcHO.h"
#include "pltNOOP.cpp"
#include "meshGeneratingFunctions.h"
// ----------------------- Copy from neighbour ------------------------------------

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


// ---------------------- Read in STL geometry ---------------------------------

void getFlagMatrixFromSTL(std::string meshFileName, plint extendedEnvelopeWidth, plint refDirLength, plint refDir,
                          VoxelizedDomain3D<T> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix, plint blockSize) {
    plint extraLayer = 0;   // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
    plint borderWidth = 1;  // Because the Guo boundary condition acts in a one-cell layer.
    
    // Requirement: margin>=borderWidth.
    plint margin = 1;  // Extra margin of allocated cells around the obstacle.

    TriangleSet<T> *triangleSet = new TriangleSet<T>(meshFileName, DBL);

    DEFscaledMesh<T> *defMesh =
            new DEFscaledMesh<T>(*triangleSet, refDirLength, refDir, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    voxelizedDomain = new VoxelizedDomain3D<T>(
            boundary, voxelFlag::inside, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
    
    // Print out some info
    pcout << "(main) Voxelisation is done. Resulting domain parameters are: " << endl;
    pcout << getMultiBlockInfo(voxelizedDomain->getVoxelMatrix()) << std::endl;


    flagMatrix = new MultiScalarField3D<int>((MultiBlock3D &) voxelizedDomain->getVoxelMatrix());

    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix->getBoundingBox(), 1);
    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix->getBoundingBox(), 1);


	// Since the domain is closed, open up the two ends by copying the slice before it.
    Box3D domainBox = flagMatrix->getBoundingBox();
    plint nx = domainBox.getNx();
    plint ny = domainBox.getNy();
    plint nz = domainBox.getNz();

    Box3D domain(0, 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor<int>(Array<plint, 3>(1, 0, 0)), domain, *flagMatrix);

    domain = Box3D(nx - 2, nx - 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor<int>(Array<plint, 3>(-1, 0, 0)), domain, *flagMatrix);

}


// --------------------- main --------------------------------------------------

int main(int argc, char *argv[]) {

    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
        return -1;
    }

    plbInit(&argc, &argv);

    printHeader();

    global::timer("hemocell_init").start();

    global::directories().setOutputDir("./tmp/");
    global::directories().setLogOutDir("./log/");
    global::directories().setInputDir("./");

    global::IOpolicy().activateParallelIO(true);
    global::IOpolicy().setStlFilesHaveLowerBound(false);
//    global::IOpolicy().setLowerBoundForStlFiles(-1.);

    std::string outputDir = global::directories().getOutputDir();
    std::string inputDir = global::directories().getInputDir();
    std::string logOutDir = global::directories().getLogOutDir();

    mkpath((outputDir + "/hdf5/").c_str(), 0777);
    mkpath(logOutDir.c_str(), 0777);


    // ------------------------ Variable declarations here -----------------------------------------

    bool checkpointed = 0;
    //bool isAdaptive = false;
    plint maxInnerIterSize = 1;
    plint minInnerIterSize = 1;
    const plint probeMaterialForceMinPeriod = 10;
    T minForce, maxForce;
    T ppForceDistance;
    plint initIter = 0;
    plint cellStep;
    plint nx, ny, nz;
    T hematocrit;
    T radius;
    plint minNumOfTriangles;
    plint tmax, tmeas;
    plint refDir, refDirN;
    std::string meshFileName;
    string particlePosFile;


    // ------------------------- Read in config file ------------------------------------------------

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);

    pcout << "(main) reading " <<  paramXmlFileName << "..." << std::endl;

    XMLreader documentXML(paramXmlFileName);
    Config cfg(paramXmlFileName);

    // ---- Particle material properties

    cfg["ibm"]["radius"].read(radius);
    cfg["ibm"]["minNumOfTriangles"].read(minNumOfTriangles);
    cfg["ibm"]["ppForceDistance"].read(ppForceDistance);

    cfg["domain"]["geometry"].read(meshFileName);
;
    cfg["domain"]["refDir"].read(refDir);
    cfg["domain"]["refDirN"].read(refDirN);
    cfg["domain"]["timeStepSize"].read(cellStep);

    cfg["sim"]["tmax"].read(tmax);
    cfg["sim"]["tmeas"].read(tmeas);
    cfg["sim"]["hematocrit"].read(hematocrit);
    cfg["sim"]["particlePosFile"].read(particlePosFile);


    hematocrit /= 100; // put it between [0;1]
    
    // ---------------------------- Check if an adaptive run is requested ---------------------------------
    if (cellStep < 1){
        //isAdaptive = true;
        maxInnerIterSize = abs(cellStep);        
        cfg["domain"]["minForce"].read(minForce);
        cfg["domain"]["maxForce"].read(maxForce);
        cfg["domain"]["minTimeStepSize"].read(minInnerIterSize);
        cellStep = minInnerIterSize;

        pcout << "(main) Adaptive membrane model integration is enabled! Membrane integration step-size in [" << minInnerIterSize << "; " << maxInnerIterSize << "]." << endl;
    }

    // ---------------------------- Read in geometry (STL) ------------------------------------------------

    pcout << "(main) reading geometry stl..." << std::endl;

    plint extendedEnvelopeWidth = 1;  // Depends on the requirements of the ibmKernel. 4 or even 2 might be enough (also depends on dx)

    plb::MultiScalarField3D<int> *flagMatrix = 0;
    VoxelizedDomain3D<T> *voxelizedDomain = 0;
    if (cfg.checkpointed) {
        getFlagMatrixFromSTL(meshFileName, extendedEnvelopeWidth, refDirN, refDir, voxelizedDomain, flagMatrix, cfg["domain"]["blockSize"].read<plint>());
    } else {
        getFlagMatrixFromSTL(meshFileName, extendedEnvelopeWidth, refDirN, refDir, voxelizedDomain, flagMatrix, -1);
    }
    Box3D domainBox = flagMatrix->getBoundingBox();
    
    nx = domainBox.getNx();
    ny = domainBox.getNy();
    nz = domainBox.getNz();

    // ---------------------------- Calc. LBM parameters -------------------------------------------------
    param::lbm_parameters(cfg, domainBox);
    // ------------------------ Init lattice --------------------------------

    pcout << "(main) init lattice structure..."  << std::endl;
    extendedEnvelopeWidth = 1;
        MultiBlockLattice3D<T, DESCRIPTOR> lattice(
            voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
#if HEMOCELL_CFD_DYNAMICS == 1
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));
#elif HEMOCELL_CFD_DYNAMICS == 2
            new GuoExternalForceMRTdynamics<T, DESCRIPTOR>(1.0/tau)); // Use with MRT dynamics!
#endif

    defineDynamics(lattice, *flagMatrix, lattice.getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

    lattice.periodicity().toggleAll(false);
    lattice.periodicity().toggle(0,true);
    lattice.toggleInternalStatistics(false);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 1., Array<T, 3>(0., 0., 0.));


    // ----------------------- Define external driving force ---------------

    T rPipe = refDirN/2.0 ; // -1 for the wall width is not needed, BB nodes seem to be forced as well
    double poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
   
    setExternalVector(lattice, lattice.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    
    
    lattice.initialize();

    extendedEnvelopeWidth = 20;
    // ----------------------- Init cell models --------------------------

    pcout << "(main) init cell structures..."  << std::endl;
    T persistenceLengthFine = 7.5e-9; // In meters

    CellFields3D  cellFields = CellFields3D(lattice,extendedEnvelopeWidth);
    std::vector<TriangularSurfaceMesh<T> *> meshes;
    std::vector<T> eqVolumes;

    //pcout << "(main) IBM kernel: " << HEMOCELL_KERNEL << endl;
    
    // ----------------------- Init RBCs ---------------------------------

    pcout << "(main)   * init RBC structures..."  << std::endl;
    Array<T,3> eulerAngles(0., 0., 0.);

    plint shape = 1;    // shape: Sphere[0], RBC from sphere[1], Cell(defined)[2], RBC from file [3] RBC from Octahedron [4] Sphere from Octahedron [5]
    std::string cellPath = " "; // If particle is loaded from stl file.

    TriangleBoundary3D<T> RBCCells = constructMeshElement(shape, (radius / param::dx) , minNumOfTriangles, param::dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = RBCCells.getMesh();
    meshes.push_back(&meshElement);
    eqVolumes.push_back(MeshMetrics<T>(meshElement).getVolume());

    MeshMetrics<T> meshmetric(meshElement);
    meshmetric.write();
    plint numVerticesPerCell = meshElement.getNumVertices();
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    
    HemoCellField * rbcfield = cellFields.addCellType(meshElement, hematocrit, "RBC");
    
    RbcHO rbcmechanics(cfg,*rbcfield);
    rbcfield->mechanics = &rbcmechanics;
    
    // ----------------------- Init platelets ----------------------------
    // TODO - Collect material properties, don't just derive them from RBC
    // TODO - We can use a numerically cheaper model

    pcout << "(main)   * init PLT structures..."  << std::endl;
    T pltRadius = 1.15e-6 / param::dx;
    T aspectRatio = 1.0 / 2.3;//(2 * pltRadius);
    TriangleBoundary3D<T> PLTCells = constructMeshElement(6, pltRadius, (plint)ceil(minNumOfTriangles/9.0), param::dx, cellPath, eulerAngles, aspectRatio);
    TriangularSurfaceMesh<T> pltMeshElement = PLTCells.getMesh();
    eqVolumes.push_back(MeshMetrics<T>(pltMeshElement).getVolume());
    HemoCellField * pltfield = cellFields.addCellType(pltMeshElement, 0.0025 * hematocrit, "PLT");
    PltNOOP pltmechanics = PltNOOP();
    pltfield->mechanics = &pltmechanics;
    
    
    vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_INPLANE,OUTPUT_FORCE_AREA};
    cellFields[0]->setOutputVariables(outputs); //SET RBC OUTPUT
    vector<int> outputs2 = {OUTPUT_POSITION,OUTPUT_TRIANGLES};

    cellFields[1]->setOutputVariables(outputs2); //SET PLT OUTPUT

    cellFields[0]->statistics();

    // ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

    //FcnCheckpoint<T, DESCRIPTOR> checkpointer(documentXML);

    if (not cfg.checkpointed) {
        pcout << "(main) initializing particle positions from " << particlePosFile << "..." << std::endl;

        //orderedPositionMultipleCellField3D(cellFields);
        //randomPositionMultipleCellField3D(cellFields, hematocrit, dx, maxPackIter);
       readPositionsBloodCellField3D(cellFields, param::dx, particlePosFile.c_str());
       cellFields.syncEnvelopes();
       cellFields.deleteIncompleteCells();
       pcout << "(main) saving checkpoint..." << std::endl;
       cellFields.save(&documentXML, initIter);
       
        pcout << "Initial output ..." << std::endl;
            //Repoint surfaceparticle forces for output
            cellFields.separate_force_vectors();
            
            //Recalculate the forces
            cellFields.applyConstitutiveModel();
            
            writeCellField3D_HDF5(cellFields,param::dx,param::dt,initIter);
            
            //Repoint surfaceparticle forces for speed
            cellFields.unify_force_vectors();
    }
    else {
        cellFields.load(&documentXML, initIter);
    	pcout << "(main) particle positions read from checkpoint." << std::endl;
    }
   
    // ------------------------- Output flow parameters ----------------------
    plint domainVol = computeSum(*flagMatrix);

    /*
    pcout << "(main) material model parameters: model: " << HEMOCELL_MATERIAL_MODEL <<  ", update scheme: " << HEMOCELL_MATERIAL_INTEGRATION << ", step-size: " << cellStep << " lt" << endl;

    pcout << "(main) Number of fluid cells: " << domainVol << " / " << nx*ny*nz << " (" << domainVol * (dx * dx * dx * 1e18) << " um^3)" << std::endl;
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        plint nCells = cellFields[iCell]->getNumberOfCells_Global();

        // Get the dynamic volume from the radius increased by the particle-particle force cutoff
        Array<T,2> xRange, yRange, zRange;
        meshes[iCell]->computeBoundingBox (xRange, yRange, zRange);
        Array<T, 3> cellBounds (xRange[1]-xRange[0], yRange[1]-yRange[0], zRange[1]-zRange[0]);
        Array<T, 3> cellRatio;
        T dynVolume = eqVolumes[iCell];
        for(int iDir = 0; iDir < 3; iDir++)
            dynVolume *= 1.0 + (ppForceDistance * 0.5 * 1e-6 / dx) / cellBounds[iDir];

        pcout << "(main) Species: #" << iCell << " (" << nCells * eqVolumes[iCell] * 100.0 / domainVol << " (" << nCells * dynVolume * 100.0 / domainVol <<")" << std::endl;
        pcout << "(main)   nCells (global) = " << nCells << ", pid: " << global::mpi().getRank();
        pcout << ", Cell volume = " << eqVolumes[iCell] * (dx * dx * dx * 1e18) << " um^3 (" <<  dynVolume * (dx * dx * dx * 1e18) << "  um^3)" << std::endl;
    }
    */
    // ------------------------- Warming up fluid domain ------------------    

    if (initIter == 0)
    {
        pcout << "(main) fresh start: warming up cell-free fluid domain for "  << cfg["parameters"]["warmup"].read<plint>() << " iterations..." << std::endl;
        for (plint itrt = 0; itrt < cfg["parameters"]["warmup"].read<plint>(); ++itrt) { lattice.collideAndStream(); }
    }
	T meanVel = computeSum(*computeVelocityNorm(lattice)) / domainVol;
	pcout << "(main) Mean velocity: " << meanVel  * (param::dx/param::dt) << " m/s; Apparent rel. viscosity: " << (param::u_lbm_max*0.5) / meanVel << std::endl;

    // ------------------------ Starting main loop --------------------------
    pcout << std::endl << "(main) * starting simulation at " << initIter << " of tmax=" << tmax << " iterations (" << tmax * param::dt << " s)..." << std::endl;
    for (plint iter = initIter; iter < tmax; iter++) {
        cellFields.applyConstitutiveModel();    // Calculate Force on Vertices
            
        // ##### Particle Force to Fluid ####
        cellFields.spreadParticleForce();
        
        // ##### 3 ##### LBM
        lattice.collideAndStream();
        
        // Reset Forces on the lattice, TODO do own efficient implementation
        setExternalVector(lattice, lattice.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

        // ##### 4 ##### IBM interpolation
        cellFields.interpolateFluidVelocity();

        //### 6 ### Might be together with interpolation
        cellFields.syncEnvelopes();
        
        //### 7 ### 
        cellFields.advanceParticles();
 
 
        //Output and checkpoint
        if ((iter + 1 % tmeas) == 0) {
            //Repoint surfaceparticle forces for output
            cellFields.separate_force_vectors();
            
            //Recalculate the forces
            cellFields.applyConstitutiveModel();
            
            writeCellField3D_HDF5(cellFields,param::dx,param::dt,iter);
            
            //Repoint surfaceparticle forces for speed
            cellFields.unify_force_vectors();
                    pcout << "(main) saving checkpoint..." << std::endl;

            cellFields.save(&documentXML, iter);
            T meanVel = computeSum(*computeVelocityNorm(lattice)) / domainVol;
            pcout << "(main) Iteration:" << iter << "(" << iter * param::dt << " s)" << std::endl;
            pcout << "(main) Mean velocity: " << meanVel * (param::dx/param::dt) << " m/s; Apparent rel. viscosity: " << (param::u_lbm_max*0.5) / meanVel << std::endl;  
        }

    }
    
    pcout << "(main) * simulation finished. :)" << std::endl;
}
