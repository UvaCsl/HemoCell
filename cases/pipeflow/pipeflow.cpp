//-----------------
//Add Definition overrides here

//-----------------
#include "hemocell.h"

typedef double T;

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
                          VoxelizedDomain3D<T> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix) {
    plint extraLayer = 0;   // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
    plint blockSize = -1;   // Zero means: no sparse representation.
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
    T dx, dt, dm, dNewton;
    T tau, nu_lbm, u_lbm_max;
    T nu_p, rho_p;
    T cell_dt;
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
    cfg["domain"]["rhoP"].read(rho_p);
    cfg["domain"]["nuP"].read(nu_p);
    cfg["domain"]["dx"].read(dx);
    cfg["domain"]["dt"].read(dt);
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
    getFlagMatrixFromSTL(meshFileName, extendedEnvelopeWidth, refDirN, refDir, voxelizedDomain, flagMatrix);
    Box3D domainBox = flagMatrix->getBoundingBox();
    
    nx = domainBox.getNx();
    ny = domainBox.getNy();
    nz = domainBox.getNz();

    // ---------------------------- Calc. LBM parameters -------------------------------------------------

    if(dt < 0.0) { // e.g. == -1, set tau = 1 and calc. dt
        tau = 1.0;
        nu_lbm = 1./3. * (tau - 0.5);
        dt = nu_lbm / nu_p * (dx * dx);
        pcout << "(main) Tau is set to unity, and dt is derived from that! (For the fluid, this is the most numerically stable settings.)" << endl;
    }
    else{  // set dt directly and calculate corresponding tau
        nu_lbm = nu_p * dt / (dx*dx); 
        tau = 3.0 * nu_lbm + 0.5;
    }

    u_lbm_max = cfg["parameters"]["Re"].read<double>() * nu_lbm / ny;  // Approximate max. occuring numerical velocity >0.1 should never happen    
    dm = rho_p * (dx * dx * dx);
    dNewton = (dm * dx / (dt * dt));
    //kBT = kBT_p / ( dNewton * dx );
    //shearRate = shearRate_p * dt;
    //stretchForceScalar = stretchForce_p / dNewton;


    pcout << "(main) dx = " << dx << ", " <<
    "dt = " << dt << ", " <<
    "cell_dt = " << cellStep * dt << ", " <<
    "dm = " << dm << ", " <<
    "dN = " << dNewton << std::endl;
    pcout << "(main) tau = " << tau << " Re = " << cfg["parameters"]["Re"].read<double>() << " u_lb_max(based on Re) = " << u_lbm_max << " nu_lb = " << nu_lbm << endl;
    pcout << "(main) Re corresponds to u_max = " << (cfg["parameters"]["Re"].read<double>() * nu_p)/(ny*dx) << " [m/s]" << endl;

    // Do a general check for suspiciosly off values.
    checkParameterSanity(nu_lbm, u_lbm_max);
    

    // ------------------------ Init lattice --------------------------------

    pcout << "(main) init lattice structure..."  << std::endl;

    #if HEMOCELL_CFD_DYNAMICS == 1
        MultiBlockLattice3D<T, DESCRIPTOR> lattice(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/tau));
    #elif HEMOCELL_CFD_DYNAMICS == 2
        MultiBlockLattice3D<T, DESCRIPTOR> lattice(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceMRTdynamics<T, DESCRIPTOR>(1.0/tau)); // Use with MRT dynamics!
    #endif

    defineDynamics(lattice, *flagMatrix, lattice.getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

    lattice.periodicity().toggleAll(true);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 1., Array<T, 3>(0., 0., 0.));


    // ----------------------- Define external driving force ---------------

    T rPipe = refDirN/2.0 ; // -1 for the wall width is not needed, BB nodes seem to be forced as well
    T poiseuilleForce =  8 * nu_lbm * (u_lbm_max * 0.5) / rPipe / rPipe;
    
    // If no saving was set up give some sane default
    if (tmeas == 0) {
        tmeas = plint(ny * ny * 1.0 / (4 * nu_lbm * cfg["parameters"]["Re"].read<double>())) / 40;
    }

    setExternalVector(lattice, lattice.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));


    lattice.initialize();

   
    // ----------------------- Init cell models --------------------------

    pcout << "(main) init cell structures..."  << std::endl;
    T persistenceLengthFine = 7.5e-9; // In meters

    std::vector<ShellModel3D<T> *> cellModels;
    CellFields3D  cellFields = CellFields3D(lattice,extendedEnvelopeWidth);
    std::vector<TriangularSurfaceMesh<T> *> meshes;
    std::vector<T> eqVolumes;

    pcout << "(main) IBM kernel: " << HEMOCELL_KERNEL << endl;
    
    // ----------------------- Init RBCs ---------------------------------

    pcout << "(main)   * init RBC structures..."  << std::endl;
    Array<T,3> eulerAngles(0., 0., 0.);

    plint shape = 1;    // shape: Sphere[0], RBC from sphere[1], Cell(defined)[2], RBC from file [3] RBC from Octahedron [4] Sphere from Octahedron [5]
    std::string cellPath = " "; // If particle is loaded from stl file.

    TriangleBoundary3D<T> RBCCells = constructMeshElement(shape, (radius / dx) , minNumOfTriangles, dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = RBCCells.getMesh();
    meshes.push_back(&meshElement);
    MeshMetrics<T> meshmetric(meshElement);
    meshmetric.write();
    eqVolumes.push_back(meshmetric.getVolume());
    plint numVerticesPerCell = meshElement.getNumVertices();
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    cellModels.push_back(
             ShapeMemoryModel3D::RBCShapeMemoryModel3D(&cfg, dx, dt, dm, meshElement));
    cellFields.addCellType(meshElement, hematocrit, cellModels[0], "RBC");

    
    // ----------------------- Init platelets ----------------------------
    // TODO - Collect material properties, don't just derive them from RBC
    // TODO - We can use a numerically cheaper model

    pcout << "(main)   * init PLT structures..."  << std::endl;
    T pltRadius = 1.15e-6 / dx;
    T aspectRatio = 1.0 / 2.3;//(2 * pltRadius);
    TriangleBoundary3D<T> PLTCells = constructMeshElement(6, pltRadius, (plint)ceil(minNumOfTriangles/9.0), dx, cellPath, eulerAngles, aspectRatio);
    TriangularSurfaceMesh<T> pltMeshElement = PLTCells.getMesh();
    meshes.push_back(&pltMeshElement);
    eqVolumes.push_back(MeshMetrics<T>(pltMeshElement).getVolume());
    cellModels.push_back(
            ShapeMemoryModel3D::PlateletShapeMemoryModel3D(&cfg, dx, dt, dm, pltMeshElement));
    cellFields.addCellType(pltMeshElement, 0.0025 * hematocrit,
                                                        cellModels[cellModels.size() - 1], "PLT");


    // ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

    //FcnCheckpoint<T, DESCRIPTOR> checkpointer(documentXML);
    cellFields.load(&documentXML, initIter);

    if (not cfg.checkpointed) {
        pcout << "(main) initializing particle positions from " << particlePosFile << "..." << std::endl;

        //orderedPositionMultipleCellField3D(cellFields);
        //randomPositionMultipleCellField3D(cellFields, hematocrit, dx, maxPackIter);
       readPositionsBloodCellField3D(cellFields, dx, particlePosFile.c_str());
       
        pcout << "(main) saving checkpoint..." << std::endl;
        cellFields.save(&documentXML, initIter);
    }
    else {
    	pcout << "(main) particle positions read from checkpoint." << std::endl;
    }
    writeCellField3D_HDF5(cellFields,dx,dt,initIter);
#if 0
    // ---------------------- Set integration scheme and time step amplification for cell fields ---------------

    //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        cellFields.setParticleUpdateScheme((T)cellStep);
    //}
    
    // ------------------------- Output flow parameters ----------------------
    pcout << "(main) material model parameters: model: " << HEMOCELL_MATERIAL_MODEL <<  ", update scheme: " << HEMOCELL_MATERIAL_INTEGRATION << ", step-size: " << cellStep << " lt" << endl;

    plint domainVol = computeSum(*flagMatrix);
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

    // ------------------------- Sync all quantities ----------------------

    pcout << "(main) synchronising quantities..."  << std::endl;
    SyncRequirements everyCCR(allReductions);
    //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        cellFields.synchronizeCellQuantities(everyCCR);
    //}

    // -------------------------- Initial output --------------------------

    pcout << "(main) saving initial output..."  << std::endl;
    global::timer("HDFOutput").start();
    bool invertXZ_for_XDMF = true;
    //writeHDF5(lattice, dx, dt, initIter, invertXZ_for_XDMF);
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        //writeCellField3D_HDF5(*cellFields[iCell], dx, dt, initIter);
        //writeCell3D_HDF5(*cellFields[iCell], dx, dt, initIter);
    }
    global::timer("HDFOutput").stop();


    // ------------------------- Create boundary particles ----------------
    
    pcout << "(main) creating boundary particles..."  << std::endl;
    MultiParticleField3D<DenseParticleField3D<T, DESCRIPTOR> > *boundaryParticleField3D =
            createBoundaryParticleField3D(lattice);
    //writeParticleField3D_HDF5(*boundaryParticleField3D, dx, dt, 0, "BoundaryParticles");

    
    // ------------------------- Add particle-wall repulsion force -------------
    // F > 0 means repulsion. F < 0 is used for adhesion
    // This represents glycocalyx layer as well
	
	// For D < 32 um decrease this to recover proper apparent viscosity
    T glycocalyxWidth = 0.75; 	// in [um]
    T k_int = 1.0e-12 / dNewton, DeltaX = glycocalyxWidth * 2.0 * 1e-6 / dx, R = glycocalyxWidth * 2.0 * 1e-6 / dx, k=1.5;  // Scale force with simulation units
    if (DeltaX > 2.0) DeltaX = 2.0; if (R > 2.0) R = 2.0;  // Should not have effect further than 2 lu -> might go out of domain envelope
    PowerLawForce<T> repWP(k_int, DeltaX, R, k);

    // ------------------------- Add particle-particle repulsion force ---------

    k_int = 1.0e-12 / dNewton; DeltaX=ppForceDistance * 1e-6 / dx; T R2=ppForceDistance * 1e-6 / dx; k=1.5;
    if (DeltaX > 2.0) DeltaX = 2.0; if (R2 > 2.0) R2 = 2.0;
    PowerLawForce<T> repPP(k_int, DeltaX, R2, k);


    // ------------------------- Warming up fluid domain ------------------    

    if (initIter == 0)
    {
        pcout << "(main) fresh start: warming up cell-free fluid domain for "  << cfg["parameters"]["warmup"].read<plint>() << " iterations..." << std::endl;
        for (plint itrt = 0; itrt < cfg["parameters"]["warmup"].read<plint>(); ++itrt) { cellFields.setFluidExternalForce(poiseuilleForce); lattice.collideAndStream(); }
    }
	T meanVel = computeSum(*computeVelocityNorm(lattice)) / domainVol;
	pcout << "(main) Mean velocity: " << meanVel  * (dx/dt) << " m/s; Apparent rel. viscosity: " << (u_lbm_max*0.5) / meanVel << std::endl;


    // ------------------------ Starting main loop --------------------------
    pcout << std::endl << "(main) * starting simulation at " << initIter << " of tmax=" << tmax << " iterations (" << tmax * dt << " s)..." << std::endl;
    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, -1, numVerticesPerCell);

    pluint lastCellUpdateSince = 0; // Counts when we need to advance the material model
    T fMax = 0;  // maximal force encountered
    plint adaptiveScaleStep = ceil(maxInnerIterSize/10.);
    plint lastForceProbeSince = 0; // If time-step is small, dont probe force at every step

    
    global::timer("mainLoop").start();
    for (plint iter = initIter; iter < tmax + 1; iter++) {

		// Check if we need to advance the material model during this iteration
        bool stepCells = false; // Whether to advance the material model
        if(0 == (lastCellUpdateSince % cellStep)){
            stepCells = true; lastCellUpdateSince =0;
        }; lastCellUpdateSince++;

        // Check if it is time to probe the forces in the system and allow some adaptivity
        bool forceProbe = false;
        if(stepCells && (lastForceProbeSince > probeMaterialForceMinPeriod)) {
        	forceProbe = true; lastForceProbeSince = 0;
        }; lastForceProbeSince++;
    
    	// Reduce later code complexity, and calculate how much iterations left before saving
    	plint stepsToSave = (iter % tmeas);
    	if (stepsToSave > 0) stepsToSave = tmeas - stepsToSave;

		// ##### 1 ##### Membrane Model + Repulsion of surfaces
        if(stepCells) {
                cellFields.applyConstitutiveModel();    // This zeroes forces at the beginning, always calculate this first
            
           // // Particle-particle force between RBCs and PLTs
            //cellFields[1]->applyDifferentCellForce(repPP, R2, &cellFields[0]->getParticleField3D()); // TODO: it might not sync the argument cell field
            // TODO: extend this in case of other cellFields (e.g. WBC)
	    	//cellFields[0]->syncParticleFieldBulk(); // in case we need a sync point
	    	// TODO: Fix different cell collision and collect them into a single sync
        }

        
        // ##### 2 ##### IBM Spreading
        cellFields.setFluidExternalForce(poiseuilleForce);
        //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
            cellFields.spreadForceIBM();
        //}

        // ##### 3 ##### LBM
        global::timer("LBM").start();
        lattice.collideAndStream();
        global::timer("LBM").stop();

        // ##### 4 & 5 ##### IBM interpolation + Particle position update
        if(stepCells){
            // Gather velocities from the fluid
            //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                //  IBM Interpolation
                cellFields.interpolateVelocityIBM();
                //  Position Update
                cellFields.advanceParticles();
            //}
        }

 		// ###### 6 ###### Probe maximal force & allow time-step adaptivity
        if(forceProbe) {

            fMax = 0;
            for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                T fMax_field = cellFields.getMaximumForce_Global();
                if( fMax_field > fMax)
                    fMax = fMax_field;                
            }

            // If adaptivity enabled, check if we need to modify time-steps
            if(isAdaptive) {

                // If force is too large decrease time step
                if(fMax > maxForce) 
                {
                    if(cellStep > minInnerIterSize){
                        cellStep-=adaptiveScaleStep; if(cellStep < minInnerIterSize) cellStep = minInnerIterSize;

                        pcout << "(main) Large force encountered (" << fMax << ", it: " << iter <<"): reducing inner time-step to: " << cellStep << endl;

                     //   for (pluint iCell = 0; iCell < cellFields.size(); ++iCell)
                            cellFields.setParticleUpdateScheme((T)cellStep);   

                    }
                } 
                else if((fMax < minForce)) // try to increase time-step
                {
                    if(cellStep < maxInnerIterSize){  // Don't go over maxInnerIterSize even if forces are small!
                        cellStep+=adaptiveScaleStep; if(cellStep > maxInnerIterSize) cellStep = maxInnerIterSize;

                        pcout << "(main) Forces are small (" << fMax << ", it: " << iter <<"): increasing inner time-step to: " << cellStep << endl;

                        //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell)
                            cellFields.setParticleUpdateScheme((T)cellStep);

                    }
                }

                // Check if we need to slow down to have an update point at saving
                if(stepsToSave != cellStep && stepsToSave != 0){ 	// If saving will NOT fall to cell update iteration
                	if(stepsToSave < cellStep)	// If we would over-step it
                	{
                		cellStep = stepsToSave;
            	    	pcout << "(main) Inner-step was reset to coincide with checkpointing (" << cellStep << ", it: " << iter <<")" <<endl;
        	        	//for (pluint iCell = 0; iCell < cellFields.size(); ++iCell)
    	                        cellFields.setParticleUpdateScheme((T)cellStep);
	                }
                	else if (stepsToSave < (probeMaterialForceMinPeriod-lastForceProbeSince)) {  // If the next check is too late and we would over-step it (can only happen, if cellStep < probeMaterialForceMinPeriod)
            	    	cellStep = 1;
        	        	pcout << "(main) Inner-step was reset to coincide with checkpointing (" << cellStep << ", it: " << iter <<")" <<endl;
    	            	//for (pluint iCell = 0; iCell < cellFields.size(); ++iCell)
	                            cellFields.setParticleUpdateScheme((T)cellStep);
                	}
            	}
            }
        }

        // #### 7 ##### Output & Sync particles
        if (stepsToSave == 0) {
            SyncRequirements everyCCR(allReductions);
            //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                cellFields.synchronizeCellQuantities(everyCCR);
            //}
            global::timer("HDFOutput").start();
            bool invertXZ_for_XDMF = true;
            //writeHDF5(lattice, dx, dt, iter, invertXZ_for_XDMF);
            //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
            //    writeCellField3D_HDF5(*cellFields[iCell], dx, dt, iter);
            //    writeCell3D_HDF5(*cellFields[iCell], dx, dt, iter);
            //}
            //global::timer("HDFOutput").stop();
            if ((iter % (2 * tmeas)) == 0) {
                cellFields.save(&documentXML, iter);
            }
            T meanVel = computeSum(*computeVelocityNorm(lattice)) / domainVol;
            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter * cellStep);
            pcout << "(main) Iteration:" << iter << "(" << iter * dt << " s)" << "; Wall time / iter. = " << dtIteration / tmeas << " s; Last largest force [lm*lu/lt^2] = " << fMax <<  std::endl;
          	pcout << "(main) Mean velocity: " << meanVel * (dx/dt) << " m/s; Apparent rel. viscosity: " << (u_lbm_max*0.5) / meanVel << std::endl;  
        } else {
            if(stepCells) // Sync only if we changed anything
            //for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
                cellFields.synchronizeCellQuantities(everyCCR);
            //}
        }

    }

    // Be nice and clean up the cell fields
    for (pluint iCell = 0; iCell < cellFields.size(); ++iCell) {
        delete cellFields[iCell];
    }

    simpleProfiler.writeIteration(tmax + 1);

    pcout << "(main) * simulation finished. :)" << std::endl;
#endif
}
