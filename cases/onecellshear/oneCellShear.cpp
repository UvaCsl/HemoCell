#include "ficsion.h"
//#include "cellStretching3D.hh"

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor


int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
        return -1;
    }

    plbInit(&argc, &argv);
    global::timer("ficsion_init").start();

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
    plint initIter = 0;
    plint warmup = 0;
    T Re;
    T dx, dt, dm, dNewton;
    T tau, nu_lbm, u_lbm_max;
    T nu_p, rho_p;
    plint nx, ny, nz;
    T lx, ly, lz;
    T shearRate_p, shearRate;
    T shellDensity, eqLengthRatio, radius;
    T k_WLC, k_rep, k_rest, k_elastic, k_bend, k_volume, k_surface, k_shear, k_stretch, eta_m;
    plint minNumOfTriangles;
    plint rbcModel;
    plint ibmKernel, ibmScheme;
    Array<T, 3> eulerAngles;
    plint tmax, tmeas;
    std::string meshFileName;


    // ------------------------- Read in config file ------------------------------------------------

    pcout << "(main) reading config xml..." << std::endl;

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader documentXML(paramXmlFileName);

    // Check if it is a fresh start or a checkpointed run
    std::string firstField = (*(documentXML.getChildren(
            documentXML.getFirstId())[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    if (firstField == "ficsion") { checkpointed = 0; }
    else { checkpointed = 1; }

    XMLreaderProxy document = checkpointed ? documentXML["Checkpoint"]["ficsion"] : documentXML["ficsion"];
   
    document["parameters"]["shearrate"].read(shearRate_p);
    document["parameters"]["warmup"].read(warmup);

    std::vector<T> ea;
    document["parameters"]["eulerAngles"].read(ea);
    if (ea.size() != 3) { ea.resize(3, 0.0); }
    eulerAngles = Array<T,3>(ea[0], ea[1], ea[2]);
    for(int i =0; i < 3; i++)
        eulerAngles[i] *= pi/180.;

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

    document["ibm"]["ibmKernel"].read(ibmKernel);
    document["ibm"]["ibmScheme"].read(ibmScheme);
    document["ibm"]["radius"].read(radius);
    document["ibm"]["minNumOfTriangles"].read(minNumOfTriangles);

    // Set lx, ly, lz --or nx, ny, nz
    lx = 20 * radius;
    ly = 20 * radius;
    lz = 20 * radius;

    document["domain"]["rhoP"].read(rho_p);
    document["domain"]["nuP"].read(nu_p);
    document["domain"]["tau"].read(tau);
    document["domain"]["dx"].read(dx);
    document["domain"]["dt"].read(dt);

    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);

// ---------------------------- Read in geometry (STL) ------------------------------------------------

    Re = (ly * (shearRate_p * (lx*0.5)) ) / nu_p;
    nx = (int)(lx / dx);
    ny = (int)(ly / dx);
    nz = (int)(lz / dx);

    nu_lbm = nu_p * dt / (dx*dx); 
    u_lbm_max = Re * nu_lbm / nx;

    tau = 3.0 * nu_lbm + 0.5;
    dm = rho_p * (dx * dx * dx);
    dNewton = (dm * dx / (dt * dt));
    //kBT = kBT_p / ( dNewton * dx );
    shearRate = shearRate_p * dt;
    
    
    pcout << "(main) dx = " << dx << ", " <<
        "dt = " << dt << ", " <<
        "dm = " << dm << ", " <<
        "dN = " << dNewton << ", " <<
        "shear rate = " << shearRate << ", " <<
        //"kT = " << kBT <<
        std::endl;

    pcout << "(main) tau = " << tau << " Re = " << Re << " u_lb = " << u_lbm_max << " nu_lb = " << nu_lbm << endl;
    pcout << "(main) Re corresponds to u_max = " << (Re * nu_p)/(ny*dx) << " [m/s]" << endl;


    // Because structures use it. Kind of nonsense. TODO: To be changed later.
    plint resolution = (int)(1.0/dx);
    IncomprFlowParam<T> parameters(
            u_lbm_max,
            Re*(resolution/nx),     // Because Palabos calculates Re in a uniform way (which is not a good thing for CFD)
            resolution,
            lx,
            ly,
            lz
    );

    checkParameterSanity(parameters);

    // Debug line
    // pcout << "(main) tau = " << parameters.getTau() << " Re = " << parameters.getRe() << " u_lb = " << u_lbm_max << " nu_lb = " << parameters.getLatticeNu() << endl;

    // ------------------------ Init lattice --------------------------------

    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;

    plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));


    // -------------------------- Define boundary conditions ---------------------

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    lattice.toggleInternalStatistics(false);

    iniLatticeSquareCouette(lattice, parameters, *boundaryCondition, shearRate);


    // ----------------------- Init cell models --------------------------

    pcout << "(main) init cell structures..."  << std::endl;
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    T persistenceLengthFine = 7.5e-9; // In meters
    T eqVolume = 0;
    plint shape = 1;    // shape: Sphere[0], RBC from sphere[1], Cell(defined)[2], RBC from file [3] RBC from Octahedron [4] Sphere from Octahedron [5]
    std::string cellPath = " "; // If particle is loaded from stl file.

    ConstitutiveModel<T, DESCRIPTOR> *cellModel;
    std::vector<CellField3D<T, DESCRIPTOR> *> cellFields;


    
    // ----------------------- Init RBCs ---------------------------------

    pcout << "(main)   init RBC structure..."  << std::endl;

    TriangleBoundary3D<T> Cells = constructMeshElement(shape, radius / dx, minNumOfTriangles, dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = Cells.getMesh();
    MeshMetrics<T> meshmetric(meshElement);    
    meshmetric.write();
    eqVolume = meshmetric.getVolume();
    plint numVerticesPerCell = meshElement.getNumVertices();

    if (rbcModel == 0) {
        pcout << "(main) Using ShapeMemoryModel3D. " << std::endl;
        cellModel = new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
        persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);
    } else if (rbcModel==1) {
        pcout << "(main) Using CellModel3D. " << std::endl;
        cellModel = new CellModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
         persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);
    } else {
        pcout << "(main) Using IntermediateModel3D. " << std::endl;
        cellModel = new IntermediateModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
             persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);
    }

    // Giving artificial hematocrit of 0.01, since it does not matter here
    CellField3D<T, DESCRIPTOR> RBCField(lattice, meshElement, 0.01, cellModel, ibmKernel, "RBC");
    cellFields.push_back(&RBCField);


// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

    FcnCheckpoint<T, DESCRIPTOR> checkpointer(documentXML);
    checkpointer.load(documentXML, lattice, cellFields, initIter);
    if (not checkpointer.wasCheckpointed()) {
        pcout << "(main) initializing"<< std::endl;
        std::vector<Array<T,3> > cellsOrigin;
        cellsOrigin.push_back( Array<T,3>(nx*0.5, ny*0.5, nz*0.5) );
        RBCField.initialize(cellsOrigin);
        checkpointer.save(lattice, cellFields, initIter);
    }


    plint nCells = RBCField.getNumberOfCells_Global();
    pcout << std::endl ;
    pcout << "(main) Volume ratio [x100%]: " << nCells*eqVolume*100.0/(nx*ny*nz) << std::endl;
    pcout << "(main) nCells (global) = " << nCells << ", pid: " << global::mpi().getRank() << std::endl;

//    MultiParticleField3D<LightParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
//                                                        createBoundaryParticleField3D(lattice);

    /* Repulsive force */
    //T k_int = 0.00025, DeltaX=1.0, R=0.75, k=1.5;
    //PowerLawForce<T> PLF(k_int, DeltaX, R, k);


    // --------------------- Warming up the fluid field ---------------------

    if (initIter == 0)
    {
        pcout << "(main) fresh start: warming up fluid domain for "  << warmup << " terations..." << std::endl;
        for (plint itrt = 0; itrt < warmup; ++itrt) { lattice.collideAndStream(); }
    }

    // -------------------------- Initial output --------------------------

    pcout << std::endl << "(main) Starting simulation i=" << initIter  << ", tmeas = " << tmeas << std::endl;

    global::timer("HDFOutput").start();
    writeHDF5(lattice, parameters, 0);
    writeCellField3D_HDF5(RBCField, dx, dt, 0);
    global::timer("HDFOutput").stop();

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, nCells, numVerticesPerCell);


    // ------------------------ Starting main loop --------------------------
    
    global::timer("mainLoop").start();
    //global::profiler().turnOn();

    for (pluint iter=initIter; iter<tmax+1; ++iter) {
        // #1# Membrane Model
        RBCField.applyConstitutiveModel();
        
        // #2# IBM Spreading
        //RBCField.setFluidExternalForce(0); // TODO: Can this line be omitted?
        RBCField.spreadForceIBM();
        
        // #3# LBM
        global::timer("LBM").start();
        lattice.collideAndStream();
        global::timer("LBM").stop();
        
        // #4# IBM Interpolation
        RBCField.interpolateVelocityIBM();
        
        // #5# Position Update
        RBCField.advanceParticles();
        
        // #6# Output
        if ((iter+1)%tmeas==0) {
            SyncRequirements everyCCR(allReductions);
            RBCField.synchronizeCellQuantities(everyCCR);

            global::timer("HDFOutput").start();
            // writeHDF5(lattice, parameters, iter+1);
            writeCellField3D_HDF5(RBCField, dx, dt, iter+1);
            writeCell3D_HDF5(RBCField, dx, dt, iter+1);
            global::timer("HDFOutput").stop();

            if ((iter+1)%(2*tmeas)==0) {
                global::timer("Checkpoint").start();
                checkpointer.save(lattice, cellFields, iter+1);
                global::timer("Checkpoint").stop();
            }

            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter+1);
            //global::profiler().writeReport();
            pcout << "(main) Iteration:" << iter + 1 << "; time "<< dtIteration*1.0/tmeas ;

        } else {
            RBCField.synchronizeCellQuantities();
        }

    }
    RBCField[0]->saveMesh("stretchedCell.stl");
    simpleProfiler.writeIteration(tmax+1);
    //global::profiler().writeReport();
    pcout << "(main) Simulation finished." << std::endl;
}
