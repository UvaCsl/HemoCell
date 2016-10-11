#include "ficsion.h"
#include "cellStretching3D.hh"

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
    T Re;
    T dx, dt, dm, dNewton;
    plint cellStep;
    T tau, nu_lbm, u_lbm_max;
    T nu_p, rho_p;
    plint nx, ny, nz;
    T lx, ly, lz;
    T stretchForce_p;
    T stretchForceScalar;
    T shellDensity, eqLengthRatio, radius;
    T k_WLC, k_rep, k_rest, k_elastic, k_bend, k_volume, k_surface, k_shear, k_stretch, eta_m;
    plint minNumOfTriangles;
    plint rbcModel;
    plint materialModel;
    plint ibmKernel, ibmScheme;
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
   
    document["parameters"]["stretchForce"].read(stretchForce_p); // In picoNewton
    stretchForce_p *= 1e-12;  // Change to Newton (SI)

    document["cellModel"]["rbcModel"].read(rbcModel);
    document["cellModel"]["materialModel"].read(materialModel);
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
    document["domain"]["timeStepSize"].read(cellStep);

    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);

// ---------------------------- Calculate parameters  ------------------------------------------------

    Re = 1;
    nx = (int)(lx / dx);
    ny = (int)(ly / dx);
    nz = (int)(lz / dx);

    nu_lbm = nu_p * dt / (dx*dx); 
    u_lbm_max = Re * nu_lbm / nx;

    tau = 3.0 * nu_lbm + 0.5;
    dm = rho_p * (dx * dx * dx);
    dNewton = (dm * dx / (dt * dt));
    stretchForceScalar = stretchForce_p / dNewton;
    //kBT = kBT_p / ( dNewton * dx );
    //shearRate = shearRate_p * dt;
    
    
    pcout << "(main) dx = " << dx << ", " <<
        "dt = " << dt << ", " <<
        "dm = " << dm << ", " <<
        "dN = " << dNewton << ", " <<
        "stretchF = " << stretchForceScalar <<
        //"kT = " << kBT <<
        std::endl;


    pcout << "(main) tau = " << tau << " Re = " << Re << " u_lb = " << u_lbm_max << " nu_lb = " << nu_lbm << endl;
    pcout << "(main) Re corresponds to u_max = " << (Re * nu_p)/(ny*dx) << " [m/s]" << endl;

    // Debug line
    // pcout << "(main) tau = " << parameters.getTau() << " Re = " << parameters.getRe() << " u_lb = " << u_lbm_max << " nu_lb = " << parameters.getLatticeNu() << endl;

    // ------------------------ Init lattice --------------------------------

    plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(1./tau));
    
    lattice.periodicity().toggleAll(true);

    

    // -------------------------- Define initial conditions ---------------------

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;
    lattice.toggleInternalStatistics(false);

    lattice.initialize();

    // Drag force specific
    //T wallVelocity =  0.0; //Re*parameters.getLatticeNu()/( 2 * (radius/dx));
    //pcout << "(main)  wall velocity: " << wallVelocity << endl;
    //iniLattice_ForGalileanInvariance(lattice, parameters, *boundaryCondition, wallVelocity);

    //util::ValueTracer<T> forceConvergeX(1, 20, 1.0e-4);


    // ----------------------- Init cell models --------------------------

    pcout << "(main) init cell structures..."  << std::endl;
    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    T persistenceLengthFine = 7.5e-9; // In meters
    T eqVolume = 0;
    T eqSurface = 0;
    plint shape = 1;    // shape: Sphere[0], RBC from sphere[1], Cell(defined)[2], RBC from file [3] RBC from Octahedron [4] Sphere from Octahedron [5]
    Array<T,3> eulerAngles(0., 0., 0.);
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
    eqSurface = meshmetric.getSurface();
    plint numVerticesPerCell = meshElement.getNumVertices();

    if (rbcModel == 0) {
        pcout << "(main) Using ShapeMemoryModel3D. " << std::endl;
        cellModel = new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
        persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement, materialModel);
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

    // Define stretch force
    CellStretch<T, DESCRIPTOR> stretchedCell(RBCField, stretchForceScalar, 0.1);

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

    // Update integration scheme and time step for surfaceParticles
    RBCField.setParticleUpdateScheme(ibmScheme, (T)cellStep);

    plint nCells = RBCField.getNumberOfCells_Global();
    pcout << std::endl ;
    pcout << "(main) Volume ratio [x100%]: " << nCells*eqVolume*100.0/(nx*ny*nz) << std::endl;
    pcout << "(main) nCells (global) = " << nCells << ", pid: " << global::mpi().getRank() << std::endl;
    pcout << std::endl << "(main) Starting simulation i=" << initIter  << ", tmeas = " << tmeas << std::endl;

//    MultiParticleField3D<LightParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
//                                                        createBoundaryParticleField3D(lattice);

    /* Repulsive force */
    //T k_int = 0.00025, DeltaX=1.0, R=0.75, k=1.5;
    //PowerLawForce<T> PLF(k_int, DeltaX, R, k);


    // -------------------------- Initial output --------------------------

    global::timer("HDFOutput").start();
    writeHDF5(lattice, dx, dt, 0);
    writeCellField3D_HDF5(RBCField, dx, dt, 0);
    global::timer("HDFOutput").stop();

    RBCField[0]->saveMesh("initialCell.stl");

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, nCells, numVerticesPerCell);

    plb_ofstream fOut;
    if(checkpointer.wasCheckpointed())
        fOut.open("stretch.log", std::ofstream::app);
    else
        fOut.open("stretch.log");
    
    Array<T,3> stretch = stretchedCell.measureStretch();
    fOut << 1 << " " << stretch[0] << " " << stretch[1] << " " << stretch[2] << endl;
    
    // ------------------------ Starting main loop --------------------------
    
    global::timer("mainLoop").start();
    //global::profiler().turnOn();

    for (pluint iter=initIter; iter<tmax+1; ++iter) {
        // #1# Membrane Model
        RBCField.applyConstitutiveModel();
        stretchedCell.stretch();
        
        // Inner iteration cycle of fluid
        for(int innerIt = 0; innerIt < cellStep; innerIt++){
        
            // #2# IBM Spreading
            RBCField.setFluidExternalForce(0); // TODO: Can this line be omitted?
            RBCField.spreadForceIBM();
        
            // #3# LBM
            global::timer("LBM").start();
            lattice.collideAndStream();
            global::timer("LBM").stop();
        }
        
        // #4# IBM Interpolation
        RBCField.interpolateVelocityIBM();
        
        // #5# Position Update
        RBCField.advanceParticles();
        
        // #6# Output
        if ((iter%tmeas)==0) {
            SyncRequirements everyCCR(allReductions);
            RBCField.synchronizeCellQuantities(everyCCR);

            global::timer("HDFOutput").start();
            // writeHDF5(lattice, parameters, iter+1);
            writeCellField3D_HDF5(RBCField, dx, dt, iter);
            writeCell3D_HDF5(RBCField, dx, dt, iter);
            global::timer("HDFOutput").stop();

            if ((iter%(2*tmeas))==0) {
                global::timer("Checkpoint").start();
                checkpointer.save(lattice, cellFields, iter);
                global::timer("Checkpoint").stop();
            }

            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter);
            //global::profiler().writeReport();
            pcout << "(main) Iteration:" << iter << "; time "<< dtIteration*1.0/tmeas ;
            //pcout << "; Volume (" << RBCField[0]->getVolume() << ")";

            /*
            if (RBCField.count(0) > 0) {
                pcout << "; Vertex_MAX-MIN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MAX) -  RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MIN) << "";
                pcout << "; Vertex_MAX " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MAX) << "";
                pcout << "; Vertex_MIN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MIN) << "";
                pcout << "; Vertex_MEAN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MEAN) << "";
            }
            */

            Array<T,3> stretch = stretchedCell.measureStretch();
			pcout << "; Stretch (" << stretch[0] <<", " << stretch[1]<<", " << stretch[2]<<"); Volume  " << RBCField[0]->getVolume() / eqVolume * 100.0 << "%; Surface " << RBCField[0]->getSurface() / eqSurface * 100.0 << "%" << std::endl;
            fOut << iter << " " << stretch[0] << " " << stretch[1] << " " << stretch[2] << " " << RBCField[0]->getVolume() / eqVolume * 100.0 << " " << RBCField[0]->getSurface() / eqSurface * 100.0 << endl;
            RBCField[0]->saveMesh("stretchedCell.stl");
        } else {
            RBCField.synchronizeCellQuantities();
        }

    }
    fOut.close();

    RBCField[0]->saveMesh("stretchedCell.stl");
    simpleProfiler.writeIteration(tmax);
    //global::profiler().writeReport();
    pcout << "(main) Simulation finished." << std::endl;
}
