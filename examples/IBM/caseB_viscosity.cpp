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

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

void readFicsionXML(XMLreader & documentXML,std::string & caseId, plint & rbcModel, T & shellDensity, T & k_rest,
        T & k_shear, T & k_bend, T & k_stretch, T & k_WLC, T & eqLengthRatio, T & k_rep, T & k_elastic, T & k_volume, T & k_surface, T & eta_m,
        T & rho_p, T & u, plint & flowType, T & Re, T & shearRate, T & stretchForce, Array<T,3> & eulerAngles, T & Re_p, T & N, T & lx, T & ly, T & lz,
        plint & forceToFluid, plint & ibmKernel, plint & ibmScheme, plint & shape, std::string & cellPath, T & radius, T & deflationRatio, pluint & relaxationTime,
        plint & minNumOfTriangles, pluint & tmax, plint & tmeas, T & hct, plint & npar, plint & flowParam, bool & checkpointed)
    {
    T nu_p, tau, dx;
    T dt, nu_lb;
    std::string firstField = (*(documentXML.getChildren( documentXML.getFirstId() )[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    if (firstField=="ficsion") { checkpointed = 0; }
    else { checkpointed = 1; }

    XMLreaderProxy document = checkpointed?documentXML["Checkpoint"]["ficsion"]:documentXML["ficsion"];
    document["caseId"].read(caseId);
    rbcModel = 3; // Rest Model
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
    k_rest = 0.0;
    document["cellModel"]["kShear"].read(k_shear);
    document["cellModel"]["kStretch"].read(k_stretch);
    document["parameters"]["flowType"].read(flowType);
    document["parameters"]["Re"].read(Re);
    Re = 0;
    document["parameters"]["shearRate"].read(shearRate);
    shearRate = 1000;
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
    T nx, ny, nz;
    try {
        document["domain"]["lx"].read(lx);
        nx = 20;
    } catch(const plb::PlbIOException & message) {
        document["domain"]["nx"].read(nx);
    }
    lx = nx * radius;
    ly = nx * radius;
    lz = nx * radius;

    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    try {
        document["sim"]["npar"].read(npar);
        hct = 0;
    } catch(const plb::PlbIOException & message) {
        document["sim"]["hct"].read(hct);
        npar = 0;
    }
    if (hct>1.0) { hct /= 100.0; }

    radius = radius*1.0/dx; // Transform from [m] to [LU]
    nu_lb = (tau-0.5)/3.0;
    dt = (nu_lb/nu_p)*dx*dx;
    u = dt*1.0/dx;
    Re_p = 1.0/nu_p;
    N = int(1.0/dx);
    flowParam = 0;
    tmax = 0.07/dt; // 0.05 seconds.
    tmeas = ceil(tmeas * 9.803921568235293e-08/dt);

    if (minNumOfTriangles <= 66) { shape = 5; minNumOfTriangles = 100; }
    else if (minNumOfTriangles <= 162) { shape = 0; minNumOfTriangles = 100; }
    else if (minNumOfTriangles <= 258) { shape = 5; minNumOfTriangles = 400; }
    else if (minNumOfTriangles <= 642) { shape = 0; minNumOfTriangles = 400; }
    else if (minNumOfTriangles <= 1026) { shape = 5; minNumOfTriangles = 800; }
    else if (minNumOfTriangles <= 2562) { shape = 0; minNumOfTriangles = 2000; }
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
    plb_ofstream ofile(  (outputDir + "viscosity.dat").c_str() );
    ofile << setprecision(10) << endl;

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

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);
    XMLreader document(paramXmlFileName);
    pcout << "(main) reading.." <<std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
            k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface, eta_m,
            rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz,  forceToFluid, ibmKernel, ibmScheme, shape, cellPath, radius, deflationRatio, relaxationTime,
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
    plint extendedEnvelopeWidth = 4;  // Because we might use ibmKernel with with 2.
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    lattice.periodicity().toggleAll(true);
    /*
     * Choose case (Square Poiseuille, Couette etc) *
     */
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;
    lattice.toggleInternalStatistics(false);

    // Drag force specific
    Array<plint,3> forceIds;
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();
    forceIds[2] = lattice.internalStatSubscription().subscribeSum();

    pcout << "(main) Using iniLatticeSquareCouetteMeasureStress. " << std::endl;
    lattice.toggleInternalStatistics(true);
    plint nMomentumExchangeCells=0;
    iniLatticeSquareCouetteMeasureStress<T, DESCRIPTOR>(lattice, parameters, *boundaryCondition, shearRate, forceIds, nMomentumExchangeCells);
    lattice.toggleInternalStatistics(false);
    pcout << "Re_p = " << shearRate * 4*radius*radius / parameters.getLatticeNu() << std::endl;
    util::ValueTracer<T> dr_ConvergeX(1, 20, 1.0e-7);


    /*
     * Initialize model *
     */

    /* CREATE MESH */
//    Array<T,3> eqShapeRotationEulerAngles = eulerAngles;
    if (flowParam == 9) { eulerAngles= Array<T,3>(0.,0.,0.); }
    // Radius in LU
    TriangleBoundary3D<T> Cells = constructMeshElement(shape, radius, cellNumTriangles, dx, cellPath, eulerAngles);
    TriangularSurfaceMesh<T> meshElement = Cells.getMesh();
    MeshMetrics<T> meshmetric(meshElement);    meshmetric.write();
    T eqVolume = meshmetric.getVolume();
    plint numVerticesPerCell = meshElement.getNumVertices();



    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    ConstitutiveModel<T,DESCRIPTOR> *cellModel;
    T persistenceLengthFine = 7.5e-9 ; // In meters
    cellModel = new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
         persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);

    CellField3D<T, DESCRIPTOR> RBCField(lattice, meshElement, hct, cellModel, ibmKernel, "RBC");
    std::vector<CellField3D<T, DESCRIPTOR>* > cellFields;
    cellFields.push_back(&RBCField);
//    CellStretch<T, DESCRIPTOR> cellStretch(RBCField, stretchForceScalar, 0.1);

    FcnCheckpoint<T, DESCRIPTOR> checkpointer(document);
    plint initIter=0;
    checkpointer.load(document, lattice, cellFields, initIter);
    if (not checkpointer.wasCheckpointed()) {
        pcout << "(main) initializing"<< std::endl;
        std::vector<Array<T,3> > cellsOrigin;
        cellsOrigin.push_back( Array<T,3>(nx*0.5+1, ny*0.5+1, nz*0.5+1) );
        RBCField.initialize(cellsOrigin);
        checkpointer.save(lattice, cellFields, initIter);
    }
    RBCField.setParticleUpdateScheme(ibmScheme);


    plint nCells = RBCField.getNumberOfCells_Global();
    pcout << std::endl ;
    pcout << "(main) Hematocrit [x100%]: " << nCells*eqVolume*100.0/(nx*ny*nz) << std::endl;
    pcout << "(main) nCells (global) = " << nCells << ", pid: " << global::mpi().getRank() << std::endl;
    pcout << std::endl << "(main) Starting simulation i=" << initIter  << ", tmeas = " << tmeas << std::endl;

//    MultiParticleField3D<LightParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
//                                                        createBoundaryParticleField3D(lattice);

    /* Repulsive force */
    T k_int = 0.00025, DeltaX=1.0, R=0.75, k=1.5;
    PowerLawForce<T> PLF(k_int, DeltaX, R, k);

    ofile << "# dt " << dt << "; dx = " << dx << std::endl;
    ofile << "# Nx " << parameters.getNx() << "; Ny = " << parameters.getNy() -2.5 << "; Nz = " << parameters.getNz() << std::endl;
    ofile << "# simulationVolume " << nMomentumExchangeCells * ( parameters.getNy()-2.5) << std::endl;
    ofile << "# nMomentumExchangeCells " << nMomentumExchangeCells << std::endl;
    pcout << "# simulationVolume " << nMomentumExchangeCells * ( parameters.getNy()-2.5) << std::endl;
    /*            I/O              */
    global::timer("HDFOutput").start();
    writeHDF5(lattice, parameters, 0);
    writeCellField3D_HDF5(RBCField, dx, dt, 0);
    global::timer("HDFOutput").stop();

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, nCells, numVerticesPerCell);
    /* --------------------------- */
    global::timer("mainLoop").start();
    global::profiler().turnOn();


//#ifdef PLB_DEBUG // Palabos has this bug. It's missing the "envelope-update" is the profiler.
    if (flowType==2) { global::profiler().turnOff(); }
//#endif
    for (pluint iter=initIter; iter<tmax+1; ++iter) {
        // #1# Membrane Model
       RBCField.applyConstitutiveModel();
//       if (flowType==3) { cellStretch.stretch(); }
        // #2# IBM Spreading
        RBCField.setFluidExternalForce(poiseuilleForce);
        RBCField.spreadForceIBM();
        // #3# LBM
        if ((iter+1)%tmeas==0 ) { lattice.toggleInternalStatistics(true); }
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
            writeCell3D_HDF5(RBCField, dx, dt, iter+1);
            if ((iter+1)%(20*tmeas)==0) {
                global::timer("Checkpoint").start();
                checkpointer.save(lattice, cellFields, iter+1);
                global::timer("Checkpoint").stop();
            }
            if ((iter+1)%(100*tmeas)==0) {
                global::timer("HDFOutput").start();
                writeHDF5(lattice, parameters, iter+1);
                writeCellField3D_HDF5(RBCField, dx, dt, iter+1);
                global::timer("HDFOutput").stop();
            }
            T dtIteration = global::timer("mainLoop").stop();
            simpleProfiler.writeIteration(iter+1);
            global::profiler().writeReport();

            T nu_lb = parameters.getLatticeNu();
            T coeff = nu_lb * nMomentumExchangeCells * shearRate; // * nMomentumExchangeCells;
            T nu_app =  lattice.getInternalStatistics().getSum(forceIds[0]) / coeff;
            T volumeFraction = (nu_app - 1.0)*2.0/5.0;
            T simulationVolume = nMomentumExchangeCells * ( parameters.getNy()-2.5);
//            pcout << "(main) Iteration:" << iter + 1 << "; time "<< dtIteration*1.0/tmeas ;
//            pcout << "; volumeFraction " << volumeFraction << ", ";
//            pcout << "; Force (" << RBCField[0]->getForce()[0]  << ", ";
//            pcout << RBCField[0]->getForce()[1] << ", ";
//            pcout << RBCField[0]->getForce()[2] << ") ";
//            pcout << "; ForceNorm (" << RBCField[0]->get3D(CCR_FORCE_NORMALIZED)[0]  << ", ";
//            pcout << RBCField[0]->get3D(CCR_FORCE_NORMALIZED)[1] << ", ";
//            pcout << RBCField[0]->get3D(CCR_FORCE_NORMALIZED)[2] << ") ";
//            pcout << "; Vertex_MAX-MIN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MAX) -  RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MIN) << "";
//            pcout << "; Vertex_MAX " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MAX) << "";
//            pcout << "; Vertex_MIN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MIN) << "";
//            pcout << "; Vertex_MEAN " << RBCField[0]->get1D(CCR_CELL_CENTER_DISTANCE_MEAN) << "";
            pcout << std::scientific << std::setprecision(15);
            pcout << "dt " << (iter+1)*dt << "; nu_app " << nu_app << "; volumeFraction " << volumeFraction;
            ofile  << "dt " << (iter+1)*dt << "; nu_app " << nu_app << "; volumeFraction " << volumeFraction;
            if (volumeFraction > 0.0) {
                T measureRadius = pow( (volumeFraction * simulationVolume) * 3.0 / (4.0 * 3.14159),        1.0/3.0 );
                pcout << "; dr (LU) " << measureRadius -radius;
                ofile << "; dr (LU) " << measureRadius -radius;
                dr_ConvergeX.takeValue(measureRadius -radius, true);
            }
            if (dr_ConvergeX.hasConverged()) {
                pcout << "Converged!" << std::endl;
                exit(0);
            }
            pcout << std::endl;
            ofile << std::endl;
        } else {
            RBCField.synchronizeCellQuantities();
        }
        if ((iter+1)%tmeas==0 ) { lattice.toggleInternalStatistics(false); }
    }
    simpleProfiler.writeIteration(tmax+1);
    global::profiler().writeReport();
    pcout << "Simulation finished." << std::endl;
}


