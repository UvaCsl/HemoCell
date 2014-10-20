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
        plint & forceToFluid, plint & ibmKernel, plint & shape, std::string & cellPath, T & radius, T & deflationRatio, pluint & relaxationTime,
        plint & minNumOfTriangles, pluint & tmax, plint & tmeas, plint & npar, plint & flowParam, bool & checkpointed)
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
        document["domain"]["ny"].read(nz);
        lx = nx * dx;
        ly = ny * dx;
        lz = nz * dx;
    }
    document["sim"]["tmax"].read(tmax);
    document["sim"]["tmeas"].read(tmeas);
    document["sim"]["npar"].read(npar);

    radius = radius*1.0/dx; // Transform from [m] to [LU]
    nu_lb = (tau-0.5)/3.0;
    dt = (nu_lb/nu_p)*dx*dx;
    u = dt*1.0/dx;
    Re_p = 1.0/nu_p;
    N = int(1.0/dx);
    flowParam = flowType/10;
    flowType = flowType%10;
    if ( (flowType == 3) or (flowType == 4) or (flowType == 5) or (flowType == 8) ) { // Cell Stretching Analysis
        shearRate = 0;
    } else if (flowType == 6) { // Tumbling Tank Treading Measurements
        if (shearRate > 0) {
            tmax  = 100.0/(shearRate*dt); // Measurements for 100 Strain rates
            tmeas = 0.02/(shearRate*dt); // 50 Measurements per Strain rates
        }
    }
    if (flowParam == 7) {
        tmax += 20.0/dt; // Add 20 seconds to the total time for Fischer2004.
    }
}



int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
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

    global::timer("simulation").start();
    std::string performanceLogFileName = logOutDir + "plbPerformance.log";
    plb_ofstream performanceLogFile(performanceLogFileName.c_str());

    plint forceToFluid, shape, cellNumTriangles, ibmKernel;
    plint rbcModel;
    std::string caseId;
    std::string cellPath;
    pluint tmax;
    plint tmeas, npar;
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
    pcout << "reading.." <<std::endl;
    readFicsionXML(document, caseId, rbcModel, shellDensity,
            k_rest, k_shear, k_bend, k_stretch, k_WLC, eqLengthRatio, k_rep, k_elastic, k_volume, k_surface, eta_m,
            rho_p, u, flowType, Re, shearRate_p, stretchForce_p, eulerAngles, Re_p, N, lx, ly, lz,  forceToFluid, ibmKernel, shape, cellPath, radius, deflationRatio, relaxationTime,
            cellNumTriangles, tmax, tmeas, npar, flowParam, checkpointed);
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
    pcout << "dx = " << dx << ", " <<
             "dt = " << dt << ", " <<
             "dm = " << dt << ", " <<
             "kT = " << kBT <<
             std::endl;

    /* ------------------ *
     * Initialize Lattice *
     * ------------------ */
    plint extendedEnvelopeWidth = 1;  // Because Guo needs 2-cell neighbor access, but we don't use Guo
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    lattice.periodicity().toggleAll(true);
    lattice.toggleInternalStatistics(false);
    /*
     * Choose case (Square Poiseuille, Couette etc) *
     */
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << std::endl << "Initializing lattice: " << nx << "x" << ny << "x" << nz << ": tau=" << tau << std::endl;
    if (flowType == 0) {
        T L_tmp = parameters.getNy();
        T nu_tmp = parameters.getLatticeNu();
        poiseuilleForce = 8 * (nu_tmp*nu_tmp) * Re / (L_tmp*L_tmp*L_tmp) ;
        iniLatticePoiseuilleWithBodyForce<T, DESCRIPTOR>(lattice, parameters, *boundaryCondition, poiseuilleForce);
    }
    else { iniLatticeSquareCouette(lattice, parameters, *boundaryCondition, shearRate); }


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
    performanceLogFile << "# Nx*Ny*Nz; " << nx * ny * nz << std::endl;
    performanceLogFile << "# Nparticles; " << meshElement.getNumVertices() << std::endl;
    performanceLogFile << "# Nx; Ny; Nz; " << nx << "; " << ny << "; "<< nz << std::endl;
    performanceLogFile << "# Ncells; " << npar << std::endl;
//    performanceLogFile << "# particleEnvelopeWidth; " << particleEnvelopeWidth << std::endl;
#ifdef PLB_MPI_PARALLEL
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    performanceLogFile << "# Nprocs; " << ntasks << std::endl;
#endif



    /* The Maximum length of two vertices should be less than 2.0 LU (or not)*/
    ConstitutiveModel<T,DESCRIPTOR> *cellModel;
    T persistenceLengthFine = 7.5e-9 ; // In meters
    if (rbcModel == 0) {
       cellModel = new ShapeMemoryModel3D<T, DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
            persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);
    } else  { // if (rbcModel == 1) {
       cellModel = new CellModel3D<T,DESCRIPTOR>(shellDensity, k_rest, k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_volume, k_surface, eta_m,
               persistenceLengthFine, eqLengthRatio, dx, dt, dm,meshElement);
    }


    std::string logFileName = logOutDir + "plbCells.log";

    T hct = (npar*eqVolume) * 1.0 / ( (nx-1)*(ny-1)*(nz-1) );
    pcout << "Hematocrit [x100%]: " << hct*100 << std::endl;
    pcout << "mu_0 = " << cellModel->getMembraneShearModulus()*dNewton/dx << std::endl;
    pcout << "K = " << cellModel->getMembraneElasticAreaCompressionModulus()*dNewton/dx << std::endl;
    pcout << "YoungsModulus = " << cellModel->getYoungsModulus()*dNewton/dx << std::endl;
    pcout << "Poisson ratio = " << cellModel->getPoissonRatio() << std::endl;

    /* Deflate if I say so */
//    MorsePotential<T> interCellularForce(dx, numVerticesPerCell, kBT, true);
    MorseAndPowerLawForce<T> interCellularForce(dx, numVerticesPerCell, kBT, true);
//    if (flowType == 2) {
//        applyProcessingFunctional (
//           new ComputeCellCellForces3D<T,DESCRIPTOR> (interCellularForce, 1.1e-6/dx),
//           immersedParticles.getBoundingBox(), particleArg );
//    }

    CellField3D<T, DESCRIPTOR> RBCField(lattice, meshElement, 0.4, cellModel, "RBC");
    std::vector<CellField3D<T, DESCRIPTOR>* > cellFields;
    cellFields.push_back(&RBCField);


    FcnCheckpoint<T, DESCRIPTOR> checkpointer(document);
    global::timer("Checkpoint").restart();
    plint initIter=0;
    checkpointer.load(document, lattice, cellFields, initIter);
    if (not checkpointer.wasCheckpointed()) {
        checkpointer.save(lattice, cellFields, initIter);
    }
    global::timer("Checkpoint").stop();

    if (not checkpointer.wasCheckpointed()) {
        pcout << "initializing"<< std::endl;
        RBCField.initialize();
//        RBCField.grow();
    }


    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > * boundaryParticleField3D =
                                                        createBoundaryParticleField3D(lattice);


    pcout << std::endl << "Starting simulation i=" << initIter << std::endl;
    /*            I/O              */
    global::timer("HDFOutput").restart();
    writeHDF5(lattice, parameters, 0);
    writeCellField3D_HDF5(RBCField, dx, dt, 0);
    global::timer("HDFOutput").stop();


    /* --------------------------- */
    for (pluint iter=initIter; iter<tmax+1; ++iter) {
        // #1# Membrane Model
       RBCField.applyConstitutiveModel();
        // #2# IBM Spreading
        RBCField.setFluidExternalForce(poiseuilleForce);
        RBCField.spreadForceIBM();
        // #3# LBM
        global::timer("LBM").start();
        lattice.collideAndStream();
        global::timer("LBM").stop();
        // #4# IBM Interpolation
        RBCField.interpolateVelocityIBM();
        // #5# Position Update
        RBCField.advanceParticles();
        if ((iter+1)%tmeas==0) {
            pcout << "Iteration:" << iter + 1<< std::endl;
            SyncRequirements everyCCR(allReductions);
            RBCField.synchronizeCellQuantities(everyCCR);
            global::timer("HDFOutput").restart();
            writeHDF5(lattice, parameters, iter+1);
            writeCellField3D_HDF5(RBCField, dx, dt, iter+1);
            global::timer("HDFOutput").stop();

            global::timer("Checkpoint").restart();
            checkpointer.save(lattice, cellFields, iter+1);
            global::timer("Checkpoint").stop();
        } else {
            RBCField.synchronizeCellQuantities();
        }
    }


    pcout << "Simulation finished." << std::endl;
}
