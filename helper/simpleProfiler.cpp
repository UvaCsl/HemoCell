#ifndef FCN_SIMPLE_PROFILER_CPP
#define FCN_SIMPLE_PROFILER_CPP

#include "palabos3D.h"
#include "palabos3D.hh"
#include <string>
using namespace std;
using namespace plb;

class SimpleFicsionProfiler {
public:
    SimpleFicsionProfiler(plint interval_):
        interval(interval_)
    {
        std::string logOutDir = global::directories().getLogOutDir();
        std::string performanceLogFileName = logOutDir + "fcnPerformance.log";
        performanceLogFile.open(performanceLogFileName.c_str());
    }
    ~SimpleFicsionProfiler() { };

    void writeInitial(plint nx, plint ny, plint nz, plint nCells, plint nVerticesPerCell) {
        performanceLogFile << "## mainLoop, LBM, IBM, Quantities and Model are calculated per iteration ##" << std::endl;
        performanceLogFile << "# Nx*Ny*Nz; " << nx * ny * nz << std::endl;
        performanceLogFile << "# Ncells; " << nCells << std::endl;
        performanceLogFile << "# Nparticles; " << nCells * nVerticesPerCell << std::endl;
        performanceLogFile << "# Nx; Ny; Nz; " << nx << "; " << ny << "; "<< nz << std::endl;
//        performanceLogFile << "# particleEnvelopeWidth; " << particleEnvelopeWidth << std::endl;
        performanceLogFile << "# Nprocs; " << global::mpi().getSize() << std::endl;

        double dtIteration = global::timer("ficsion_init").stop(); global::timer("ficsion_init").reset();
        performanceLogFile << "ficsion_init" << "; " << 0 << "; "<< dtIteration << std::endl;
        dtIteration = global::timer("CellInit").stop(); global::timer("CellInit").reset();
        performanceLogFile << "CellInit" << "; " << 0 << "; "<< dtIteration << std::endl;
        dtIteration = global::timer("Checkpoint").stop(); global::timer("Checkpoint").reset();
        performanceLogFile << "Checkpoint" << "; " << 0 << "; "<< dtIteration << std::endl;

    }
    void writeIteration(plint iter) {
        // WRITE PERFORMANCE OUTPUT
        double dtIteration = global::timer("mainLoop").stop(); global::timer("mainLoop").reset();
        performanceLogFile << "Iteration" << "; " << iter << "; "<< dtIteration*1.0/interval << std::endl;
        dtIteration = global::timer("LBM").stop(); global::timer("LBM").reset();
        performanceLogFile << "LBM" << "; " << iter << "; "<< dtIteration*1.0/interval  << std::endl;
        dtIteration = global::timer("IBM").stop(); global::timer("IBM").reset();
        performanceLogFile << "IBM" << "; " << iter << "; "<< dtIteration*1.0/interval  << std::endl;
        dtIteration = global::timer("Quantities").stop(); global::timer("Quantities").reset();
        performanceLogFile << "Quantities" << "; " << iter << "; "<< dtIteration*1.0/interval  << std::endl;
        dtIteration = global::timer("Model").stop(); global::timer("Model").reset();
        performanceLogFile << "Model" << "; " << iter << "; "<< dtIteration*1.0/interval  << std::endl;
        dtIteration = global::timer("CellCellForce").stop(); global::timer("CellCellForce").reset();
        performanceLogFile << "CellCellForce" << "; " << iter << "; "<< dtIteration*1.0/interval  << std::endl;
        dtIteration = global::timer("HDFOutput").stop(); global::timer("HDFOutput").reset();
        performanceLogFile << "HDFOutput" << "; " << iter << "; "<< dtIteration << std::endl;
        dtIteration = global::timer("Checkpoint").stop(); global::timer("Checkpoint").reset();
        if (dtIteration>0) {
            performanceLogFile << "Checkpoint" << "; " << iter << "; "<< dtIteration << std::endl;
        }
        global::timer("mainLoop").start();
    }
private:
    plb_ofstream performanceLogFile;
    plint interval;
};


#endif // FCN_SIMPLE_PROFILER_H
