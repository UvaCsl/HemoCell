/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "hemocell.h"
#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <limits>
#include <mpi.h>

#include "readPositionsBloodCells.h"
#include "hemoCellFunctional.h"
#include "hemoCellParticle.h"
#include "hemoCellField.h"
#include "ParticleHdf5IO.h"
#include "FluidHdf5IO.h"
#include "writeCellInfoCSV.h"
#include "genericFunctions.h"

#include "palabos3D.h"
#include "palabos3D.hh"
#include "preInlet.h"
#include "bindingField.h"
#include "interiorViscosity.h"
#include "cellInfo.h"

using namespace std;
using namespace plb;

using namespace hemo;

namespace {
volatile sig_atomic_t interrupted = 0;
bool interrupt_handler_set = false;

// Signals that request a graceful stop (checkpoint + exit). SIGABRT is
// deliberately not caught: abort() must terminate the process.
const int exit_signals[] = {SIGINT, SIGTERM, SIGHUP, SIGQUIT, SIGUSR1, SIGUSR2};
const int n_exit_signals = sizeof(exit_signals) / sizeof(exit_signals[0]);
struct sigaction previous_actions[n_exit_signals];

void set_interrupt(int signum) {
  interrupted = 1;
}

void set_interrupt_handler() {
  struct sigaction sa;
  memset(&sa, 0, sizeof(struct sigaction));
  sa.sa_handler = set_interrupt;
  sa.sa_flags = SA_RESTART; // do not interrupt system calls (e.g. inside MPI)
  sigemptyset(&sa.sa_mask);
  for (int i = 0; i < n_exit_signals; i++) { sigaction(exit_signals[i], &sa, &previous_actions[i]); }
  interrupt_handler_set = true;
}

void restore_interrupt_handler() {
  if (!interrupt_handler_set) { return; }
  for (int i = 0; i < n_exit_signals; i++) { sigaction(exit_signals[i], &previous_actions[i], 0); }
  interrupt_handler_set = false;
  interrupted = 0;
}

// Give every block of `blocks` one of the ranks firstRank .. firstRank +
// nRanks - 1: in ascending block ID, block i gets firstRank + i when there are
// at most as many blocks as ranks, otherwise consecutive blocks share a rank
// (each rank gets one or more). A block without an entry would get the rank
// equal to its block ID (Palabos' default), which can be a rank of the other
// lattice (with a pre-inlet) or no rank at all.
void attributeBlocksToRanks(const SparseBlockStructure3D & blocks, int firstRank, int nRanks,
                            std::map<plint, plint> & blockToMpi) {
  const plint nBlocks = blocks.getNumBlocks();
  plint i = 0;
  for (const auto & bulk : blocks.getBulks()) {
    blockToMpi[bulk.first] = firstRank + (nBlocks <= nRanks ? i : i * nRanks / nBlocks);
    i++;
  }
}

// Put the process-wide state that a simulation sets back to its initial
// values, so that a next HemoCell object in the same process (as in the test
// suite) starts as it would in a new process. Process-wide state of HemoCell:
// the Parameters (`param::`), `global`, the log file, the signal handlers,
// HemoCellParticleField::cellFields and CellInformationFunctionals::info_per_cell.
// (HemoCellField::default_output is a user setting and is kept; Palabos'
// global settings, such as the output directory, are set by every simulation.)
void resetGlobalState() {
  Parameters::reset();
  resetGlobalConfigValues();
  HemoCellParticleField::cellFields = 0;
  CellInformationFunctionals::clear_list();
  // Give the process its original signal behaviour back.
  restore_interrupt_handler();
  // Close this simulation's log file, so that a next HemoCell object logs to
  // its own file.
  hlog.close();
  hemo::global.hemoCellInitialized = false;
}
} // namespace

HemoCell::HemoCell(const std::string & configFileName, int argc, char * argv[], HemoCell::MPIHandle mpi_handler) {
  if (mpi_handler == HemoCell::MPIHandle::Internal) { plb::plbInit(&(argc), &(argv)); }

  if (global.hemoCellInitialized) { fatal("A HemoCell object already exists; only one may exist at a time."); }
  global.hemoCellInitialized = true;
  // Log (and buffer until the log file is opened) on the main process only.
  hlog.mainProcess = global::mpi().isMainProcessor();

  pcout << "(HemoCell) (Config) reading " << configFileName << endl;
  cfg = new Config(configFileName);
  documentXML = new XMLreader(configFileName);
  if (cfg->checkpointed) {
    pcout
        << "(HemoCell) (Config) Checkpointed config, deferring the loading of the directories (out,log,checkpoint) until loadCheckpoint is called"
        << endl;
  } else {
    loadDirectories(cfg);
  }
  loadGlobalConfigValues(cfg);
  printHeader();

  //Start statistics
  global.statistics.start();
}

HemoCell::~HemoCell() {
  if (cellfields) { delete cellfields; }
  if (cfg) { delete cfg; }
  if (documentXML) {
    delete documentXML;
    documentXML = 0;
  }
  // With a pre-inlet, `lattice` is one of preinlet_lattice and domain_lattice.
  if (preinlet_lattice && preinlet_lattice != lattice) { delete preinlet_lattice; }
  if (domain_lattice && domain_lattice != lattice) { delete domain_lattice; }
  if (lattice) { delete lattice; }
  if (domain_lattice_management) { delete domain_lattice_management; }
  if (loadBalancer) { delete loadBalancer; }
  if (preInlet) { delete preInlet; }
  if (preinlet_lattice_management) { delete preinlet_lattice_management; }
  // Allows another HemoCell object to be created afterwards (e.g. in tests).
  resetGlobalState();
}

void HemoCell::latticeEquilibrium(T rho, hemo::Array<T, 3> vel) {
  hlog << "(HemoCell) (Fluid) Setting Fluid Equilibrium" << endl;
  plb::Array<T, 3> vel_plb = {vel[0], vel[1], vel[2]};
  plb::initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel_plb);
}

void HemoCell::initializeCellfield() {
  if (!domain_lattice) { domain_lattice = lattice; }
  cellfields = new HemoCellFields(*lattice, (*cfg)["domain"]["particleEnvelope"].read<int>(), *this);

  //Set envelope of fluid to 1 again, while maintaining outer one for correct force distribution
  lattice->getMultiBlockManagement().changeEnvelopeWidth(1);
  lattice->signalPeriodicity();

  // Correct place for init
  loadBalancer = new LoadBalancer(*this);
}

void HemoCell::setOutputs(string name, vector<int> outputs) {
  hlog << "(HemoCell) (CellField) Setting output variables for " << name << " cells" << endl;
  vector<int> outputs_c = outputs;
  (*cellfields)[name]->setOutputVariables(outputs_c);
}

void HemoCell::setFluidOutputs(vector<int> outputs) {
  hlog << "(HemoCell) (Fluid) Setting output variables for fluid field" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredFluidOutputVariables = outputs_c;
}

void HemoCell::setCEPACOutputs(vector<int> outputs) {
  hlog << "(HemoCell) (Fluid) Setting CEPAC output variables for fluid field" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredCEPACfieldOutputVariables = outputs_c;
}


void HemoCell::setSystemPeriodicity(unsigned int axis, bool bePeriodic) {
  if (lattice == 0) { fatal("(Periodicity) Create the lattice before setting the periodicity."); }
  if (cellfields == 0 || cellfields->immersedParticles == 0) {
    fatal("(Periodicity) Call initializeCellfield() before setting the periodicity.");
  }
  lattice->periodicity().toggle(axis, bePeriodic);
  cellfields->immersedParticles->periodicity().toggle(axis, bePeriodic);
  cellfields->InitAfterLoadCheckpoint();
}

void HemoCell::setSystemPeriodicityLimit(unsigned int axis, int limit) {
  hlog << "(HemoCell) (Periodicity) Setting periodicity limit of axis " << axis << " to " << limit << endl;
  cellfields->periodicity_limit[axis] = limit;

  //recalculate offsets :
  cellfields->periodicity_limit_offset_y = cellfields->periodicity_limit[0];
  cellfields->periodicity_limit_offset_z = cellfields->periodicity_limit[0] * cellfields->periodicity_limit[1];
}

void HemoCell::loadParticles() {
  hlog << "(HemoCell) (CellField) Loading particle positions " << endl;
  loadParticlesIsCalled = true;
  readPositionsBloodCellField3D(*cellfields, param::dx, *cfg);
  cellfields->syncEnvelopes();
  cellfields->deleteIncompleteCells(false);

  // Report how many cells were actually placed: a cell is dropped as a whole
  // when one of its vertices cannot be placed.
  std::vector<long> placed(cellfields->size(), 0);
  std::vector<plint> ids;
  for (plint bid : cellfields->immersedParticles->getLocalInfo().getBlocks()) {
    for (pluint ctype = 0; ctype < cellfields->size(); ctype++) {
      cellfields->immersedParticles->getComponent(bid).getOutputCells(ctype, ids);
      placed[ctype] += ids.size();
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, placed.data(), int(placed.size()), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  bool dropped = false;
  for (pluint ctype = 0; ctype < cellfields->size(); ctype++) {
    const HemoCellField & type = *(*cellfields)[ctype];
    hlog << "(HemoCell) (CellField) " << type.name << ": placed " << placed[ctype] << " of " << type.cellsInPositionFile
         << " cells" << endl;
    dropped = dropped || placed[ctype] < type.cellsInPositionFile;
  }
  if (dropped) {
    hlog << "(HemoCell) (CellField) WARNING: not all cells were placed. A cell is placed only if all its vertices lie "
            "in the fluid, at least the setInitialMinimumDistanceFromSolid() distance from a wall, and not across a "
            "periodic boundary or outside the domain."
         << endl;
  }
}

void HemoCell::loadCheckPoint() {
  hlog << "(HemoCell) (Saving Functions) Loading Checkpoint" << endl;
  cellfields->load(documentXML, iter, cfg);
  if (global.enableSolidifyMechanics) { bindingFieldHelper::restore(*cellfields); }
  if (global.enableInteriorViscosity) { InteriorViscosityHelper::restore(*cellfields); }
}

void HemoCell::saveCheckPoint() {
  hlog << "(HemoCell) (Saving Functions) Saving Checkpoint at timestep " << iter << endl;
  cellfields->save(documentXML, iter, cfg);
  if (global.enableSolidifyMechanics) { bindingFieldHelper::get(*cellfields).checkpoint(); }
  if (global.enableInteriorViscosity) { InteriorViscosityHelper::get(*cellfields).checkpoint(); }
  // Last, so that an interrupted save leaves no checkpoint.xml (the previous
  // checkpoint is complete in the .old files).
  cellfields->saveCheckPointDescription(documentXML, iter);
}

void HemoCell::writeOutput() {
  global.statistics["output"].start();
  std::string tpi =
      ((iter != lastOutputAt) ? Profiler::toString((global.statistics.elapsed() - lastOutput) / (iter - lastOutputAt))
                              : "0.00");

  //save Residence time
  //Needs to be used before LastOutputAT is updated
  cellfields->updateResidenceTime((iter - lastOutputAt));


  lastOutput = global.statistics.elapsed();
  lastOutputAt = iter;
  // Very naive performance approximation
  pcout << "(HemoCell) (Output) writing output at timestep " << iter << " (" << param::dt * iter
        << " s). Approx. performance: " << tpi << " s / iteration." << endl;
  // Output must not change the simulation, otherwise the output interval
  // would alter the trajectory: it works on a copy of the particles (see
  // CompleteCellsScope), and repulsion is not recomputed (the output shows the
  // repulsion force currently applied, computed at its last update step).

  // update larger lattice envelopes
  lattice->getMultiBlockManagement().changeEnvelopeWidth(2);
  lattice->signalPeriodicity();
  lattice->getBlockCommunicator().duplicateOverlaps(*lattice, modif::staticVariables);
  lattice->getMultiBlockManagement().changeEnvelopeWidth(1);

  {
    // Complete cells in every block, on a copy of the particles that is
    // discarded at the end of this block.
    HemoCellFields::CompleteCellsScope completeCells(*cellfields);

    // Store the total force that will be spread next (force_total) and repoint
    // the per-component force pointers to separate output arrays.
    cellfields->separate_force_vectors();

    // Recompute the force components into the separate arrays only. For cell
    // types with a material timescale separation > 1 the components reflect
    // the current shape, while the total force is the one computed at the last
    // material update.
    cellfields->applyConstitutiveModel(true);

    // Creating a new directory per save
    if (global::mpi().isMainProcessor()) {
      string folder = global::directories().getOutputDir() + "/hdf5/" + zeroPadNumber(iter);
      mkpath(folder.c_str(), 0777);
      folder = global::directories().getOutputDir() + "/csv";
      mkpath(folder.c_str(), 0777);
    }
    global::mpi().barrier();

    // Write Output
    global.statistics.getCurrent()["writeOutput"].start();
    writeCellField3D_HDF5(*cellfields, param::dx, param::dt, iter);
    writeFluidField_HDF5(*cellfields, param::dx, param::dt, iter);
    if (global.enableCEPACfield) { writeCEPACField_HDF5(*cellfields, param::dx, param::dt, iter); }
    writeCellInfo_CSV(*this);
    global.statistics.getCurrent().stop();

    // Free the separate output arrays.
    cellfields->unify_force_vectors();
  }
  // Continue with performance measurement
  global.statistics.getCurrent().stop();
}

void HemoCell::checkExitSignals() {
  if (!interrupt_handler_set) { set_interrupt_handler(); }
  // A signal may reach only some ranks, at different moments. Agree on all
  // ranks before acting, otherwise some ranks would leave while the others
  // block in the next collective operation.
  int local_interrupt = interrupted ? 1 : 0;
  int any_interrupt = 0;
  MPI_Allreduce(&local_interrupt, &any_interrupt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any_interrupt) {
    pcout << endl << "(HemoCell) Caught exit signal at iteration " << iter;
    if (checkpointOnExitSignal) {
      pcout << ", saving a checkpoint and quitting." << endl << std::flush;
      saveCheckPoint();
    } else {
      pcout << ", quitting without checkpoint." << endl << std::flush;
    }
    exit(1);
  }
}

void HemoCell::iterate() {
  if (!interrupt_handler_set) { set_interrupt_handler(); }
  // Collective, so only every exitSignalCheckInterval iterations (same
  // decision on all ranks, since iter is identical everywhere).
  if (exitSignalCheckInterval > 0 && iter % exitSignalCheckInterval == 0) { checkExitSignals(); }
  if (!sanityCheckDone) { sanityCheck(); }
  // (Re)build the large-envelope particle communication structure on the first
  // step and after the particle field has been re-created (load balancing).
  if (!cellfields->large_communicator) { cellfields->calculateCommunicationStructure(); }
  global.statistics.getCurrent()["iterate"].start();
  // ### 1 ### Particle Force to Fluid
  const bool repulsionStep = repulsionEnabled && iter % cellfields->repulsionTimescale == 0;
  const bool boundaryRepulsionStep = boundaryRepulsionEnabled && iter % cellfields->boundaryRepulsionTimescale == 0;
  if (repulsionStep) { cellfields->applyRepulsionForce(); }
  if (boundaryRepulsionStep) {
    // Particle-particle repulsion resets force_repulsion; without it the wall
    // contribution must reset it itself, or it would accumulate every update.
    cellfields->applyBoundaryRepulsionForce(!repulsionEnabled);
  }
  // The repulsion force grows without bound as two points approach; cap it
  // (once per update) like the membrane force in spreadParticleForce().
  if (repulsionStep || boundaryRepulsionStep) { cellfields->capRepulsionForce(); }
  cellfields->spreadParticleForce();

  // #### 2 #### LBM
  global.statistics.getCurrent()["collideAndStream"].start();
  lattice->collideAndStream();
  global.statistics.getCurrent().stop();

  if (global.enableCEPACfield) {
    global.statistics.getCurrent()["CEPACcollideAndStream"].start();
    cellfields->CEPACfield->collideAndStream();
    global.statistics.getCurrent().stop();
  }

  if (iter % cellfields->particleVelocityUpdateTimescale == 0) {
    // #### 3 #### IBM interpolation
    cellfields->interpolateFluidVelocity();
    // ### 4 ### sync the particles
    cellfields->syncEnvelopes();
  }

  if (global.enableSolidifyMechanics && !(iter % cellfields->solidifyTimescale)) {
    global.statistics.getCurrent()["solidifyCells"].start();
    cellfields->prepareSolidification();
    cellfields->syncEnvelopes();
    cellfields->solidifyCells();
    global.statistics.getCurrent().stop();
  }
  // ### 5 ###
  cellfields->advanceParticles();

  // ### 6 ###
  cellfields->applyConstitutiveModel(); // Calculate Force on Vertices

  if (global.enableInteriorViscosity && iter % cellfields->interiorViscosityEntireGridTimescale == 0) {
    cellfields->deleteIncompleteCells(); // Must be done, next function expects whole cells
    global.statistics.getCurrent()["internalParticleGridPoints"].start();
    cellfields->findInternalParticleGridPoints();
    global.statistics.getCurrent().stop();
  }
  if (global.enableInteriorViscosity && iter % cellfields->interiorViscosityTimescale == 0) {
    global.statistics.getCurrent()["internalGridPointsMembrane"].start();
    cellfields->internalGridPointsMembrane();
    global.statistics.getCurrent().stop();
  }

  // We can safely delete non-local cells here, assuming model timestep is divisible by velocity timestep
  if (iter % cellfields->particleVelocityUpdateTimescale == 0) {
    // Cells that lost a vertex (it entered a wall or left the domain, or its
    // force was not finite) are deleted here, at a fixed point of the step, on
    // every block keeping vertices of them. Output and statistics do not
    // delete cells, so that they do not change the trajectory.
    const plint keptEnvelope = 3;
    cellfields->deleteBrokenCells(keptEnvelope, global.cellsDeletedInfo);
    cellfields->deleteNonLocalParticles(keptEnvelope);
  }

  global.statistics.getCurrent()["setExternalVector"].start();
  // Reset Forces on the lattice, TODO do own efficient implementation
  setExternalVector(*lattice, (*lattice).getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  global.statistics.getCurrent().stop();

  iter++;
  global.statistics.getCurrent().stop();
}

T HemoCell::calculateFractionalLoadImbalance() {
  hlog << "(HemoCell) (LoadBalancer) Calculating Fractional Load Imbalance at timestep " << iter << endl;
  return loadBalancer->calculateFractionalLoadImbalance();
}

void HemoCell::setMaterialTimeScaleSeparation(string name, unsigned int separation) {
  hlog << "(HemoCell) (Timescale Seperation) Setting seperation of " << name << " to " << separation << " timesteps"
       << endl;
  (*cellfields)[name]->timescale = separation;
}

void HemoCell::setParticleVelocityUpdateTimeScaleSeparation(unsigned int separation) {
  hlog << "(HemoCell) (Timescale separation) Setting update separation of all particles to " << separation
       << " timesteps" << endl;
  hlogfile << "(HemoCell) WARNING time-scale separation can introduce numerical error! " << endl;
  cellfields->particleVelocityUpdateTimescale = separation;
}

void HemoCell::setRepulsionTimeScaleSeperation(unsigned int separation) {
  hlog << "(HemoCell) (Repulsion Timescale Seperation) Setting seperation to " << separation << " timesteps" << endl;
  cellfields->repulsionTimescale = separation;
}

void HemoCell::setSolidifyTimeScaleSeperation(unsigned int separation) {
  hlog << "(HemoCell) (Solidify Timescale Seperation) Setting seperation to " << separation << " timesteps" << endl;
  cellfields->solidifyTimescale = separation;
}

void HemoCell::setInteriorViscosityTimeScaleSeperation(unsigned int separation, unsigned int separation_entire_grid) {
  hlog << "(HemoCell) (Interior Viscosity Timescale Seperation) Setting seperation to " << separation << " timesteps"
       << endl;
  hlog << "(HemoCell) (Interior Viscosity Timescale Seperation) Setting entire grid refresh to "
       << separation_entire_grid << " timesteps" << endl;

  cellfields->interiorViscosityTimescale = separation;
  cellfields->interiorViscosityEntireGridTimescale = separation_entire_grid;
}

void HemoCell::setInitialMinimumDistanceFromSolid(string name, T distance) {
  hlog << "(HemoCell) (Set Distance) Setting minimum distance from solid to " << distance << " micrometer for " << name
       << endl;
  if (loadParticlesIsCalled) {
    pcout
        << "(HemoCell) (Set Distance) WARNING: this function is called after the particles are loaded, so it has no effect!"
        << endl;
  }
  (*cellfields)[name]->minimumDistanceFromSolid = distance;
}

namespace {
// Repulsion constant [pN] in lattice force units; checked once at setup.
T repulsionConstantToLattice(T constant_pN, T cutoff_um, const char * what) {
  if (!(param::df > 0) || !(param::dx > 0)) {
    fatal("(Repulsion) Compute the lattice parameters (param::lbm_*_parameters()) before enabling ", what, ".");
  }
  if (!(constant_pN >= 0) || !(cutoff_um > 0)) {
    fatal("(Repulsion) The ", what, " constant (", constant_pN, " pN) must be >= 0 and its cutoff (", cutoff_um,
          " um) positive.");
  }
  const T constant = constant_pN * 1e-12 / param::df;
  hlog << "(HemoCell) (Repulsion) " << what << ": k = " << constant_pN << " pN (" << constant
       << " lattice force units), cutoff " << cutoff_um << " um (" << cutoff_um * (1e-6 / param::dx)
       << " lattice units)" << endl;
  return constant;
}
} // namespace

void HemoCell::setRepulsion(T repulsionConstant, T repulsionCutoff) {
  hlogfile << "(HemoCell) (Repulsion) Enabling repulsion." << endl;
  cellfields->repulsionConstant = repulsionConstantToLattice(repulsionConstant, repulsionCutoff, "repulsion");
  cellfields->repulsionCutoff = repulsionCutoff * (1e-6 / param::dx);
  repulsionEnabled = true;
}

void HemoCell::enableBoundaryParticles(T boundaryRepulsionConstant, T boundaryRepulsionCutoff, unsigned int timestep) {
  cellfields->populateBoundaryParticles();
  hlogfile << "(HemoCell) (Repulsion) Enabling boundary repulsion" << endl;
  cellfields->boundaryRepulsionConstant =
      repulsionConstantToLattice(boundaryRepulsionConstant, boundaryRepulsionCutoff, "boundary repulsion");
  cellfields->boundaryRepulsionCutoff = boundaryRepulsionCutoff * (1e-6 / param::dx);
  cellfields->boundaryRepulsionTimescale = timestep;
  boundaryRepulsionEnabled = true;
}

void HemoCell::initializeLattice(MultiBlockManagement3D const & management) {
  if (lattice) { delete lattice; }

  if (!preInlet) {
    try {
      SparseBlockStructure3D sb =
          createRegularDistribution3D(management.getBoundingBox(), (*cfg)["domain"]["mABx"].read<int>(),
                                      (*cfg)["domain"]["mABy"].read<int>(), (*cfg)["domain"]["mABz"].read<int>());

      hlog << "(HemoCell) Domain management overwritten from config file." << endl;
      attributeBlocksToRanks(sb, 0, global::mpi().getSize(), BlockToMpi);
      ExplicitThreadAttribution * eta = new ExplicitThreadAttribution(BlockToMpi);
      domain_lattice_management =
          new MultiBlockManagement3D(sb, eta, management.getEnvelopeWidth(), management.getRefinementLevel());

      lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
          *domain_lattice_management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
          new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

    } catch (const std::invalid_argument & e) {
      hlog << "(HemoCell) Using default domain management." << endl;

      lattice =
          new MultiBlockLattice3D<T, DESCRIPTOR>(management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                 defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                 defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
                                                 new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));
    }

    domain_lattice = lattice;
    return;
  }
  if (global::mpi().getSize() <= 1) { fatal("(PreInlet) A pre-inlet needs at least 2 MPI processes."); }

  plint totalNodes = 0;
  totalNodes += preInlet->getNumberOfNodes();
  totalNodes += cellsInBoundingBox(management.getBoundingBox());

  try { // Look for block management info in the config file
    plint preInlet_pABx = (*cfg)["preInlet"]["parameters"]["pABx"].read<plint>();
    plint preInlet_pABy = (*cfg)["preInlet"]["parameters"]["pABy"].read<plint>();
    plint preInlet_pABz = (*cfg)["preInlet"]["parameters"]["pABz"].read<plint>();
    preInlet->nProcs = preInlet_pABx * preInlet_pABy * preInlet_pABz;
  } catch (const std::invalid_argument & e) {
    preInlet->nProcs = global::mpi().getSize() * (preInlet->getNumberOfNodes() / (T)totalNodes);
    if (preInlet->nProcs == 0) { preInlet->nProcs = 1; }
  }

  const int size = global::mpi().getSize();
  if (preInlet->nProcs >= size) {
    fatal("(PreInlet) The pre-inlet uses ", preInlet->nProcs, " of the ", size,
          " MPI processes; the domain needs at least one more.");
  }
  const int nProcs = size - preInlet->nProcs; // processes of the domain

  // The first preInlet->nProcs ranks hold the pre-inlet, the others the domain.
  partOfpreInlet = global::mpi().getRank() < preInlet->nProcs;

  try { // Look for block management info in the config file
    SparseBlockStructure3D sb_preinlet = createRegularDistribution3D(
        preInlet->location, (*cfg)["preInlet"]["parameters"]["pABx"].read<int>(),
        (*cfg)["preInlet"]["parameters"]["pABy"].read<int>(), (*cfg)["preInlet"]["parameters"]["pABz"].read<int>());
    hlog << "(HemoCell) (PreInlet) Domain management overwritten from config file." << endl;
    attributeBlocksToRanks(sb_preinlet, 0, preInlet->nProcs, preInlet->BlockToMpi);
    ExplicitThreadAttribution * eta_preinlet = new ExplicitThreadAttribution(preInlet->BlockToMpi);
    preinlet_lattice_management = new MultiBlockManagement3D(sb_preinlet, eta_preinlet, management.getEnvelopeWidth(),
                                                             management.getRefinementLevel());
  } catch (const std::invalid_argument & e) { // If nonexistent, fall back to the default distribution
    hlog << "(HemoCell) (PreInlet) Using default domain management." << endl;
    SparseBlockStructure3D sb_preinlet = createRegularDistribution3D(preInlet->location, preInlet->nProcs);
    attributeBlocksToRanks(sb_preinlet, 0, preInlet->nProcs, preInlet->BlockToMpi);
    ExplicitThreadAttribution * eta_preinlet = new ExplicitThreadAttribution(preInlet->BlockToMpi);
    preinlet_lattice_management = new MultiBlockManagement3D(sb_preinlet, eta_preinlet, management.getEnvelopeWidth(),
                                                             management.getRefinementLevel());
  }

  try { // Look for block management info in the config file
    SparseBlockStructure3D sb =
        createRegularDistribution3D(management.getBoundingBox(), (*cfg)["domain"]["mABx"].read<int>(),
                                    (*cfg)["domain"]["mABy"].read<int>(), (*cfg)["domain"]["mABz"].read<int>());
    hlog << "(HemoCell) Domain management overwritten from config file." << endl;
    attributeBlocksToRanks(sb, preInlet->nProcs, nProcs, BlockToMpi);
    ExplicitThreadAttribution * eta = new ExplicitThreadAttribution(BlockToMpi);
    domain_lattice_management =
        new MultiBlockManagement3D(sb, eta, management.getEnvelopeWidth(), management.getRefinementLevel());
  } catch (const std::invalid_argument & e) { // If nonexistent, fall back to the default distribution
    hlog << "(HemoCell) Using default sparse domain management." << endl;
    SparseBlockStructure3D sb = createRegularDistribution3D(management.getBoundingBox(), nProcs);
    attributeBlocksToRanks(sb, preInlet->nProcs, nProcs, BlockToMpi);
    ExplicitThreadAttribution * eta = new ExplicitThreadAttribution(BlockToMpi);
    domain_lattice_management =
        new MultiBlockManagement3D(sb, eta, management.getEnvelopeWidth(), management.getRefinementLevel());
  }

  preinlet_lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      *preinlet_lattice_management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));
  domain_lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      *domain_lattice_management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

  if (!partOfpreInlet) {
    lattice = domain_lattice;
  } else {
    lattice = preinlet_lattice;
    preInlet->partOfpreInlet = true;
  }
}

#ifdef HEMO_PARMETIS
void HemoCell::doLoadBalance() {
  pcout << "(HemoCell) (LoadBalancer) Balancing Atomic Block over mpi processes" << endl;
  loadBalancer->doLoadBalance();
}
#endif

void HemoCell::doRestructure(bool checkpoint_avail) {
  hlog << "(HemoCell) (LoadBalancer) Restructuring Atomic Blocks on processors" << endl;
  loadBalancer->restructureBlocks(checkpoint_avail);
}

void HemoCell::sanityCheck() {
  hlog << "(HemoCell) (SanityCheck) Performing Sanity check on simulation parameters and setup" << endl;
  //Lattice sanity
  if (param::dx != 5e-7) {
    hlog << "(HemoCell) (SanityCheck) WARNING: Fluid dx is not 5e-7 but " << param::dx << " This is unvalidated!"
         << endl;
  }

  int env_min_width = (12e-6 / param::dx) + 1;
  int env_width = cellfields->immersedParticles->getMultiBlockManagement().getEnvelopeWidth();
  if (env_width < env_min_width) {
    hlog << "(HemoCell) (SanityCheck) WARNING: Envelope width is very small: " << env_width << " ("
         << env_width * param::dx << "µm) Instead of " << env_min_width << "!" << endl;
  }

  // Material sanity
  if (boundaryRepulsionEnabled) {
    if (cellfields->boundaryRepulsionTimescale % cellfields->particleVelocityUpdateTimescale != 0) {
      fatal("The boundary repulsion interval (", cellfields->boundaryRepulsionTimescale,
            ") is not a multiple of the particle velocity update interval (",
            cellfields->particleVelocityUpdateTimescale, ").");
    }
  }

  if (repulsionEnabled) {
    if (cellfields->repulsionTimescale % cellfields->particleVelocityUpdateTimescale != 0) {
      fatal("The repulsion interval (", cellfields->repulsionTimescale,
            ") is not a multiple of the particle velocity update interval (",
            cellfields->particleVelocityUpdateTimescale, ").");
    }
  }

  // Both repulsion kinds share sv.force_repulsion, which particle-particle
  // repulsion resets at its updates; with different intervals the wall term
  // would accumulate between those resets.
  if (repulsionEnabled && boundaryRepulsionEnabled &&
      cellfields->repulsionTimescale != cellfields->boundaryRepulsionTimescale) {
    fatal("Particle and boundary repulsion must use the same interval (", cellfields->repulsionTimescale, " vs ",
          cellfields->boundaryRepulsionTimescale, ").");
  }

  // The repulsion neighbour searches only visit adjacent lattice nodes (see
  // computePairRepulsion and applyBoundaryRepulsionForce), so pairs further
  // apart than about one lattice spacing are only partially detected.
  if (repulsionEnabled && cellfields->repulsionCutoff > 1.0) {
    hlog
        << "(HemoCell) (SanityCheck) WARNING: repulsion cutoff (" << cellfields->repulsionCutoff
        << " lu) exceeds one lattice spacing; the neighbour search only covers adjacent lattice nodes, so pairs further apart are partially missed."
        << endl;
  }
  if (boundaryRepulsionEnabled && cellfields->boundaryRepulsionCutoff > 1.0) {
    hlog << "(HemoCell) (SanityCheck) WARNING: boundary repulsion cutoff (" << cellfields->boundaryRepulsionCutoff
         << " lu) exceeds one lattice spacing; only particles in lattice nodes adjacent to a wall node are checked."
         << endl;
  }

  if (global.enableInteriorViscosity) {
    if (cellfields->interiorViscosityEntireGridTimescale % cellfields->particleVelocityUpdateTimescale != 0 ||
        cellfields->interiorViscosityTimescale % cellfields->particleVelocityUpdateTimescale != 0) {
      fatal("The interior viscosity intervals (", cellfields->interiorViscosityTimescale, ", ",
            cellfields->interiorViscosityEntireGridTimescale,
            ") are not multiples of the particle velocity update interval (",
            cellfields->particleVelocityUpdateTimescale, ").");
    }
  }

  for (unsigned int i = 0; i < cellfields->size(); i++) {
    if ((*cellfields)[i]->timescale % cellfields->particleVelocityUpdateTimescale != 0) {
      fatal("The material update interval of ", (*cellfields)[i]->name, " (", (*cellfields)[i]->timescale,
            ") is not a multiple of the particle velocity update interval (",
            cellfields->particleVelocityUpdateTimescale, ").");
    }
  }

  // Cellfields Sanity
  // Check number of neighbours
  if (cellfields->max_neighbours > 30) {
    hlog << "(HemoCell) WARNING: The number of atomic neighbours is suspiciously high: " << cellfields->max_neighbours
         << " Usually it should be < 30 ! Please check the domain decomposition structure!" << endl;
  }

  if (global.enableInteriorViscosity) {
    if (cellfields->interiorViscosityEntireGridTimescale == 1) {
      hlog
          << "(HemoCell) WARNING: Interior viscosity (entire grid) timescale is 1, did you forget to call hemocell->setInteriorViscosityTimescaleSeparation()?"
          << endl;
    }
    if (cellfields->interiorViscosityTimescale == 1) {
      hlog
          << "(HemoCell) WARNING: Interior viscosity timescale is 1, did you forget to call hemocell->setInteriorViscosityTimescaleSeparation()?"
          << endl;
    }
  }

  //Parameter Sanity
#ifdef FORCE_LIMIT
  hlog << "(HemoCell) WARNING: Force limit active at " << FORCE_LIMIT
       << " pN. Results can be inaccurate due to force capping." << endl;
#endif
  if (sizeof(T) == sizeof(float)) {
    hlog << "(HemoCell) WARNING: Running with single precision, you might want to switch to double precision" << endl;
  }

  // Check lattice viscosity [0.01, 0.45]
  if (param::nu_lbm < 0.01 || param::nu_lbm > 0.45) {
    hlog << "(WARNING!!!) lattice viscosity [" << param::nu_lbm << "] is not in the stable range for LBM [0.01, 0.45]!"
         << std::endl;
  }

  // Check for lattice velocity to ensure low Courant number (LBM is explicit afterall...)
  if (param::u_lbm_max > 0.1) {
    hlog << "(WARNING!!!) lattice velocity [" << param::u_lbm_max << "] is too high [>0.1]!" << std::endl;
  }


  if ((*cfg)["sim"]["tmax"].read<long int>() > 100000000000) {
    hlog
        << "(HemoCell) (SanityChecking) More than 100000000000 iterations requested, this means that the zero padding will not be consistent, therefore string sorting output will not work!"
        << endl;
  };

  // A cell ID grows by up to (1 + offset_y + offset_z) x number_of_cells at
  // every periodic crossing and is never folded back
  // (HemoCellParticleDataTransfer::getOffset()); check once that the 64-bit
  // IDs allow plenty of crossings.
  if (cellfields->number_of_cells) { // running with cells
    const plint maxShift =
        (1 + plint(cellfields->periodicity_limit_offset_y) + plint(cellfields->periodicity_limit_offset_z)) *
        plint(cellfields->number_of_cells);
    const plint crossings = std::numeric_limits<plint>::max() / maxShift;
    if (crossings < 1000000) {
      fatal("(SanityCheck) Periodic cell IDs would overflow after ", crossings,
            " periodic crossings of a cell: too many cells (", cellfields->number_of_cells,
            ") for the periodicity limits (setSystemPeriodicityLimit()).");
    }
  }

  bool printed = false;
  plint nx, ny, nz, numBlocks, numAllocatedCells;
  Box3D smallest, largest;
  // Check if atomic block sizes are in an acceptable range, not as fancy as the check in voxelizeDomain, maybe for later
  getMultiBlockInfo(*lattice, nx, ny, nz, numBlocks, smallest, largest, numAllocatedCells);
  // check largest
  if (largest.getNx() > 25 || largest.getNy() > 25 || largest.getNz() > 25) {
    hlog << getMultiBlockInfo(*lattice);
    hlog
        << "(SanityCheck) one of the dimensions of the largest atomic block is more than 25.\n  This can be inefficient. The best performance is below 26x26x26 blocks.\n  It is recommended to adjust the number of processors or the sparseBlockStructure accordingly."
        << endl;
    printed = true;
  }
  if (smallest.getNx() < 16 || smallest.getNy() < 16 || smallest.getNz() < 16) {
    if (!printed) { hlog << getMultiBlockInfo(*lattice); }
    hlog
        << "(SanityCheck) one of the dimensions of the smallest atomic block is less than 16.\n  This can be inefficient. For the best performance we recommend between 20^3 to 25^3 blocks.\n  It is recommended to adjust the number of processors or the sparseBlockStructure, and making sure that the particleEnvelope is smaller than the block size."
        << endl;
    printed = true;
  }


  sanityCheckDone = true;
}