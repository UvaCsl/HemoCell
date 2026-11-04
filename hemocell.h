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
#ifndef HEMOCELL_H
#define HEMOCELL_H

//Load Constants
#include "constant_defaults.h"
#include "config.h"

/* CORE libs */
#include "hemoCellFields.h"

/* IO */
#include "loadBalancer.h"
#include "profiler.h"

/* MECHANICS */
#include "cellMechanics.h"
#include "constantConversion.h"

/* Helpers */
#include "preInlet.h"
// #include "leesEdwardsBC.h"

/* Always used palabos functions in case files*/
#ifndef COMPILING_HEMOCELL_LIBRARY
#include "voxelizeDomain.h"
#include "palabos3D.h"
#include "palabos3D.hh"
using namespace hemo;
using namespace plb;
#endif

namespace hemo { 

/*!
 * The HemoCell class contains all the information, data and methods to set up a
 * basic HemoCell simulation.
 *  
 * Most importantly:
 *    - It contains one Palabos FluidField
 *    - It contains one CellFields
 *
 * After everything is set HemoCell::iterate() is used to iterate the
 * simulation.
 */
class HemoCell {
  public:
  /** \brief Switch to indicate if the MPI instance is managed externally.
   * The MPI instance is typically invoked when initialising the HemoCell class
   * through the `plb::plbInit` routine. However, in some scenarios we might want
   * to manually control the MPI environment and do not want HemoCell to setup
   * MPI environment.
   */
  enum class MPIHandle {Internal, External};

  /**
   * Creates an hemocell object
   *
   * @param configFileName the location of the main config file
   *
   * Unfortunately, due to palabos regulations, it is required to pass the
   * commandline arguments
   */
  HemoCell(char * configFileName, int argc, char* argv[]) : HemoCell{ configFileName, argc, argv, HemoCell::MPIHandle::Internal } {};

  /// Creates an HemoCell object and allows to indicate the MPI environment
  /// is managed externally. If so, this constructor avoids plb::plbInit().
  ///
  /// @param mpi_handle the MPIHanndle indicator.
  HemoCell(char * configFileName, int argc, char* argv[], HemoCell::MPIHandle mpi_handle);

  /*
   * Clean up hemocell
   */
  ~HemoCell();
  
  /**
   * @brief Initialise every fluid node at equilibrium.
   * @param rho density [lattice units, typically 1]
   * @param vel velocity [lattice units]
   */
  void latticeEquilibrium(T rho, hemo::Array<T, 3> vel);

  /**
   * @brief Create the cell fields and the particle field on the current
   * lattice. Collective. Call after the lattice is created and initialised,
   * before adding cell types.
   */
  void initializeCellfield();

  /**
   * @brief Add a cell type with constitutive model `Mechanics`.
   *
   * Example: `addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE)`.
   * @param name identifier of the type; its material is read from
   *        `<name>.xml` and its initial positions from `<name>.pos`
   *        (relative to the working directory). Use the same name in the
   *        per-type setters below.
   * @param constructType shape of the reference mesh: RBC_FROM_SPHERE,
   *        ELLIPSOID_FROM_SPHERE (platelet), WBC_SPHERE or MESH_FROM_STL
   *        (constants in constant_defaults.h).
   */
  template<class Mechanics>
  void addCellType(string name, int constructType) {
    HemoCellField * cellfield = cellfields->addCellType(name,constructType);
    Mechanics * mechanics = new Mechanics(*cellfield->materialCfg, *cellfield);
    cellfield->mechanics = mechanics;
    cellfield->statistics();
  }

  /**
   * @brief Select the cell fields written by writeOutput() for a cell type.
   * @param name cell type name (as given to addCellType())
   * @param outputs OUTPUT_* constants from constant_defaults.h, e.g.
   *        {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE}
   */
  void setOutputs(string name, vector<int> outputs);

  bool repulsionEnabled = false;
  bool boundaryRepulsionEnabled = false;
  /**
   * @brief Enable repulsion between vertices of different cells.
   *
   * Force `k * r_c / d` along the connecting line for `d < r_c`. The
   * neighbour search covers adjacent lattice nodes only (cutoffs above one
   * lattice spacing are partially effective).
   * @param repulsionConstant k [lattice force units, used as given]
   * @param repulsionCutoff r_c [micrometre]
   */
  void setRepulsion(T repulsionConstant, T repulsionCutoff);

  // Lees-Edwards boundary condition
  bool leesEdwardsBC = false;
  double * LEcurrentDisplacement;

  /**
   * @brief Evaluate the constitutive model of a cell type every `separation`
   * iterations (default 1). Must be a multiple of the velocity update interval
   * (setParticleVelocityUpdateTimeScaleSeparation()).
   */
  void setMaterialTimeScaleSeparation(string name, unsigned int separation);
  
  //Enable solidify mechanics of a celltype
  void enableSolidifyMechanics(string name) {
    hlog << "(HemoCell) Enabling Solidify Mechanics for " << name << " mechanical model" << endl;
    (*cellfields)[name]->doSolidifyMechanics = true;
  }
  
  /**
   * @brief Interpolate the fluid velocity to the vertices (and synchronise the
   * particle envelopes) every `separation` iterations (default 1). It must
   * divide the material, repulsion and interior-viscosity intervals.
   * Separations > 1 introduce a numerical error.
   */
  void setParticleVelocityUpdateTimeScaleSeparation(unsigned int separation);

  /// Recompute repulsion every `separation` iterations (default 1); a multiple
  /// of the velocity update interval and, if wall repulsion is enabled too,
  /// equal to its interval.
  void setRepulsionTimeScaleSeperation(unsigned int separation);

  /// Run the solidification step every `separation` iterations (default 1).
  void setSolidifyTimeScaleSeperation(unsigned int separation);

  /// Interior viscosity: update the nodes near membranes every `separation`
  /// iterations and redetect all interior nodes (expensive) every
  /// `separation_entire_grid` iterations; both multiples of the velocity
  /// update interval.
  void setInteriorViscosityTimeScaleSeperation(unsigned int separation, unsigned int separation_entire_grid);

  /**
   * @brief Enable repulsion between vertices and wall (solid) nodes.
   * Collective (detects the wall nodes next to fluid).
   * @param boundaryRepulsionConstant k [lattice force units, used as given]
   * @param boundaryRepulsionCutoff r_c [micrometre]
   * @param timestep update interval [iterations]; equal to the particle
   *        repulsion interval if that is enabled too
   */
  void enableBoundaryParticles(T boundaryRepulsionConstant, T boundaryRepulsionCutoff, unsigned int timestep = 1);

  /**
   * @brief Reject initial vertices of a cell type closer than `distance` to a
   * solid node: vertices within floor(distance / dx) lattice nodes of a wall
   * are not placed and their cells are removed. Call before loadParticles().
   * @param distance [micrometre]
   */
  void setInitialMinimumDistanceFromSolid(string name, T distance);

  /// Select the fluid fields written by writeOutput() (OUTPUT_* constants).
  void setFluidOutputs(vector<int> outputs);

  /// Select the CEPAC fields written by writeOutput() (OUTPUT_* constants).
  void setCEPACOutputs(vector<int> outputs);

  /**
   * @brief Set periodicity of the fluid and the particle field together.
   * Call after initializeCellfield() and before loadParticles().
   * @param axis 0 (x), 1 (y) or 2 (z)
   */
  void setSystemPeriodicity(unsigned int axis, bool bePeriodic);

  /**
   * @brief Number of periodic images of a cell distinguishable along `axis`
   * (default 100 for each axis). Periodic copies of a cell get IDs shifted by
   * multiples of the cell count; the shifts per axis are 1, limit[0] and
   * limit[0]*limit[1], and the largest shift times the number of cells must
   * fit in an int (checked before the first step). Lower the limits for runs
   * with more than about 200,000 cells.
   */
  void setSystemPeriodicityLimit(unsigned int axis, int limit);
  
  //Load the particles
private:
  bool loadParticlesIsCalled = false;
public:
  /**
   * @brief Place the cells listed in the `<name>.pos` files (positions in
   * micrometre, rotations in degrees, applied in X, Y, Z order). Collective.
   * Cells that do not fit the domain are removed.
   */
  void loadParticles();

  /// Load the checkpoint given as configuration (see the case's
  /// `cfg->checkpointed` branch). Collective.
  void loadCheckPoint();

  /// Save a checkpoint to `<outputDirectory>/<checkpointDirectory>/`. Collective.
  void saveCheckPoint();

  /// Write output in SI units (default) or in lattice units.
  bool outputInSiUnits = true;

  /// Write the selected cell, fluid and CEPAC fields (HDF5) and the cell CSV
  /// summary for the current iteration. Collective; does not change the
  /// simulation state.
  void writeOutput();

  /**
   * @brief Advance the coupled fluid and cells by one time step. Collective.
   * The lattice force is reset at the end: a case driven by a body force must
   * set it again after every call.
   */
  void iterate();

  /**
   * @brief Act on a caught exit signal (SIGINT, SIGTERM, SIGHUP, SIGQUIT,
   * SIGUSR1, SIGUSR2).
   *
   * Collective: must be called by all ranks. The ranks agree whether any of
   * them caught a signal; if so a checkpoint is saved (unless
   * checkpointOnExitSignal is false) and the program exits with code 1.
   * iterate() calls this every exitSignalCheckInterval iterations.
   */
  void checkExitSignals();

  /// Iterations between collective exit-signal checks in iterate() (0 disables them).
  unsigned int exitSignalCheckInterval = 100;
  /// Save a checkpoint before exiting on a caught signal.
  bool checkpointOnExitSignal = true;

  //Load balancing library functions
  /// Calculate and return the fractional load imbalance 
  T calculateFractionalLoadImbalance();
  
  ///Load balance the domain (only necessary with nAtomic blocks > nMpi processors, also checkpoints
  void doLoadBalance();
  
  ///Restructure the grid, has an optional argument to specify whether a checkpoint from this iteration is available, default is YES!
  void doRestructure(bool checkpoint_avail = true);
  
  ///Initialize the fluid field with the given management, should be done after specifing the pre inlets and before initializing the cellfields
  void initializeLattice(MultiBlockManagement3D const & management);
 
  PreInlet * preInlet = 0;
  
  bool partOfpreInlet = false;
  
  map<plint,plint> BlockToMpi;
  
  LoadBalancer * loadBalancer = 0;
  ///The fluid lattice
  MultiBlockLattice3D<T, DESCRIPTOR> * lattice = 0, *preinlet_lattice = 0, * domain_lattice = 0;
  
  MultiBlockManagement3D * preinlet_lattice_management = 0, * domain_lattice_management = 0;
  
  Config * cfg = 0;
  ///The cellfields contains the particle field and all celltypes
  HemoCellFields * cellfields = 0;
  unsigned int iter = 0;
  
  XMLreader * documentXML = 0; //Needed for legacy checkpoint reading TODO fix

  /// Iteration of the last writeOutput() call. Residence times are advanced by
  /// the number of iterations since then, so it is stored in checkpoints.
  unsigned int getLastOutputIteration() const { return lastOutputAt; }
  /// Restore the last output iteration (used when loading a checkpoint).
  void setLastOutputIteration(unsigned int iteration) { lastOutputAt = iteration; }
  private:
  /// Store the last time (iteration) output occured
  unsigned int lastOutputAt = 0;
  std::chrono::high_resolution_clock::duration lastOutput = std::chrono::high_resolution_clock::duration::zero();
  
  /// To be run right before the first iteration, all checking should move here
  void sanityCheck();
  /// Checked in iteration, do sanity check when not yet done
  bool sanityCheckDone = false;
};
}
#endif // HEMOCELL_H
