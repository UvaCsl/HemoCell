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
#ifndef HEMOCELLFIELDS_H
#define HEMOCELLFIELDS_H
namespace hemo {
class HemoCellFields;
class bindingFieldHelper;
class InteriorViscosityHelper;
}
#include "hemoCellParticleField.h"
#include "genericFunctions.h"
#include "hemoCellFunctional.h"
#include "hemoCellField.h"
#include "hemoCellParticle.h"
#include "config.h"
#include <unistd.h>

#include "latticeBoltzmann/advectionDiffusionLattices.hh"
#include "multiBlock/multiBlockLattice3D.hh"
#include "offLattice/triangularSurfaceMesh.hh"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "particles/multiParticleField3D.hh"
#include "parallelism/parallelBlockCommunicator3D.h"

namespace hemo {
class HemoCell;

/*!
 * This class can contain many cellTypes, it keeps track of all the particles
 * of all types. The option exists to get exclusive access to a single
 * cellType if necessary. all the particles are stored in a single particlefield
 */
class HemoCellFields {
public:
  ///Default constructor, needs an palabos lattice, envelope width (lbm units), and hemocell reference
  HemoCellFields(plb::MultiBlockLattice3D<T, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth, HemoCell &);

  /*
   * Create the particle field seperately, takes the arguments set in the constructor
   * Is called in the constructor as well
   */
  void createParticleField(plb::SparseBlockStructure3D * sbStructure_ = 0, plb::ThreadAttribution * tAttribution_ = 0);

  void createCEPACfield();

  ///Used to set variables inside the celltypes for correct access, called through createParticleField
  void InitAfterLoadCheckpoint();

public:
  ///Generic Destructor
  ~HemoCellFields();

  ///Add a cell type with a certain mesh; its name also selects `name_.xml` and `name_.pos`
  HemoCellField * addCellType(std::string name_, int constructType);

  ///Easy access to contained celltypes
  HemoCellField * operator[](unsigned int index);

  ///Easy access to contained celltypes
  HemoCellField * operator[](std::string name);

  ///Get the number of celltypes
  unsigned int size();

  /*Checkpoint functions*/
private:
  void copyXMLreader2XMLwriter(plb::XMLreader const & reader, plb::XMLwriter & writer);
  void copyXMLreader2XMLwriter(plb::XMLreaderProxy readerProxy, plb::XMLwriter & writer);

public:
  /// Load a checkpoint, store the current iteration in `iter`. Checks the
  /// build fingerprint written by saveCheckPointDescription() (fatal error on
  /// a mismatch; older checkpoints without one load with a warning).
  /// Collective.
  void load(plb::XMLreader * documentXML, unsigned int & iter, Config * cfg = NULL);
  /// Write the checkpoint data files, after renaming the previous ones
  /// (including checkpoint.xml) to `.old`. checkpoint.xml itself is written
  /// last, by saveCheckPointDescription(), so that a checkpoint interrupted
  /// while saving has no checkpoint.xml instead of a mismatching one.
  /// Collective.
  void save(plb::XMLreader * documentXML, unsigned int iter, Config * cfg = NULL);
  /// Write checkpoint.xml (configuration, iteration, number of cells and the
  /// build fingerprint) under a temporary name and rename it into place.
  /// Call after all checkpoint data files are written. Collective.
  void saveCheckPointDescription(plb::XMLreader * documentXML, unsigned int iter);

private:
  /// Fatal error if the build fingerprint of a checkpoint does not match this
  /// build and setup; warning if the checkpoint has none (older versions).
  void checkCheckPointFingerprint(plb::XMLreader & documentXML);
  /// Load the CEPAC field from a checkpoint directory if it is enabled and present.
  void loadCEPACfield(const std::string & checkpointDirectory);

public:
  ///Legacy Helper function to get the particle field, mostly unused as direct access is available
  plb::MultiParticleField3D<HemoCellParticleField> & getParticleField3D();
  ///Legacy reads in only RBC an PLT from a single pos file
  void readPositionsCellFields(std::string particlePosFile);

  //Class functionals
  ///Advance the particles in an iteration
  void advanceParticles();

  // Find interior lattice points and set omega
  void findInternalParticleGridPoints();

  // Look for lattice nodes near the membrane and update those
  void internalGridPointsMembrane();

  ///Interpolate the velocity of the fluid to the individual particles
  void interpolateFluidVelocity();

  ///Spread the force of all particles over the fluid in this iteration
  void spreadParticleForce();

  /// Separate the force vectors of particles so it becomes clear what the vector for each separate force is
  void separate_force_vectors();

  /// Unify the force vectors of particles to point to a single force
  void unify_force_vectors();

  /// Apply (and calculate) the repulsion force between particles
  void applyRepulsionForce();

  /// Apply (and calculate) the repulsion force between particles and the boundary.
  /// It is added to sv.force_repulsion; with @p resetFirst that field is
  /// zeroed first (needed when particle-particle repulsion, which normally
  /// resets it, is not enabled).
  void applyBoundaryRepulsionForce(bool resetFirst = false);

  /// Limit the repulsion force (sv.force_repulsion, particle and wall
  /// repulsion together) of every particle to FORCE_LIMIT (param::f_limit);
  /// no-op if FORCE_LIMIT is not defined. Local to each block (no
  /// communication); called after a repulsion update.
  void capRepulsionForce();

  /// Delete any incomplete cells on a block
  void deleteIncompleteCells(bool verbose = true);

  /// Delete the cells that lost vertices, on every block that holds one of
  /// their vertices within @p envelope nodes of its bulk (see
  /// HemoCellParticleField::deleteBrokenCells()). Called by HemoCell::iterate()
  /// after the particles moved. Not collective.
  void deleteBrokenCells(plint envelope, bool verbose);

  /**
   * @brief Complete cells for reading, without changing the simulation.
   *
   * Output and statistics need every cell complete in the blocks holding it.
   * While an object of this class exists, every local block holds a copy of
   * its particles, synchronised with syncEnvelopes() and with incomplete
   * cells deleted. The destructor puts back the particles as they were, so
   * that reading does not change the trajectory (the order of the particles,
   * which sets the summation order of force spreading, included). Changes to
   * the particles inside the scope are discarded.
   *
   * Scopes nest: only the outermost one copies and restores. Collective: all
   * ranks must create and destroy it at the same point of the program.
   */
  class CompleteCellsScope {
  public:
    explicit CompleteCellsScope(HemoCellFields & fields);
    ~CompleteCellsScope();
    CompleteCellsScope(const CompleteCellsScope &) = delete;
    CompleteCellsScope & operator=(const CompleteCellsScope &) = delete;

  private:
    HemoCellFields & fields;
  };

  /// Apply the material model of the cells to the particles, updating their force
  void applyConstitutiveModel(bool forced = false);

  /// Sync the particle envelopes between domains
  void syncEnvelopes();

  /// Get particles in a given domain
  void getParticles(std::vector<HemoCellParticle *> & particles, plb::Box3D & domain);

  /// Add particles to local processors
  void addParticles(std::vector<HemoCellParticle> & particles);

  /// Add boundary particles on the fluid-solid boundary
  void populateBoundaryParticles();

  /// Add bindingSites
  void populateBindingSites(plb::Box3D * box = 0);

  /// Delete non local particles (do not delete in envelopesize)
  void deleteNonLocalParticles(int envelope);

  /// Conditionally solidify cells if requested
  void solidifyCells();
  void prepareSolidification();

  /// increment cell residence time
  void updateResidenceTime(unsigned int rtime);

  //Class Variables

  ///the fluid lattice
  plb::MultiBlockLattice3D<T, DESCRIPTOR> * lattice;
  ///A vector specifying the output variables (from const_defaults.h)
  std::vector<int> desiredFluidOutputVariables;

  std::vector<int> desiredCEPACfieldOutputVariables;
  ///Reference to parent
  HemoCell & hemocell;
  ///Vector containing the cellTypes
  std::vector<HemoCellField *> cellFields;
  ///The envelopeSize for the particles
  pluint envelopeSize;
  /// palabos field storing the particles
  plb::MultiParticleField3D<HemoCellParticleField> * immersedParticles = 0;
  /// seperate preinlet and domain pointers whenever necessary
  plb::MultiParticleField3D<HemoCellParticleField> *preinlet_immersedParticles = 0, *domain_immersedParticles = 0;

  /// palabos field for storing the CPAC scalar field if used
  plb::MultiBlockLattice3D<T, CEPAC_DESCRIPTOR> * CEPACfield = 0;

  /// Repulsion cutoff [lattice units], set by HemoCell::setRepulsion()
  T repulsionCutoff = 0.0;
  /// Repulsion constant [lattice force units], set by HemoCell::setRepulsion() (in pN)
  T repulsionConstant = 0.0;
  ///Timescale seperation for repulsion, set through hemocell.h
  pluint repulsionTimescale = 1;

  /// Wall repulsion cutoff [lattice units], set by HemoCell::enableBoundaryParticles()
  T boundaryRepulsionCutoff = 0.0;
  /// Wall repulsion constant [lattice force units], set by HemoCell::enableBoundaryParticles() (in pN)
  T boundaryRepulsionConstant = 0.0;
  ///Timescale seperation for boundary repulsion, set through hemocell.h
  pluint boundaryRepulsionTimescale = 1;

  ///Timescale seperation for the velocity interpolation from the fluid to the particle
  pluint particleVelocityUpdateTimescale = 1;

  pluint solidifyTimescale = 1;

  pluint interiorViscosityTimescale = 1;
  pluint interiorViscosityEntireGridTimescale = 1;

  ///Limit of cycles in a direction (xyz); offsets below must match these defaults
  int periodicity_limit[3] = {100, 100, 100};
  //Internal calculation of offset created by flattening limits
  int periodicity_limit_offset_y = 100;
  int periodicity_limit_offset_z = 10000;

private:
  /// Fill the particle envelopes of all blocks with the particles in the bulk
  /// of their neighbours (normal and enabled periodic overlaps of the particle
  /// field). Called by syncEnvelopes(); collective over all ranks.
  void exchangeEnvelopes();

  std::vector<std::vector<NoInitChar>> sendBuffers, recvBuffers;
  /// Partner ranks (ascending) and per-partner message sizes, data and receive
  /// overlaps of exchangeEnvelopes(), kept between calls to avoid reallocation.
  std::vector<int> envelopeSendRanks, envelopeRecvRanks;
  std::vector<std::vector<int>> envelopeSendSizes, envelopeRecvSizes;
  std::vector<std::vector<NoInitChar>> envelopeSendData, envelopeRecvData;
  std::vector<std::vector<const plb::CommunicationInfo3D *>> envelopeRecvInfos;
  /// Private duplicate of MPI_COMM_WORLD for the particle exchanges in
  /// syncEnvelopes(), so their messages cannot be matched by Palabos or other
  /// HemoCell traffic. Created in the constructor (collective over all ranks),
  /// freed in the destructor.
  MPI_Comm particleComm = MPI_COMM_NULL;
  /// Number of CompleteCellsScope objects that exist (they nest).
  unsigned int completeCellsScopes = 0;

public:
  /**
   * Total number of cells in the simulation, should be constant
   * Is set through hemocell::loadParticles() is added (TODO: also when removed)
   */
  int number_of_cells = 0;
  /// True once number_of_cells was restored from a checkpoint; it is then no
  /// longer re-derived from the .pos files (which may have changed since).
  bool number_of_cells_from_checkpoint = false;
  /// The ID a cell was loaded with (0 .. number_of_cells - 1) from the ID of
  /// one of its periodic images, which grows by a multiple of number_of_cells
  /// at every periodic crossing (HemoCellParticleDataTransfer::getOffset()).
  inline plint base_cell_id(plint wrapped) const {
    return ((wrapped % number_of_cells) + number_of_cells) % number_of_cells;
  }
  unsigned int max_neighbours = 0;

  /// Helpers created on first use through bindingFieldHelper::get() and
  /// InteriorViscosityHelper::get(); owned (and deleted) by this object.
  bindingFieldHelper * bindingHelper = 0;
  InteriorViscosityHelper * interiorViscosityHelper = 0;

  plb::CommunicationStructure3D * large_communicator = 0;

  void calculateCommunicationStructure();

  /*
   * Functionals needed for access of the cellfields
   */
  class HemoSeperateForceVectors : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoSeperateForceVectors * clone() const;
  };
  class HemoUnifyForceVectors : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoUnifyForceVectors * clone() const;
  };
  class HemoSpreadParticleForce : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoSpreadParticleForce * clone() const;
  };
  class HemoFindInternalParticleGridPoints : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoFindInternalParticleGridPoints * clone() const;
  };

  class HemoInternalGridPointsMembrane : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoInternalGridPointsMembrane * clone() const;
  };


  class HemoInterpolateFluidVelocity : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoInterpolateFluidVelocity * clone() const;
  };
  class HemoAdvanceParticles : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoAdvanceParticles * clone() const;
  };
  class HemoApplyConstitutiveModel : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoApplyConstitutiveModel * clone() const;

  public:
    bool forced = false;
  };
  class HemoRepulsionForce : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoRepulsionForce * clone() const;
  };
  class HemoBoundaryRepulsionForce : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoBoundaryRepulsionForce * clone() const;

  public:
    bool resetFirst = false;
  };
  class HemoDeleteIncompleteCells : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoDeleteIncompleteCells * clone() const;

  public:
    bool verbose;
  };
  class HemoGetParticles : public HemoCellFunctional {
    std::vector<HemoCellParticle *> & particles;
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoGetParticles * clone() const;

  public:
    HemoGetParticles(std::vector<HemoCellParticle *> & particles_) : particles(particles_) {}
  };
  class HemoSetParticles : public HemoCellFunctional {
    std::vector<HemoCellParticle> & particles;
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoSetParticles * clone() const;

  public:
    HemoSetParticles(std::vector<HemoCellParticle> & particles_) : particles(particles_) {}
  };
  class HemoPopulateBoundaryParticles : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoPopulateBoundaryParticles * clone() const;
  };
  class HemoPopulateBindingSites : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoPopulateBindingSites * clone() const;
  };
  class HemoDeleteNonLocalParticles : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoDeleteNonLocalParticles * clone() const;

  public:
    int envelopeSize;
    HemoDeleteNonLocalParticles(int envelope_) : envelopeSize(envelope_) {}
  };
  class HemoSolidifyCells : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoSolidifyCells * clone() const;
  };
  class HemoPrepareSolidification : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoPrepareSolidification * clone() const;
  };
  class HemoupdateResidenceTime : public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>);
    HemoupdateResidenceTime * clone() const;

  public:
    unsigned int rtime;
  };
};
}
#endif
