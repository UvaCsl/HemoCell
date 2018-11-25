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
}
#include "hemoCellParticleField.h"
#include "genericFunctions.h"
#include "hemoCellFunctional.h"
#include "hemoCellField.h"
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
class HemoCellFields
{
public:
  
  ///Default constructor, needs an palabos lattice, envelope width (lbm units), and hemocell reference
  HemoCellFields(plb::MultiBlockLattice3D<T, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth,HemoCell &);
 
  /*
   * Create the particle field seperately, takes the arguments set in the constructor
   * Is called in the constructor as well
   */
  void createParticleField(plb::SparseBlockStructure3D* sbStructure_ = 0, plb::ThreadAttribution * tAttribution_ = 0);
  
  void createCEPACfield();
  
  ///Used to set variables inside the celltypes for correct access, called through createParticleField
  void InitAfterLoadCheckpoint();

public:
  ///Generic Destructor
  ~HemoCellFields();
  
  ///Add an celltype with a certain mesh, the name also specifies <name_>.xml and <name_>.pos
  HemoCellField * addCellType(std::string name_, int constructType);
  
  ///Easy access to contained celltypes
  HemoCellField * operator[](unsigned int index);
  
  ///Easy access to contained celltypes
  HemoCellField * operator[](string name);
  
  ///Get the number of celltypes
  unsigned int size();
  
  /*Checkpoint functions*/
private:
  void copyXMLreader2XMLwriter(plb::XMLreader const& reader, plb::XMLwriter & writer);
  void copyXMLreader2XMLwriter(plb::XMLreaderProxy readerProxy, plb::XMLwriter & writer);
public:
  ///Load a checkpoint, store the current iteration in &iter
  void load(plb::XMLreader * documentXML, unsigned int & iter, Config * cfg = NULL);
  ///Save a checkpoint
  void save(plb::XMLreader * documentXML, unsigned int iter, Config * cfg = NULL);
    
  ///Legacy Helper function to get the particle field, mostly unused as direct access is available
  plb::MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & getParticleField3D();
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

  /// Apply (and calculate) the repulsion force between particles and the boundary
  void applyBoundaryRepulsionForce();
  
  /// Delete any incomplete cells on a block
  void deleteIncompleteCells(bool verbose = true);
  
  /// Apply the material model of the cells to the particles, updating their force
  void applyConstitutiveModel(bool forced = false);
  
  /// Sync the particle envelopes between domains
  void syncEnvelopes();

  /// Get particles in a given domain
  void getParticles(vector<HemoCellParticle*> & particles, plb::Box3D & domain);
  
  /// Add particles to local processors
  void addParticles(vector<HemoCellParticle> & particles);

  /// Add boundary particles on the fluid-solid boundary
  void populateBoundaryParticles();
  
  /// Delete non local particles (do not delete in envelopesize)
  void deleteNonLocalParticles(int envelope);
  
  /// Conditionally solidify cells if requested
  void solidifyCells();
  
  //Class Variables
  
  ///the fluid lattice
  plb::MultiBlockLattice3D<T, DESCRIPTOR> * lattice;
  ///A vector specifying the output variables (from const_defaults.h)
  vector<int> desiredFluidOutputVariables;
  
  vector<int> desiredCEPACfieldOutputVariables;
  ///Reference to parent
  HemoCell & hemocell;
  ///Vector containing the cellTypes
  vector<HemoCellField *> cellFields;
  ///The envelopeSize for the particles
  pluint envelopeSize;
  /// palabos field storing the particles
  plb::MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * immersedParticles = 0;
  /// palabos field for storing the CPAC scalar field if used
  plb::MultiBlockLattice3D<T,CEPAC_DESCRIPTOR> * CEPACfield = 0; 

  ///Repulsion variable set through hemocell.h
  T repulsionCutoff = 0.0;
  ///Repulsion variable set through hemocell.h
  T repulsionConstant = 0.0;
  ///Timescale seperation for repulsion, set through hemocell.h
  pluint repulsionTimescale = 1;

  ///Boundary repulsion variable set through hemocell.h
  T boundaryRepulsionCutoff = 0.0;
  ///Boundary repulsion variable set through hemocell.h
  T boundaryRepulsionConstant = 0.0;
  ///Timescale seperation for boundary repulsion, set through hemocell.h
  pluint boundaryRepulsionTimescale = 1;
  
  ///Timescale seperation for the velocity interpolation from the fluid to the particle
  pluint particleVelocityUpdateTimescale = 1;
  
  pluint solidifyTimescale = 1;
  
  pluint interiorViscosityTimescale = 1;
  pluint interiorViscosityEntireGridTimescale = 1;
  
  ///Limit of cycles in a direction (xyz)
  int periodicity_limit[3] = {100};
  //Internal calculation of offset created by flattening limits
  int periodicity_limit_offset_y = 100;
  int periodicity_limit_offset_z = 10000;
  
  /**
   * Total number of cells in the simulation, should be constant
   * Is set through hemocell::loadParticles() is added (TODO: also when removed)
   */
   int number_of_cells = 0;
   inline int base_cell_id(int wrapped) {
     return ((wrapped%number_of_cells)+number_of_cells)%number_of_cells;
   }
   unsigned int max_neighbours = 0;
   
   plb::CommunicationStructure3D * large_communicator = 0;
   plb::ParallelBlockCommunicator3D envelope_communicator;
   
   void calculateCommunicationStructure();
   
   /*
   * Functionals needed for access of the cellfields
   */
  class HemoSeperateForceVectors: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoSeperateForceVectors * clone() const;
  };
  class HemoUnifyForceVectors: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoUnifyForceVectors * clone() const;
  };
  class HemoSpreadParticleForce: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoSpreadParticleForce * clone() const;
  }; 
  class HemoFindInternalParticleGridPoints: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoFindInternalParticleGridPoints * clone() const;
  };
  
  class HemoInternalGridPointsMembrane: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoInternalGridPointsMembrane * clone() const;
  };
  
  
  class HemoInterpolateFluidVelocity: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoInterpolateFluidVelocity * clone() const;
  };
  class HemoAdvanceParticles: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoAdvanceParticles * clone() const;
  };
  class HemoApplyConstitutiveModel: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoApplyConstitutiveModel * clone() const;
  public:
   bool forced = false;
  };
  class HemoRepulsionForce: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoRepulsionForce * clone() const;
  };
  class HemoBoundaryRepulsionForce: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoBoundaryRepulsionForce * clone() const;
  };
  class HemoDeleteIncompleteCells: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoDeleteIncompleteCells * clone() const;
  public:
    bool verbose;
  };
  class HemoSyncEnvelopes: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoSyncEnvelopes * clone() const;
   void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;
  };
  class HemoGetParticles: public HemoCellFunctional {
    vector<HemoCellParticle *> & particles;
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    HemoGetParticles * clone() const;
  public:
    HemoGetParticles(vector<HemoCellParticle *> & particles_) : particles(particles_) {}
  };
  class HemoSetParticles: public HemoCellFunctional {
    vector<HemoCellParticle> & particles;
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    HemoSetParticles * clone() const;
  public:
    HemoSetParticles(vector<HemoCellParticle> & particles_) : particles(particles_) {}
  };
  class HemoPopulateBoundaryParticles: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoPopulateBoundaryParticles * clone() const;
  };
  class HemoDeleteNonLocalParticles: public HemoCellFunctional {
   void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
   HemoDeleteNonLocalParticles * clone() const;
   public:
    int envelopeSize;
    HemoDeleteNonLocalParticles(int envelope_) : envelopeSize(envelope_) {}
  };
  class HemoSolidifyCells: public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    HemoSolidifyCells * clone() const;
  };
};
}
#endif
