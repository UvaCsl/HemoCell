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
#ifndef HEMOCELLPARTICLEFIELD_H
#define HEMOCELLPARTICLEFIELD_H

namespace hemo {
class HemoCellParticleField;
}

#include "hemoCellFields.h"
#include "hemoCellParticleDataTransfer.h"
#include "hemoCellParticle.h"
#include "particleGrid.h"

#include "atomicBlock/blockLattice3D.hh"

namespace hemo {
class HemoCellParticleField : public plb::AtomicBlock3D {
public:
  HemoCellParticleField(plint nx, plint ny, plint nz);
  virtual ~HemoCellParticleField();
  HemoCellParticleField(HemoCellParticleField const & rhs);
  /// Copy the particles of @p rhs into this block (used by Palabos when a
  /// MultiParticleField3D is copied). The geometry, data transfer and
  /// lattice links of this block are kept. (The former implementation
  /// deleted `this` and returned a new object.)
  HemoCellParticleField & operator=(HemoCellParticleField const & rhs);
  HemoCellParticleField * clone() const;
  virtual void applyConstitutiveModel(bool forced = false);
  virtual void addParticle(HemoCellParticle * particle);
  void addParticle(const HemoCellParticle::serializeValues_t & sv);
  void addParticlePreinlet(const HemoCellParticle::serializeValues_t & sv);

  virtual void removeParticles(plb::Box3D domain);
  virtual void removeParticles_inverse(plb::Box3D domain);
  virtual void removeParticles(plb::Box3D domain, plint tag);
  virtual void removeParticles(plint tag);
  virtual void findParticles(plb::Box3D domain, std::vector<HemoCellParticle *> & found);
  virtual void findParticles(plb::Box3D domain, std::vector<const HemoCellParticle *> & found) const;
  void findParticles(plb::Box3D domain, std::vector<HemoCellParticle *> & found, pluint type);
  virtual void advanceParticles();
  void applyRepulsionForce(bool forced = false);
  /// Scale sv.force_repulsion down to magnitude @p limit where it is larger
  /// [lattice force units].
  void capRepulsionForce(T limit);
  virtual void interpolateFluidVelocity(plb::Box3D domain);
  virtual void spreadParticleForce(plb::Box3D domain);
  void separateForceVectors();
  void unifyForceVectors();
  void updateResidenceTime(unsigned int rtime);

  virtual void findInternalParticleGridPoints(plb::Box3D domain);
  virtual void internalGridPointsMembrane(plb::Box3D domain);

  /// Delete the particles of cells of type @p ctype that are not complete in this block.
  int deleteIncompleteCells(pluint ctype, bool verbose = true);
  /// Delete the particles of all cells that are not complete in this block.
  int deleteIncompleteCells(bool verbose = true);
  /**
   * @brief Delete the particles of cells that lost vertices.
   *
   * Deletes the cells that are incomplete in this block and have a vertex
   * within @p envelope nodes of the bulk. Right after syncEnvelopes() and
   * advanceParticles() those are the cells with a vertex missing everywhere
   * (it entered a wall, left the domain or had a non-finite force), provided
   * the particle envelope holds whole cells: cells cut by the edge of the
   * particle envelope lie further out and are left alone.
   * @param[in] envelope Distance from the bulk [lattice nodes].
   * @param[in] verbose Print a warning per deleted cell from the block that owns its vertex.
   * @return Number of particles deleted.
   */
  int deleteBrokenCells(plint envelope, bool verbose);
  /// Replace the particles by a copy, keeping the originals aside, until
  /// restoreParticles() (see HemoCellFields::CompleteCellsScope).
  void snapshotParticles();
  /// Discard the current particles and put back those kept by snapshotParticles().
  void restoreParticles();
  void populateBoundaryParticles();
  void applyBoundaryRepulsionForce(bool resetFirst = false);
  void populateBindingSites(plb::Box3D & domain);

  T eigenValueFromCell(plb::Cell<T, DESCRIPTOR> & cell);

  void solidifyCells();
  void prepareSolidification();

  void setlocalDomain(plb::Box3D & localDomain_);

  void computeGridPosition(hemo::Array<T, 3> const & position, plint & iX, plint & iY, plint & iZ) const;
  inline void computeGridPosition(hemo::Array<T, 3> const & position, plint * iX, plint * iY, plint * iZ) const;

  inline bool isContainedABS(const hemo::Array<T, 3> & pos, const plb::Box3D & box) const {
    plb::Dot3D const & location = this->getLocation();
    T x = pos[0] - location.x;
    T y = pos[1] - location.y;
    T z = pos[2] - location.z;

    return (x > box.x0 - 0.5) && (x <= box.x1 + 0.5) && (y > box.y0 - 0.5) && (y <= box.y1 + 0.5) &&
           (z > box.z0 - 0.5) && (z <= box.z1 + 0.5);
  }
  plb::Box3D & getBoundingBox() { return boundingBox; }
  // Keep AtomicBlock3D's const version visible (same box), e.g. for the
  // PLB_ASSERTs in const member functions.
  using plb::AtomicBlock3D::getBoundingBox;
  /**
   * @brief The cells of type @p ctype that this block writes to the cell output.
   *
   * These are the complete cells whose vertex 0 lies in this block's bulk, in
   * ascending cell ID, so every cell is written by exactly one block. Envelope
   * copies, including the shifted images of cells across a periodic boundary,
   * are not written. All per-vertex outputs and the triangle and inner-link
   * lists of a block use this list, so their rows correspond.
   *
   * @param[in] ctype Cell type index.
   * @param[out] cellIds Cell IDs to write.
   */
  void getOutputCells(pluint ctype, std::vector<plint> & cellIds);
  //Ugly output functions:
  void outputPositions(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputVelocities(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForces(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceVolume(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceArea(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceBending(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceLink(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceVisc(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputForceRepulsion(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);

  void outputTriangles(plb::Box3D, std::vector<std::vector<plint>> &, pluint, std::string &);
  void outputInnerLinks(plb::Box3D, std::vector<std::vector<plint>> &, pluint, std::string &);

  void outputVertexId(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  /// One row per vertex of the cells this block writes (getOutputCells()),
  /// holding the cell's ID (exact, not converted to T).
  void outputCellIds(std::vector<std::vector<plint>> & output, pluint ctype, std::string & name);
  void outputForceInnerLink(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);
  void outputResTime(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);


  void AddOutputMap();
  std::map<int, void (HemoCellParticleField::*)(plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &)>
      outputFunctionMap;
  void passthroughpass(int, plb::Box3D, std::vector<std::vector<T>> &, pluint, std::string &);

public:
  virtual HemoCellParticleDataTransfer & getDataTransfer();
  virtual HemoCellParticleDataTransfer const & getDataTransfer() const;
  static std::string getBlockName();
  static HemoCellFields * cellFields;
  pluint atomicBlockId;
  plb::BlockLattice3D<T, DESCRIPTOR> * atomicLattice = 0;
  plb::BlockLattice3D<T, CEPAC_DESCRIPTOR> * CEPAClattice = 0;

  std::vector<plint> neighbours;
  std::vector<plb::Dot3D> boundaryParticles;
  pluint envelopeSize;
  pluint getsize() { return particles.size(); }
  plint nearestCell(T const) const;
  static std::string basicType() { return std::string(plb::NativeType<T>::getName()); }
  static std::string descriptorType() { return std::string(DESCRIPTOR<T>::name); }
  std::vector<HemoCellParticle> particles;
  plb::Box3D boundingBox;
  int nFluidCells = 0;

private:
  bool lpc_up_to_date = false;
  bool ppt_up_to_date = false;
  bool ppc_up_to_date = false;
  bool preinlet_ppc_up_to_date = false;
  bool pg_up_to_date = false;

public:
  void invalidate_lpc() { lpc_up_to_date = false; };
  void invalidate_ppt() { ppt_up_to_date = false; };
  void invalidate_ppc() { ppc_up_to_date = false; };
  void invalidate_preinlet_ppc() { preinlet_ppc_up_to_date = false; };
  void invalidate_pg() { pg_up_to_date = false; };

private:
  std::vector<std::vector<unsigned int>> _particles_per_type;
  std::map<plint, std::vector<int>> _particles_per_cell;
  std::map<plint, std::vector<int>> _preinlet_particles_per_cell;
  std::map<plint, bool> _lpc;
  void update_lpc();
  void update_ppc();
  void update_preinlet_ppc();
  void update_ppt();
  void update_pg();
  void issueWarning(HemoCellParticle & p);
  /// Delete incomplete cells: only of type @p ctype if @p onlyType, only those
  /// with a vertex in @p near if it is not null.
  int deleteIncompleteCellsImpl(bool onlyType, pluint ctype, const plb::Box3D * near, bool verbose);
  /// The particles kept aside by snapshotParticles().
  std::vector<HemoCellParticle> particlesSnapshot;
  /// Zero the force of a vertex with a non-finite force and tag it (tag 1)
  /// for removal by spreadParticleForce(); warns for the first one (count 0).
  void dropNonFiniteParticle(HemoCellParticle & particle, plint & count);

  /// Particles binned per lattice node; valid only while pg_up_to_date is
  /// true (any change to `particles` invalidates it, see update_pg()).
  ParticleGrid particleGrid;

  std::vector<hemo::Array<T, 3> *> allocated_for_output;

public:
  const std::vector<std::vector<unsigned int>> & get_particles_per_type();
  const std::map<plint, std::vector<int>> & get_particles_per_cell();
  const std::map<plint, std::vector<int>> & get_preinlet_particles_per_cell();
  const std::map<plint, bool> & get_lpc();

  std::set<plb::Dot3D> internalPoints; // Store found interior points
  plb::ScalarField3D<T> * interiorViscosityField = 0;


  //vector<vector<vector<vector<HemoCellParticle*>>>> particle_grid; //maybe better to make custom data structure, But that would be slower
  void insert_ppc(HemoCellParticle * particle, unsigned int index);
  void insert_preinlet_ppc(HemoCellParticle * particle, unsigned int index);

  HemoCellParticleDataTransfer & particleDataTransfer;

public:
  plint addParticleCount = 0;
  plb::Box3D localDomain;

  //These should be edited through the helper/solidifyField.h functions
  std::set<plb::Dot3D> bindingSites;
  plb::ScalarField3D<bool> * bindingField = 0;
};

}
#endif
