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

#include "atomicBlock/blockLattice3D.hh"

namespace hemo {
using namespace std;
class HemoCellParticleField : public plb::AtomicBlock3D {
public:
    HemoCellParticleField(plint nx, plint ny, plint nz);
    virtual ~HemoCellParticleField();
    HemoCellParticleField(HemoCellParticleField const& rhs);
    HemoCellParticleField& operator=(HemoCellParticleField const& rhs);
    HemoCellParticleField* clone() const;
    void swap(HemoCellParticleField& rhs);
    virtual void applyConstitutiveModel(bool forced = false);
    virtual void addParticle(HemoCellParticle* particle);
    void addParticle(const HemoCellParticle::serializeValues_t & sv);
    void addParticlePreinlet(const HemoCellParticle::serializeValues_t & sv);

    virtual void removeParticles(plb::Box3D domain);
    virtual void removeParticles_inverse(plb::Box3D domain);
    virtual void removeParticles(plb::Box3D domain,plint tag);
    virtual void removeParticles(plint tag);
    virtual void findParticles(plb::Box3D domain,
                               std::vector<HemoCellParticle*>& found);
    virtual void findParticles(plb::Box3D domain,
                               std::vector<const HemoCellParticle*>& found) const;
    void findParticles(plb::Box3D domain,
                               std::vector<HemoCellParticle*>& found,
                               pluint type);
    virtual void advanceParticles();
    void applyRepulsionForce(bool forced = false);
    virtual void interpolateFluidVelocity(plb::Box3D domain);
    virtual void spreadParticleForce(plb::Box3D domain);
    void separateForceVectors();
    void unifyForceVectors();
    
    virtual void findInternalParticleGridPoints(plb::Box3D domain);
    virtual void internalGridPointsMembrane(plb::Box3D domain);
    
    int deleteIncompleteCells(pluint ctype, bool verbose=true);
    int deleteIncompleteCells(bool verbose=true);
    void syncEnvelopes();
    void populateBoundaryParticles();
    void applyBoundaryRepulsionForce();
    void populateBindingSites(plb::Box3D & domain);

    T eigenValueFromCell(plb::Cell<T,DESCRIPTOR> & cell);
    
    void solidifyCells();
    
    void setlocalDomain(plb::Box3D & localDomain_);

    void computeGridPosition ( hemo::Array<T,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const;
    inline void computeGridPosition ( hemo::Array<T,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const;
    
    inline bool isContainedABS(const hemo::Array<T,3> & pos, const plb::Box3D & box) const {
	    plb::Dot3D const& location = this->getLocation();
	    T x = pos[0]-location.x;
	    T y = pos[1]-location.y;
	    T z = pos[2]-location.z;

	    return (x > box.x0-0.5) && (x <= box.x1+0.5) &&
	           (y > box.y0-0.5) && (y <= box.y1+0.5) &&
	           (z > box.z0-0.5) && (z <= box.z1+0.5);

    }
    plb::Box3D & getBoundingBox() {
      return boundingBox;
    }
    //Ugly output functions:
    void outputPositions(plb::Box3D,vector<vector<T>>&, pluint, std::string&); 
    void outputVelocities(plb::Box3D,vector<vector<T>>&, pluint, std::string&); 
    void outputForces   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceVolume   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceArea   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceBending   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceLink   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceVisc    (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceRepulsion  (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    
    void outputTriangles   (plb::Box3D,vector<vector<plint>>&, pluint, std::string&);
    void outputInnerLinks   (plb::Box3D,vector<vector<plint>>&, pluint, std::string&);
    
    void outputVertexId    (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputCellId    (plb::Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceInnerLink   (plb::Box3D,vector<vector<T>>&, pluint, std::string&);


    void AddOutputMap();
    map<int,void (HemoCellParticleField::*)(plb::Box3D,vector<vector<T>>&,pluint,std::string&)> outputFunctionMap;
    void passthroughpass(int,plb::Box3D,vector<vector<T>>&,pluint,std::string&);

public:
    virtual HemoCellParticleDataTransfer& getDataTransfer();
    virtual HemoCellParticleDataTransfer const& getDataTransfer() const;
    static std::string getBlockName();
    static HemoCellFields* cellFields;
    pluint atomicBlockId;
    plb::BlockLattice3D<T, DESCRIPTOR> * atomicLattice = 0;
    plb::BlockLattice3D<T, CEPAC_DESCRIPTOR> * CEPAClattice = 0;

    vector<plint> neighbours;
    vector<plb::Dot3D> boundaryParticles;
    pluint envelopeSize;
    pluint getsize() { return particles.size();}
    plint nearestCell(T const) const;
    static std::string basicType() {return std::string(plb::NativeType<T>::getName());}
    static std::string descriptorType() {
      return std::string(DESCRIPTOR<T>::name);
    }
    vector<HemoCellParticle> particles;
    plb::Box3D boundingBox; 
    int nFluidCells = 0;
    
private:
  bool lpc_up_to_date = false;
  bool ppt_up_to_date = false;
  bool ppc_up_to_date = false;
  bool preinlet_ppc_up_to_date = false;
  bool pg_up_to_date = false;
public:
  void invalidate_lpc() { lpc_up_to_date = false;};
  void invalidate_ppt() { ppt_up_to_date = false;};
  void invalidate_ppc() { ppc_up_to_date = false;};
  void invalidate_preinlet_ppc() { preinlet_ppc_up_to_date = false;};
  void invalidate_pg() { pg_up_to_date = false;};
private:
  vector<vector<unsigned int>> _particles_per_type;
  map<int,vector<int>> _particles_per_cell;
  map<int,vector<int>> _preinlet_particles_per_cell;
  map<int,bool> _lpc;
  void update_lpc();
  void update_ppc();
  void update_preinlet_ppc();
  void update_ppt();
  void update_pg();
  void issueWarning(HemoCellParticle & p);
  
  hemo::Array<unsigned int,10> * particle_grid = 0;
  unsigned int * particle_grid_size = 0;
  unsigned int grid_index(int & nx,int & ny,int & nz) {
    return nz+this->atomicLattice->getNz()*(ny+(this->atomicLattice->getNy()*nx));
  }
  
  vector<hemo::Array<T,3>*> allocated_for_output;
  
public:
  const vector<vector<unsigned int>> & get_particles_per_type(); 
  const map<int,vector<int>> & get_particles_per_cell();
  const map<int,vector<int>> & get_preinlet_particles_per_cell();
  const map<int,bool> & get_lpc();
  
  set<plb::Dot3D> internalPoints; // Store found interior points
  plb::ScalarField3D<T> * interiorViscosityField = 0;
  
    
    //vector<vector<vector<vector<HemoCellParticle*>>>> particle_grid; //maybe better to make custom data structure, But that would be slower
    void insert_ppc(HemoCellParticle* particle,unsigned int index);
    void insert_preinlet_ppc(HemoCellParticle* particle,unsigned int index);

    HemoCellParticleDataTransfer & particleDataTransfer;
public:
    plb::Box3D localDomain;
    
    //These should be edited through the helper/solidifyField.h functions
    std::set<plb::Dot3D> bindingSites;
    plb::ScalarField3D<bool> * bindingField = 0;
};

}
#endif
