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

#include "hemocell_internal.h"
#include "hemoCellFields.h"
#include "hemoCellParticleDataTransfer.h"
#include "hemoCellParticle.h"

class HemoCellParticleField : public AtomicBlock3D {
public:
    HemoCellParticleField(plint nx, plint ny, plint nz);
    virtual ~HemoCellParticleField();
    HemoCellParticleField(HemoCellParticleField const& rhs);
    HemoCellParticleField& operator=(HemoCellParticleField const& rhs);
    HemoCellParticleField* clone() const;
    void swap(HemoCellParticleField& rhs);
    virtual void applyConstitutiveModel(bool forced = false);
    virtual void addParticle(Box3D domain, HemoCellParticle* particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles_inverse(Box3D domain);
    virtual void removeParticles(Box3D domain,plint tag);
    virtual void removeParticles(plint tag);
    virtual void findParticles(Box3D domain,
                               std::vector<HemoCellParticle*>& found);
    void findParticles(Box3D domain,
                               std::vector<HemoCellParticle*>& found,
                               pluint type);
    virtual void advanceParticles();
    void applyRepulsionForce(bool forced = false);
    virtual void interpolateFluidVelocity(Box3D domain);
    virtual void spreadParticleForce(Box3D domain);
    void separateForceVectors();
    void unifyForceVectors();

    int deleteIncompleteCells(pluint ctype, bool verbose=true);
    int deleteIncompleteCells(bool verbose=true);
    void syncEnvelopes();

    void setlocalDomain(Box3D & localDomain_);

    void computeGridPosition ( hemo::Array<T,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const;
    inline void computeGridPosition ( hemo::Array<T,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const;
    
    inline bool isContainedABS(hemo::Array<T,3> pos, Box3D box) const {
	    Dot3D const& location = this->getLocation();
	    T x = pos[0]-location.x;
	    T y = pos[1]-location.y;
	    T z = pos[2]-location.z;
	    //if (box.z1 < -10000) {
	    //  exit(0);
	    //} 
	    //if (location.x < -10000) {
	    //  exit(0);
	    //} 

	    return (x > box.x0-0.5) && (x <= box.x1+0.5) &&
	           (y > box.y0-0.5) && (y <= box.y1+0.5) &&
	           (z > box.z0-0.5) && (z <= box.z1+0.5);

	}

    //Ugly output functions:
    void outputPositions(Box3D,vector<vector<T>>&, pluint, std::string&); 
    void outputForces   (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceVolume   (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceArea   (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceBending   (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceLink   (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceVisc    (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputTriangles   (Box3D,vector<vector<plint>>&, vector<vector<T>>&, pluint, std::string&);
    void outputLines   (Box3D,vector<vector<plint>>&, vector<vector<T>>&, plint, std::string&);
    void outputVertexId    (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputCellId    (Box3D,vector<vector<T>>&, pluint, std::string&);
    void outputForceInnerLink   (Box3D,vector<vector<T>>&, pluint, std::string&);

    void AddOutputMap();
    map<int,void (HemoCellParticleField::*)(Box3D,vector<vector<T>>&,pluint,std::string&)> outputFunctionMap;
    void passthroughpass(int,Box3D,vector<vector<T>>&,pluint,std::string&);

public:
    virtual HemoCellParticleDataTransfer& getDataTransfer();
    virtual HemoCellParticleDataTransfer const& getDataTransfer() const;
    static std::string getBlockName();
    static HemoCellFields* cellFields;
    pluint atomicBlockId;
    BlockLattice3D<T, DESCRIPTOR> * atomicLattice;
    vector<plint> neighbours;
    pluint envelopeSize;
    pluint getsize() { return particles.size();}
    plint nearestCell(T const) const;
    static std::string basicType() {return std::string(NativeType<T>::getName());}
    static std::string descriptorType() {
      return std::string(DESCRIPTOR<T>::name);
    }
    vector<HemoCellParticle> particles;
    int nFluidCells;
    
private:
  bool lpc_up_to_date = false;
  bool ppt_up_to_date = false;
  bool ppc_up_to_date = false;
  bool preinlet_ppc_up_to_date = false;
  vector<vector<unsigned int>> _particles_per_type;
  map<int,vector<int>> _particles_per_cell;
  map<int,vector<int>> _preinlet_particles_per_cell;
  map<int,bool> _lpc;
  void update_lpc();
  void update_ppc();
  void update_preinlet_ppc();
  void update_ppt();
  void issueWarning(HemoCellParticle & p);
  
public:
  const vector<vector<unsigned int>> & get_particles_per_type(); 
  const map<int,vector<int>> & get_particles_per_cell();
  const map<int,vector<int>> & get_preinlet_particles_per_cell();
  const map<int,bool> & get_lpc();
    
    //vector<vector<vector<vector<HemoCellParticle*>>>> particle_grid; //maybe better to make custom data structure, But that would be slower
    void insert_ppc(HemoCellParticle* particle,unsigned int index);
    void insert_preinlet_ppc(HemoCellParticle* particle,unsigned int index);

    HemoCellParticleDataTransfer dataTransfer;
public:
    Box3D localDomain;

};

#endif
