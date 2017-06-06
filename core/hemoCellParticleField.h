#ifndef HEMOCELLParticleField3D
#define HEMOCELLParticleField3D

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
    void syncEnvelopes() {};

    void setlocalDomain(Box3D & localDomain_);

    void computeGridPosition ( Array<double,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const;
    inline void computeGridPosition ( Array<double,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const;
    
    bool isContainedABS(Array<double,3> pos, Box3D domain) const;

    //Ugly output functions:
    void outputPositions(Box3D,vector<vector<double>>&, pluint, std::string&); 
    void outputForces   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputForceVolume   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputForceArea   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputForceBending   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputForceLink   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputForceVisc    (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputTriangles   (Box3D,vector<vector<plint>>&, vector<vector<double>>&, pluint, std::string&);
    void AddOutputMap();
    map<int,void (HemoCellParticleField::*)(Box3D,vector<vector<double>>&,pluint,std::string&)> outputFunctionMap;
    void passthroughpass(int,Box3D,vector<vector<double>>&,pluint,std::string&);

public:
    virtual HemoCellParticleDataTransfer& getDataTransfer();
    virtual HemoCellParticleDataTransfer const& getDataTransfer() const;
    static std::string getBlockName();
    static HemoCellFields* cellFields;
    pluint atomicBlockId;
    BlockLattice3D<double, DESCRIPTOR> * atomicLattice;
    vector<plint> neighbours;
    pluint envelopeSize;
    pluint getsize() { return particles.size();}
    plint nearestCell(double const) const;
    static std::string basicType() {return std::string(NativeType<double>::getName());}
    static std::string descriptorType() {
      return std::string(DESCRIPTOR<double>::name);
    }
    vector<HemoCellParticle> particles;
    private:

    vector<vector<unsigned int>> particles_per_type;
    public:
    map<int,vector<int>> particles_per_cell;
    map<int,bool> lpc;
    
    //vector<vector<vector<vector<HemoCellParticle*>>>> particle_grid; //maybe better to make custom data structure, But that would be slower
    void insert_ppc(HemoCellParticle* particle,unsigned int index);
    HemoCellParticleDataTransfer dataTransfer;
public:
    Box3D localDomain;

};

#endif
