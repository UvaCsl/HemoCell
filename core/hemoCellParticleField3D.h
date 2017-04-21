#ifndef HEMOCELLParticleField3D
#define HEMOCELLParticleField3D

#include "hemocell_internal.h"
#include "surfaceParticle3D.h"
#include "cellFields3D.h"

class HemoParticleDataTransfer3D : public BlockDataTransfer3D {
public:
    HemoParticleDataTransfer3D();
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset);
    virtual void receive( Box3D domain, std::vector<char> const& buffer,
                          modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    virtual void setBlock(AtomicBlock3D& block);
    virtual void setConstBlock(AtomicBlock3D const& block) ;
    virtual HemoParticleDataTransfer3D* clone() const {
      return new HemoParticleDataTransfer3D(*this);
    }
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset);
private:
    HemoParticleField3D* particleField;
    HemoParticleField3D const * constParticleField;
};

class HemoParticleField3D : public AtomicBlock3D {
public:
public:
    HemoParticleField3D(plint nx, plint ny, plint nz);
    virtual ~HemoParticleField3D();
    HemoParticleField3D(HemoParticleField3D const& rhs);
    HemoParticleField3D& operator=(HemoParticleField3D const& rhs);
    HemoParticleField3D* clone() const;
    void swap(HemoParticleField3D& rhs);
public:
    virtual void applyConstitutiveModel();
    virtual void addParticle(Box3D domain, SurfaceParticle3D* particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles(Box3D domain,plint tag);
    virtual void removeParticles(plint tag);
    virtual void findParticles(Box3D domain,
                               std::vector<SurfaceParticle3D*>& found);
    void findParticles(Box3D domain,
                               std::vector<SurfaceParticle3D*>& found,
                               pluint type) const;
    virtual void advanceParticles();
    void applyRepulsionForce();
    virtual void interpolateFluidVelocity(Box3D domain);
    virtual void spreadParticleForce(Box3D domain);
    void separateForceVectors();
    void unifyForceVectors();

    int deleteIncompleteCells(pluint ctype, bool twice=false);
    int deleteIncompleteCells(bool twice=false);
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
    void outputForceInPlane   (Box3D,vector<vector<double>>&, pluint, std::string&);
    void outputTriangles   (Box3D,vector<vector<plint>>&, vector<vector<double>>&, pluint, std::string&);
    void AddOutputMap();
    map<int,void (HemoParticleField3D::*)(Box3D,vector<vector<double>>&,pluint,std::string&)> outputFunctionMap;
    void passthroughpass(int,Box3D,vector<vector<double>>&,pluint,std::string&);

public:
    virtual HemoParticleDataTransfer3D& getDataTransfer();
    virtual HemoParticleDataTransfer3D const& getDataTransfer() const;
    static std::string getBlockName();
    static CellFields3D* cellFields;
    pluint atomicBlockId;
    BlockLattice3D<double, DESCRIPTOR> * atomicLattice;
    pluint envelopeSize;
    pluint getsize() { return particles.size();}
    plint nearestCell(double const) const;
    static std::string basicType() {return std::string(NativeType<double>::getName());}
    static std::string descriptorType() {
      return std::string(DESCRIPTOR<double>::name);
    }
    std::vector<SurfaceParticle3D*> particles;
    private:

    std::vector<std::vector<SurfaceParticle3D*>> particles_per_type;
    std::map<int,std::vector<SurfaceParticle3D*>> particles_per_cell;
    std::map<int,bool> lpc;
    
    vector<vector<vector<vector<SurfaceParticle3D*>>>> particle_grid; //maybe better to make custom data structure, But that would be slower
    void insert_ppc(SurfaceParticle3D* particle);
    HemoParticleDataTransfer3D dataTransfer;
    Box3D localDomain;

};

#endif
