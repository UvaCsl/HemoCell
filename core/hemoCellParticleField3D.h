#ifndef HEMOCELLParticleField3D
#define HEMOCELLParticleField3D

#include "core/globalDefs.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "particles/particle3D.h"
#include "cellFields3D.h"

class HemoParticleField3D;

class HemoParticleDataTransfer3D : public BlockDataTransfer3D {
public:
    HemoParticleDataTransfer3D(HemoParticleField3D& particleField_);
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset);
    virtual void receive( Box3D domain, std::vector<char> const& buffer,
                          modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset);
private:
    HemoParticleField3D& particleField;
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
    virtual void addParticle(Box3D domain, Particle3D<double,DESCRIPTOR>* particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles(Box3D domain,plint tag);
    virtual void findParticles(Box3D domain,
                               std::vector<Particle3D<double,DESCRIPTOR>*>& found);
    void findParticles(Box3D domain,
                               std::vector<Particle3D<double,DESCRIPTOR> *>& found,
                               pluint type) const;
    virtual void advanceParticles(Box3D domain);

    int deleteIncompleteCells(pluint ctype);

    void setlocalDomain(Box3D & localDomain_);

    void computeGridPosition ( Array<double,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const;
    
    
    bool isContainedABS(Array<double,3> pos, Box3D domain) const;

    //Ugly output functions:
    void outputPositions(Box3D,vector<vector<double>>&, pluint, std::string&); 
    void outputForces   (Box3D,vector<vector<double>>&, pluint, std::string&);
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
private:
    std::vector<Particle3D<double,DESCRIPTOR>*> particles;
    std::vector<std::vector<Particle3D<double,DESCRIPTOR>*>> particles_per_type;
    std::map<int,std::vector<SurfaceParticle3D*>> particles_per_cell;
    std::map<int,bool> lpc;
    void insert_ppc(SurfaceParticle3D* particle);
    HemoParticleDataTransfer3D dataTransfer;
    Box3D localDomain;

};
CellFields3D* HemoParticleField3D::cellFields=0;

#endif
