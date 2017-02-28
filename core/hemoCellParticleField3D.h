#ifndef HEMOCELLParticleField3D
#define HEMOCELLParticleField3D

#include "core/globalDefs.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "particles/particle3D.h"
#include "cellFields3D.h"

template<typename T, template<typename U> class Descriptor> class HemoParticleField3D;

template<typename T, template<typename U> class Descriptor>
class HemoParticleDataTransfer3D : public BlockDataTransfer3D {
public:
    HemoParticleDataTransfer3D(HemoParticleField3D<T,Descriptor>& particleField_);
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
    HemoParticleField3D<T,Descriptor>& particleField;
};

template<typename T, template<typename U> class Descriptor>
class HemoParticleField3D : public ParticleField3D<T,Descriptor> {
public:
    typedef Particle3D<T,Descriptor> ParticleT;
public:
    HemoParticleField3D(plint nx, plint ny, plint nz);
    virtual ~HemoParticleField3D();
    HemoParticleField3D(HemoParticleField3D<T,Descriptor> const& rhs);
    HemoParticleField3D<T,Descriptor>& operator=(HemoParticleField3D<T,Descriptor> const& rhs);
    HemoParticleField3D<T,Descriptor>* clone() const;
    void swap(HemoParticleField3D<T,Descriptor>& rhs);
public:
    virtual void addParticle(Box3D domain, Particle3D<T,Descriptor>* particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles(Box3D domain,plint tag);
    virtual void findParticles(Box3D domain,
                               std::vector<Particle3D<T,Descriptor>*>& found);
    virtual void findParticles(Box3D domain,
                               std::vector<Particle3D<T,Descriptor> const*>& found) const;
    void findParticles(Box3D domain,
                               std::vector<Particle3D<T,Descriptor> *>& found,
                               pluint type) const;
    virtual void velocityToParticleCoupling(Box3D domain, TensorField3D<T,3>& velocity, T scaling=0.);
    virtual void velocityToParticleCoupling(Box3D domain, NTensorField3D<T>& velocity, T scaling=0.);
    virtual void rhoBarJtoParticleCoupling(Box3D domain, NTensorField3D<T>& rhoBarJ, bool velIsJ, T scaling=0.);
    virtual void fluidToParticleCoupling(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, T scaling=0.);
    virtual void advanceParticles(Box3D domain, T cutOffValue=-1.);

    int deleteIncompleteCells();

    void setlocalDomain(Box3D & localDomain_) {localDomain = localDomain_;}

    //Ugly output functions:
    void outputPositions(Box3D,vector<vector<double>>&, int, std::string&); 
    void outputForces   (Box3D,vector<vector<double>>&, int, std::string&);
    void AddOutputMap();
    map<int,void (HemoParticleField3D<T,Descriptor>::*)(Box3D,vector<vector<double>>&,int,std::string&)> outputFunctionMap;
    void passthroughpass(int,Box3D,vector<vector<double>>&,int,std::string&);

public:
    virtual HemoParticleDataTransfer3D<T,Descriptor>& getDataTransfer();
    virtual HemoParticleDataTransfer3D<T,Descriptor> const& getDataTransfer() const;
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();
    CellFields3D* cellFields;
private:
    std::vector<Particle3D<T,Descriptor>*> particles;
    std::vector<std::vector<Particle3D<T,Descriptor>*>> particles_per_type;
    std::map<int,std::vector<SurfaceParticle3D*>> particles_per_cell;
    std::map<int,bool> lpc;
    void insert_ppc(SurfaceParticle3D* particle);
    HemoParticleDataTransfer3D<T,Descriptor> dataTransfer;
    Box3D & localDomain;

};

#endif
