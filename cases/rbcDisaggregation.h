#ifndef RBC_DISAGGREGATION_H
#define RBC_DISAGGREGATION_H

#include "cellStretchingForces3D.h"
#include "immersedCellParticleVtk3D.h"


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
class RBCDisaggregation3D {
public:
        RBCDisaggregation3D(
                MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles_,
                T scalarForce_,
                plint numParticlesPerSide_, plint flowType_,
                T dx_, T dt_, T dNewton_,
                std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_);
        virtual ~RBCDisaggregation3D();
        void applyForce();
        void fixPositions() ;
        void write(plint iter=0);
public:
        void setScalarForce(T scalarForce_);
        T getScalarForce();
private:
        std::vector<plint> outerLeftTags, outerRightTags;
        std::vector<std::vector<plint>*> lateralCellParticleTags;

        MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles;
        T scalarForce;
        plint numParticlesPerSide;
        plint flowType;
        T dx, dt, dNewton;
        std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D;
        bool checkpointed;
};


template<typename T, template<typename U> class Descriptor>
class ZeroOutForceVelocity3D : public BoxProcessingFunctional3D
{
public:
    ZeroOutForceVelocity3D (std::vector<plint> const& pTags_,
                            std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_,
                            plint zeroOut_=0); // 1: Velocity, 2: Force, -- 0: Both
    virtual ~ZeroOutForceVelocity3D() {} ;
    ZeroOutForceVelocity3D(ZeroOutForceVelocity3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ZeroOutForceVelocity3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> const& pTags;
    std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D;
    plint zeroOut;
};



#include "rbcDisaggregation.hh"

#endif  // RBC_DISAGGREGATION_HH

