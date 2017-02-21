#ifndef RBC_DISAGGREGATION_HH
#define RBC_DISAGGREGATION_HH

#include "rbcDisaggregation.h"



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
        RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::RBCDisaggregation3D(
                MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles_,
                T scalarForce_,
                plint numParticlesPerSide_, plint flowType_,
                T dx_, T dt_, T dNewton_,
                std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_) :
        particles(particles_), scalarForce(scalarForce_), numParticlesPerSide(numParticlesPerSide_), flowType(flowType_),
        dx(dx_), dt(dt_), dNewton(dNewton_), tagToParticle3D(tagToParticle3D_)
{
    //PLB_PRECONDITION( npar == 1 && MPI::COMM_WORLD.Get_size() == 1 );
    if (flowType == 2) {
        lateralCellParticleTags.push_back(&outerLeftTags);
        lateralCellParticleTags.push_back(&outerRightTags);

        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(&particles);
        applyProcessingFunctional (
            new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerLeftTags, &outerRightTags, TFL_DISAGGREGATION_UP),
            particles.getBoundingBox(), particleArg );
        applyForce();
        fixPositions();
    }
}



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::~RBCDisaggregation3D() {};



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::setScalarForce(T scalarForce_) {
        scalarForce = scalarForce_;
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
T RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::getScalarForce() {
        return scalarForce;
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::applyForce()
{
    if (flowType == 2) {
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(&particles);

        std::vector<std::vector<plint> > particleTags;
        particleTags.push_back(outerRightTags);
        std::vector<Array<T,3> > forces;
        forces.push_back(Array<T,3>(0, scalarForce, 0));
        T cellDensity = 1.0;
        applyProcessingFunctional (
                new ApplyStretchingForce3D<T,Descriptor>(particleTags, forces, cellDensity, tagToParticle3D),
                particles.getBoundingBox(), particleArg );

        applyProcessingFunctional ( // Zero Out force (==2)
                new ZeroOutForceVelocity3D<T,Descriptor>(outerLeftTags, tagToParticle3D, 2),
                particles.getBoundingBox(), particleArg );
    }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::fixPositions()
{
    if (flowType == 2) {
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(&particles);

        applyProcessingFunctional ( // Zero Out velocity (==1)
                new ZeroOutForceVelocity3D<T,Descriptor>(outerLeftTags, tagToParticle3D, 1),
                particles.getBoundingBox(), particleArg );
    }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void RBCDisaggregation3D<T,Descriptor,ParticleFieldT>::write(plint iter) { }


/* ================================================================================ */
/* ******** ZeroOutForceVelocity3D *********************************** */
/* ================================================================================ */

template<typename T, template<typename U> class Descriptor>
ZeroOutForceVelocity3D<T,Descriptor>::ZeroOutForceVelocity3D (std::vector<plint> const& pTags_,
        std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_, plint zeroOut_)
    :  pTags(pTags_), tagToParticle3D(tagToParticle3D_), zeroOut(zeroOut_)
{ }


template<typename T, template<typename U> class Descriptor>
ZeroOutForceVelocity3D<T,Descriptor>::ZeroOutForceVelocity3D (
        ZeroOutForceVelocity3D<T,Descriptor> const& rhs)
    : pTags(rhs.pTags),
      tagToParticle3D(rhs.tagToParticle3D),
      zeroOut(rhs.zeroOut)
{ }


template<typename T, template<typename U> class Descriptor>
void ZeroOutForceVelocity3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );

    pluint nParts = pTags.size();
    for (pluint iT = 0; iT < nParts; ++iT) {
        plint tag = pTags[iT];
        if (tagToParticle3D.count(tag) > 0) {
            SurfaceParticle3D * particle = dynamic_cast<SurfaceParticle3D*>( tagToParticle3D[tag]);
//                particle->get_a()     += forces[var] * (1.0/nParts)/cellDensity;
            if (zeroOut == 0 or zeroOut == 1) { particle->get_v() = particle->get_vPrevious() = Array<T,3>(0,0,0); }
            if (zeroOut == 0 or zeroOut == 2) { particle->get_force() = Array<T,3>(0,0,0); }
        } else pcout << "ImmerseCellParticle3D not found! Something is wrong here!" << std::endl;
    }
}

template<typename T, template<typename U> class Descriptor>
ZeroOutForceVelocity3D<T,Descriptor>*
    ZeroOutForceVelocity3D<T,Descriptor>::clone() const
{
    return new ZeroOutForceVelocity3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ZeroOutForceVelocity3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ZeroOutForceVelocity3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ZeroOutForceVelocity3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field.
}



#endif  // RBC_DISAGGREGATION_HH

