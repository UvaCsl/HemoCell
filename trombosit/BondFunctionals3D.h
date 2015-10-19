#ifndef BOND_FUNCTIONALS_3D_H
#define BOND_FUNCTIONALS_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include "cellCellForces3D.hh"
#include <string>

using namespace plb;


namespace trombocit {

template<typename T, template<typename U> class Descriptor>
class BondField3D ;

template<typename T, template<typename U> class Descriptor>
class UpdateBondParticles3D : public BoxProcessingFunctional3D
{
public:
    UpdateBondParticles3D (BondField3D<T, Descriptor> & bondField_) : bondField(bondField_) { } ;
    ~UpdateBondParticles3D() { };
    UpdateBondParticles3D(UpdateBondParticles3D<T,Descriptor> const& rhs) : bondField(rhs.bondField) { } ;
    /// Arguments: [0] Particle-field 0
    /// Arguments: [1] Particle-field 1
    /// Arguments: [2] BondParticle-field
    virtual UpdateBondParticles3D<T,Descriptor>* clone() const { return new UpdateBondParticles3D<T,Descriptor>(*this); };
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; } ;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::allVariables; // ParticleField 0
        modified[1] = modif::allVariables; // ParticleField 1
        modified[2] = modif::allVariables; // BondParticles
    } ;
    void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
        isWritten[1] = true;
        isWritten[2] = true;
    };
private:
    BondField3D<T, Descriptor> & bondField;
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
        PLB_PRECONDITION( blocks.size()==3 );
        ParticleField3D<T,Descriptor>& particleField0 =
            *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
        ParticleField3D<T,Descriptor>& particleField1 =
            *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);
        ParticleField3D<T,Descriptor>& bondParticleField =
            *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[2]);

        std::vector<Particle3D<T,Descriptor>*> particles, bondParticles;
        std::map<std::pair<plint,plint>, Particle3D<T,Descriptor>*> particle0Map, particle1Map;
        bondParticleField.findParticles(bondParticleField.getBoundingBox(), bondParticles);
        // Don't perform unnecessary calculations
        bool particleFieldsAreTheSame= (&particleField0) == (&particleField1);



        // Find particles from field 0 and make particle map from field 0;
        particleField0.findParticles(particleField0.getBoundingBox(), particles);
        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
            plint cellId = castParticleToICP3D(particles[iParticle])->get_cellId();
            plint vertexId = castParticleToICP3D(particles[iParticle])->getVertexId();
            particle0Map[std::make_pair<plint,plint>(cellId, vertexId)] = particles[iParticle];
        }

        // Do the same for field 1, but if it's the same just copy it.
        if (particleFieldsAreTheSame) { particle1Map = particle0Map; }
        else {
            particleField1.findParticles(particleField1.getBoundingBox(), particles);
            for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                plint cellId = castParticleToICP3D(particles[iParticle])->get_cellId();
                plint vertexId = castParticleToICP3D(particles[iParticle])->getVertexId();
                particle1Map[std::make_pair<plint,plint>(cellId, vertexId)] = particles[iParticle];
            }
        }

        // Scan all bondParticles and update their edges.
        bondParticleField.findParticles(bondParticleField.getBoundingBox(), bondParticles);
        for (pluint iParticle=0; iParticle<bondParticles.size(); ++iParticle) {
        	// TODO: TAKE CARE for multiple copies from envelopes
            BondParticle3D<T,Descriptor>* bondParticle = castParticle3DToBondParticle3D(bondParticles[iParticle]);
            std::pair<plint,plint> pair0 = bondParticle->getCellIdVertexIdPair(0);
            std::pair<plint,plint> pair1 = bondParticle->getCellIdVertexIdPair(1);
            Particle3D<T,Descriptor>* p0=0;
            Particle3D<T,Descriptor>* p1=0;
            if ( particle0Map.find(pair0) != particle0Map.end() ) { p0 = particle0Map[pair0]; }
            if ( particle1Map.find(pair1) != particle1Map.end() ) { p1 = particle1Map[pair1]; }
            bondParticle->update(p0, p1);
            if (not bondField.getBondType().breakBond(bondParticle)) { // if breakBond is true, it breaks the bond as well: bondParticle->getTag()=-1
            	bondField.getBondType().applyForce(bondParticle);
            }
        }
        bondParticleField.removeParticles(bondParticleField.getBoundingBox(), -1);
        bondParticleField.advanceParticles(bondParticleField.getBoundingBox(), -1);
    };
private:
};


} // namespace trombocit


#include "BondFunctionals3D.hh"



#endif  // BOND_FUNCTIONALS_3D_H
