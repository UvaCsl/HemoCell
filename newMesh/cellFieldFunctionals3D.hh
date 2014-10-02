#ifndef CELL_FIELD_FUNCTIONALS_3D_HH
#define CELL_FIELD_FUNCTIONALS_3D_HH
#include "cellFieldFunctionals3D.h"

/* ******** ComputeCellForce3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeCellForce3D<T,Descriptor>::ComputeCellForce3D (ConstitutiveModel<T,Descriptor>* cellModel_, std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_)
: cellModel(cellModel_), cellIdToCell3D(cellIdToCell3D_) { }

template<typename T, template<typename U> class Descriptor>
ComputeCellForce3D<T,Descriptor>::~ComputeCellForce3D()
{
    delete cellModel;
}

template<typename T, template<typename U> class Descriptor>
ComputeCellForce3D<T,Descriptor>::ComputeCellForce3D (
            ComputeCellForce3D<T,Descriptor> const& rhs)
    : cellModel(rhs.cellModel), cellIdToCell3D(rhs.cellIdToCell3D)
{ }


template<typename T, template<typename U> class Descriptor>
void ComputeCellForce3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    typename std::map<plint, Cell3D<T,Descriptor>  >::iterator iter;
    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
        cellModel->computeCellForce(iter->second);
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeCellForce3D<T,Descriptor>*
    ComputeCellForce3D<T,Descriptor>::clone() const
{
    return new ComputeCellForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellForce3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeCellForce3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellForce3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}




/* ******** FluidVelocityToImmersedCell3D *********************************** */

template<typename T, template<typename U> class Descriptor>
FluidVelocityToImmersedCell3D<T,Descriptor>::FluidVelocityToImmersedCell3D (
        plint ibmKernel_) : ibmKernel(ibmKernel_) {};

template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCell3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);
    /* Not used */
//    Dot3D offset = computeRelativeDisplacement(particleField, fluid);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    std::vector<Cell<T,Descriptor>*> cells;
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iParticle]);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        Array<T,3> velocity; velocity.resetToZero();
        cellPos = particle->getIBMcoordinates();
        weights = particle->getIBMweights();
        if (cellPos.size() == 0) {
        	interpolationCoefficients(fluid, position, cellPos, weights, ibmKernel);
        }
        particle->get_v().resetToZero();
        for (pluint iCell=0; iCell < weights.size(); ++iCell) {
            velocity.resetToZero();
            Dot3D cellPosition = cellPos[iCell];
//            Dot3D cellPosition = cellPos[iCell] + offset;
            fluid.get(cellPosition.x, cellPosition.y, cellPosition.z).computeVelocity(velocity);
            particle->get_v() += weights[iCell] * velocity;
        }
        particle->get_vPrevious() = particle->get_v();
    }
}

template<typename T, template<typename U> class Descriptor>
FluidVelocityToImmersedCell3D<T,Descriptor>* FluidVelocityToImmersedCell3D<T,Descriptor>::clone() const {
    return new FluidVelocityToImmersedCell3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCell3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FluidVelocityToImmersedCell3D<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}



/* ******** ForceToFluid3D *********************************** */
template<typename T, template<typename U> class Descriptor>
ForceToFluid3D<T,Descriptor>::ForceToFluid3D (
        plint ibmKernel_) : ibmKernel(ibmKernel_) {};

template<typename T, template<typename U> class Descriptor>
void ForceToFluid3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);
//    Dot3D offset = computeRelativeDisplacement(particleField, fluid); NOT USED

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    Cell<T,Descriptor>* cell;
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iParticle]);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        cellPos = particle->getIBMcoordinates();
        weights = particle->getIBMweights();
        interpolationCoefficients(fluid, position, cellPos, weights, ibmKernel);
        Array<T,3> elasticForce = particle->get_force();
        // pcout << "elastic force: (" << elasticForce[0] << ", "<< elasticForce[1] << ", "<< elasticForce[2] << ")\n";
        for (pluint iCell = 0; iCell < weights.size(); ++iCell) {
            Dot3D cellPosition = cellPos[iCell];
            cell = &fluid.get(cellPosition.x, cellPosition.y, cellPosition.z);
            T *locForce = cell->getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
            for (pluint iA = 0; iA < 3; ++iA) {
                locForce[iA] += (weights[iCell])*elasticForce[iA];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ForceToFluid3D<T,Descriptor>* ForceToFluid3D<T,Descriptor>::clone() const {
    return new ForceToFluid3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ForceToFluid3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
    modified[1] = modif::staticVariables; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ForceToFluid3D<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** FillCellMap *********************************** */
template<typename T, template<typename U> class Descriptor>
void FillCellMap<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    cellIdToCell3D.clear();

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(particleField.getBoundingBox(), found); // Gets the whole domain.

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle = 
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (found[iParticle]);

        plint cellId = particle->get_cellId();
        typename std::map<plint, Cell3D<T,Descriptor>  >::iterator iter = cellIdToCell3D.find(cellId);
        if (iter != cellIdToCell3D.end()) {
            iter->second = Cell3D<T,Descriptor>(mesh, cellId);
        }
        iter->second.push_back(particle);
    }

    typename std::map<plint, Cell3D<T,Descriptor> >::iterator iter;
    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
        (iter->second).close();
    }
}



template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::FillCellMap (
	TriangularSurfaceMesh<T> const& mesh_,
	std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_)
: mesh(mesh_), cellIdToCell3D(cellIdToCell3D_) { }


template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::~FillCellMap() { }


template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::FillCellMap (
            FillCellMap<T,Descriptor> const& rhs)
    : mesh(rhs.mesh), cellIdToCell3D(rhs.cellIdToCell3D) { }


template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>*
    FillCellMap<T,Descriptor>::clone() const
{   return new FillCellMap<T,Descriptor>(*this);    }


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FillCellMap<T,Descriptor>::appliesTo() const 
{    return BlockDomain::bulk;   } // This changes in processGenericBlocks with getBoundingBox();


template<typename T, template<typename U> class Descriptor>
void FillCellMap<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{    modified[0] = modif::nothing;  } // Particle field.



/* ******** ComputeRequiredQuantities *********************************** */




template<typename T, template<typename U> class Descriptor>
void ComputeRequiredQuantities<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& reductionParticleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles); // Gets particle only from the bulk

    std::vector<Particle3D<T,Descriptor>*> reductionParticles;
    reductionParticleField.removeParticles(reductionParticleField.getBoundingBox());

    typename std::map<plint, Cell3D<T,Descriptor>  >::iterator iter;
    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
        iter->second.clearReducer();
    }
    std::map<plint, pluint> particlesPerCellId;
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        for (std::vector<plint>::iterator i = ccrRequirements.begin(); i != ccrRequirements.end(); ++i)
        {
            plint ccrId = *i;
            plint cellId = castParticleToICP3D(particles[iParticle])->get_cellId();

            typename std::map<plint, Cell3D<T,Descriptor>  >::iterator iter = cellIdToCell3D.find(cellId);
            if (iter != cellIdToCell3D.end()) {
                iter->second.computeCCRQuantities(ccrId, particles[iParticle]);
            }

            if (particlesPerCellId.count(cellId) > 0) { particlesPerCellId[cellId] += 1; }
            else  { particlesPerCellId[cellId] = 1; }
        }
    }

    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
        plint cellId = iter->first;
        Cell3D<T,Descriptor> & cell = iter->second;
        cell.closeCCRQuantities();
        Array<T,3> const& vertex = cell.get3D(CCR_NO_PBC_POSITION_MEAN);
        ReductionParticle3D<T,Descriptor>* rParticle = new ReductionParticle3D<T,Descriptor>(cellId, vertex);
        rParticle->get_nParticles() = particlesPerCellId[cellId];
        rParticle->updateCQH(cell);
        reductionParticleField.addParticle(reductionParticleField.getBoundingBox(), rParticle);
    }
}



template<typename T, template<typename U> class Descriptor>
ComputeRequiredQuantities<T,Descriptor>::ComputeRequiredQuantities (
    std::vector<plint> ccrRequirements_, 
    std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_)
: ccrRequirements(ccrRequirements_), cellIdToCell3D(cellIdToCell3D_) { }


template<typename T, template<typename U> class Descriptor>
ComputeRequiredQuantities<T,Descriptor>::~ComputeRequiredQuantities() { }


template<typename T, template<typename U> class Descriptor>
ComputeRequiredQuantities<T,Descriptor>::ComputeRequiredQuantities (
            ComputeRequiredQuantities<T,Descriptor> const& rhs)
    : ccrRequirements(rhs.ccrRequirements), cellIdToCell3D(rhs.cellIdToCell3D) { }


template<typename T, template<typename U> class Descriptor>
ComputeRequiredQuantities<T,Descriptor>*
    ComputeRequiredQuantities<T,Descriptor>::clone() const
{ return new ComputeRequiredQuantities<T,Descriptor>(*this); }


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeRequiredQuantities<T,Descriptor>::appliesTo() const
{    return BlockDomain::bulk;    }


template<typename T, template<typename U> class Descriptor>
void ComputeRequiredQuantities<T,Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const  {
    modified[0] = modif::nothing; // Particle field.
    modified[1] = modif::allVariables; // Reduction particles;
}



/* ******** SyncCellQuantities *********************************** */

template<typename T, template<typename U> class Descriptor>
SyncCellQuantities<T,Descriptor>::SyncCellQuantities (std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_)
: cellIdToCell3D(cellIdToCell3D_) { }

template<typename T, template<typename U> class Descriptor>
SyncCellQuantities<T,Descriptor>::SyncCellQuantities (SyncCellQuantities<T,Descriptor> const& rhs)
: cellIdToCell3D(rhs.cellIdToCell3D) { }



template<typename T, template<typename U> class Descriptor>
void SyncCellQuantities<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
        PLB_PRECONDITION( blocks.size()==2 );
        ParticleField3D<T,Descriptor>& reductionParticleField
            = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);
        std::vector<Particle3D<T,Descriptor>*> reductionParticles;
        reductionParticleField.findParticles(reductionParticleField.getBoundingBox(), reductionParticles);
        std::map< plint, std::map<plint, bool> > particleInProcessorAndCellid; // processorCellIdMap[processor][cellId] = true
        for (pluint iA = 0; iA < reductionParticles.size(); ++iA) {
            ReductionParticle3D<T,Descriptor>* particle =
                    dynamic_cast<ReductionParticle3D<T,Descriptor>*> (reductionParticles[iA]);

            plint cellId = particle->get_cellId();
            plint processor = particle->get_processor();
            plint nParticles = particle->get_nParticles();
            particleInProcessorAndCellid[processor]; // Create "processor" entry
            typename std::map<plint, Cell3D<T,Descriptor>  >::iterator iter = cellIdToCell3D.find(cellId);
            if (iter != cellIdToCell3D.end()) {
                Cell3D<T,Descriptor> & chq = iter->second;
                if (particleInProcessorAndCellid[processor].count(cellId) == 0
                                         &&  processor != particle->getMpiProcessor()) {
                    particleInProcessorAndCellid[processor][cellId] = true;
                    std::map<plint, T > const& q1d = particle->getQuantities1D();
                    std::map<plint, Array<T,3> > const& q3d = particle->getQuantities3D();
                    std::map<plint, std::vector<T> > const& qNd = particle->getQuantitiesND();
                    typename std::map<plint, T >::const_iterator iter1D;
                    for (iter1D  = q1d.begin(); iter1D != q1d.end(); ++iter1D) {
                        chq.reduceQuantity(iter1D->first, iter1D->second, nParticles);
                    }
                    typename std::map<plint, Array<T,3> >::const_iterator iter3D;
                    for (iter3D  = q3d.begin(); iter3D != q3d.end(); ++iter3D) {
                        chq.reduceQuantity(iter3D->first, iter3D->second, nParticles);
                    }
                    typename std::map<plint, std::vector<T> >::const_iterator iterND;
                    for (iterND  = qNd.begin(); iterND != qNd.end(); ++iterND) {
                        chq.reduceQuantity(iterND->first, iterND->second, nParticles);
                    }
                    chq.getParticlesPerCellId() += nParticles;
                }
            }
        }
        reductionParticleField.removeParticles(reductionParticleField.getBoundingBox(), -1);
}

template<typename T, template<typename U> class Descriptor>
SyncCellQuantities<T,Descriptor>*
    SyncCellQuantities<T,Descriptor>::clone() const
{
    return new SyncCellQuantities<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SyncCellQuantities<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor>
void SyncCellQuantities<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
    modified[1] = modif::nothing; // Reduction Particle field.
}



#endif  // CELL_FIELD_FUNCTIONALS_3D_HH

