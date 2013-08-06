#ifndef MESH_TO_PARTICLE_FIELD_3D_HH
#define MESH_TO_PARTICLE_FIELD_3D_HH

#include "meshToParticleField3D.h"

/* ******** MapVertexToParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::MapVertexToParticle3D (
            TriangleBoundary3D<T> const& triangleBoundary_,
            std::map<plint, Particle3D<T,Descriptor>*> & iVertexToParticle3D_)
    : triangleBoundary(triangleBoundary_)
{
    	iVertexToParticle3D = &iVertexToParticle3D_;
    	iVertexToParticle3D[0].clear();
}

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::~MapVertexToParticle3D()
{
}

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::MapVertexToParticle3D (
            MapVertexToParticle3D<T,Descriptor> const& rhs)
    : triangleBoundary(rhs.triangleBoundary),
      iVertexToParticle3D(rhs.iVertexToParticle3D)
{ }


template<typename T, template<typename U> class Descriptor>
void MapVertexToParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getVertexId();
        iVertexToParticle3D[0][vertexId] = nonTypedParticle;
    }


}

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>*
    MapVertexToParticle3D<T,Descriptor>::clone() const
{
    return new MapVertexToParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void MapVertexToParticle3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT MapVertexToParticle3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void MapVertexToParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
}




/* ******** CopyParticleToMeshVertex3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CopyParticleToMeshVertex3D<T,Descriptor>::CopyParticleToMeshVertex3D (
        TriangularSurfaceMesh<T>& mesh_)
    : mesh(mesh_)
{ }

template<typename T, template<typename U> class Descriptor>
void CopyParticleToMeshVertex3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    // Manually extend the domain to the full envelope.
    domain = particleField.getBoundingBox();

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (found[iParticle]);
        Array<T,3> position(particle.getPosition());
        plint vertexId = particle.getVertexId();
        mesh.replaceVertex(vertexId, position);
    }
}

template<typename T, template<typename U> class Descriptor>
CopyParticleToMeshVertex3D<T,Descriptor>* CopyParticleToMeshVertex3D<T,Descriptor>::clone() const {
    return new CopyParticleToMeshVertex3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT CopyParticleToMeshVertex3D<T,Descriptor>::appliesTo() const {
    // The data processor acts on envelope too, but extension to the envelope
    //   is done manually in processGenericBlocks.
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void CopyParticleToMeshVertex3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
}



#endif // MESH_TO_PARTICLE_FIELD_3D_HH
