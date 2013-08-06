#ifndef MESH_TO_PARTICLE_FIELD_3D_H
#define MESH_TO_PARTICLE_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include <vector>
#include <map>
#include <stack>

/*
 * Maps the vertices of the mesh (including envelopes) to an std::map<plint, Particle3D*>
 */
template<typename T, template<typename U> class Descriptor>
class MapVertexToParticle3D : public BoxProcessingFunctional3D
{
public:
    MapVertexToParticle3D (
            TriangleBoundary3D<T> const& triangleBoundary_,
            std::map<plint, Particle3D<T,Descriptor>*> & iVertexToParticle3D_);
    ~MapVertexToParticle3D();
    MapVertexToParticle3D(MapVertexToParticle3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MapVertexToParticle3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangleBoundary3D<T> const& triangleBoundary;
    std::map<plint, Particle3D<T,Descriptor>*> * iVertexToParticle3D;
};

//* std::map<plint, pluint> cellIdToMeshId
//    + cellIdToMeshId[cellId] = cellId * numberOfMeshVertices
//* std::stack <pluint> freeMeshIds;

//if (0 == cellIdToMeshId.count(particle.cellId)) {
//    if (freeMeshIds.size() > 0)
//        cellIdToMeshId[particle.cellId] = freeMeshIds.pop()
//    else
//        addCellToMesh()


template<typename T, template<typename U> class Descriptor>
class MeshParticleFieldConnection
{
public:
	MeshParticleFieldConnection(TriangleBoundary3D<T> const& triangleBoundary_);
	Particle3D<T,Descriptor>* getParticleFromVertex(plint iVertex);
	plint getVertexFromParticle(plint particleTag);
	plint getVertexFromParticle(Particle3D<T,Descriptor> const& particle);
	plint getCellIdFromVertex(plint iVertex);
	plint getTriangleIdFromParticle(plint iVertex);
private:
	plint numVertices, numTriangles;
	std::stack <pluint> freeMeshIds;
	std::map<plint, pluint> cellIdToMeshId;
    std::map<plint, Particle3D<T,Descriptor>*> * iVertexToParticle3D;
};


template<typename T, template<typename U> class Descriptor>
class CopyParticleToMeshVertex3D : public BoxProcessingFunctional3D
{
public:
    CopyParticleToMeshVertex3D(TriangularSurfaceMesh<T>& mesh_);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CopyParticleToMeshVertex3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangularSurfaceMesh<T>& mesh;
};


#include "meshToParticleField3D.hh"

#endif // MESH_TO_PARTICLE_FIELD_3D_H
