#ifndef MESH_TO_PARTICLE_FIELD_3D_H
#define MESH_TO_PARTICLE_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include <vector>
#include <map>
#include <stack>

template<typename T, template<typename U> class Descriptor>
class CellFieldQuantityHolder
{
public:
	CellFieldQuantityHolder(TriangleBoundary3D<T> const& triangleBoundary_, plint numVerticesPerCell_) :
		numVerticesPerCell(numVerticesPerCell_) {
		plint nVertices = triangleBoundary_.getMesh().getNumVertices();
		numMeshCells = nVertices / numVerticesPerCell;
		PLB_ASSERT(numMeshCells * numVerticesPerCell == nVertices);
		for (plint i=0; i < numMeshCells; i++) {
			freeMeshCellIds.push(i);
		}
	} ;
	CellFieldQuantityHolder(plint numMeshCells_, plint numVerticesPerCell_) :
		numMeshCells(numMeshCells_), numVerticesPerCell(numVerticesPerCell_) {
		for (plint i=0; i < numMeshCells; i++) {
			freeMeshCellIds.push(i);
		}
	} ;
	plint & getNumVerticesPerCell() { return numVerticesPerCell; }
	plint & getNumMeshCells() { return numMeshCells; }
	std::stack<plint> & getFreeMeshCellIds() { return freeMeshCellIds; }
	std::map<plint, plint> & getCellIdToMeshCellId() { return cellIdToMeshCellId; }
public:
	plint getVertexId(plint particleTag) {
		plint cellId = particleTag/numVerticesPerCell;
		plint vertexIdOnCell = (particleTag%numVerticesPerCell);
		return vertexIdOnCell + numVerticesPerCell * cellIdToMeshCellId[cellId];
	}
private:
	plint numMeshCells, numVerticesPerCell;
	std::map<plint, plint> cellIdToMeshCellId;
	std::stack<plint> freeMeshCellIds;
//    std::map<plint, Particle3D<T,Descriptor>*> * iVertexToParticle3D;
};

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
    ~MapVertexToParticle3D() { };
    MapVertexToParticle3D(MapVertexToParticle3D<T,Descriptor> const& rhs) :
    	triangleBoundary(rhs.triangleBoundary), iVertexToParticle3D(rhs.iVertexToParticle3D) { }
    virtual MapVertexToParticle3D<T,Descriptor>* clone() const { return new MapVertexToParticle3D<T,Descriptor>(*this); }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const { isWritten[0] = false; }
    virtual BlockDomain::DomainT appliesTo() const  { return BlockDomain::bulkAndEnvelope; }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const  { modified[0] = modif::nothing; }
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
private:
    TriangleBoundary3D<T> const& triangleBoundary;
    std::map<plint, Particle3D<T,Descriptor>*> * iVertexToParticle3D;
};


template<typename T, template<typename U> class Descriptor>
class MeshToParticleField3D : public BoxProcessingFunctional3D
{
public:
	MeshToParticleField3D(CellFieldQuantityHolder<T,Descriptor> & cellInfoHolder_,
						   std::map<plint, Particle3D<T,Descriptor>*> & iVertexToParticle3D_) :
							   iVertexToParticle3D(&iVertexToParticle3D_),
							   verticesPerCell(cellInfoHolder_.getNumVerticesPerCell()),
							   numParticlesPerCellId(),
							   cellIdToMeshCellId(&cellInfoHolder_.getCellIdToMeshCellId()),
							   freeMeshCellIds(&cellInfoHolder_.getFreeMeshCellIds())	{	};
	MeshToParticleField3D(CellFieldQuantityHolder<T,Descriptor> & cellInfoHolder_) :
							   iVertexToParticle3D(),
							   verticesPerCell(cellInfoHolder_.getNumVerticesPerCell()),
							   numParticlesPerCellId(),
							   cellIdToMeshCellId(&cellInfoHolder_.getCellIdToMeshCellId()),
							   freeMeshCellIds(&cellInfoHolder_.getFreeMeshCellIds())	{	};
    ~MeshToParticleField3D() { } ;
    MeshToParticleField3D(MeshToParticleField3D<T,Descriptor> const& rhs) :
		   iVertexToParticle3D(rhs.iVertexToParticle3D),
		   verticesPerCell(rhs.verticesPerCell),
		   numParticlesPerCellId(rhs.numParticlesPerCellId),
		   cellIdToMeshCellId(rhs.cellIdToMeshCellId),
		   freeMeshCellIds(rhs.freeMeshCellIds) {    };
    virtual MeshToParticleField3D<T,Descriptor>* clone() const { return new MeshToParticleField3D<T,Descriptor>(*this); }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const { isWritten[0] = false; }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; }
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
private:
    std::map<plint, Particle3D<T,Descriptor>*> * iVertexToParticle3D;
	plint verticesPerCell;
    std::map<plint, plint> numParticlesPerCellId;
	std::map<plint, plint> * cellIdToMeshCellId;
	std::stack<plint> * freeMeshCellIds;
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
