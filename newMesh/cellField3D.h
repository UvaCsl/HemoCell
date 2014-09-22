#ifndef CELL_FIELD_3D_H
#define CELL_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cell3D.h"

#include "cellReductionTypes.h"
#include <vector>
#include <map>
#include <stack>

using namespace plb;
using namespace std;

template<typename T, template<typename U> class Descriptor>
class CellField3D
{
public:
    CellField3D(TriangleBoundary3D<T> const& Cells_,
            plint numVerticesPerCell_, plint numTriangles_, T maxDiameter_,
            std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D_);
    CellField3D(TriangleBoundary3D<T> const& Cells_, T dx_, T dt_,
            plint numVerticesPerCell_, plint numTriangles_, T maxDiameter_,
            std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D_);
    virtual ~CellField3D() { };
    void write(plint iter=0);
    void writeVTK(plint iter=0);
    Cell3D<T,Descriptor> operator[](plint cellId) { return CellIdToCell3D[cellId]; }

private:
    TriangleBoundary3D<T> const& Cells;
    std::map<plint, Cell3D<T,Descriptor> > CellIdToCell3D;
    T dx, dt;
    plint numMeshCells, numVerticesPerCell, numTriangles;
    T maxDiameter;
    std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D;
    std::vector<plint> cellIds;

    std::map<plint, plint> particlesPerCellId; // particlesPerCellId[cellId]= nParticles;
    std::map<plint, std::map<plint, T > > quantities1D; // quantities1D[cellId][CCR_EDGE_DISTANCE_STD]= MEAN_EDGE_DISTANCE
    std::map<plint, std::map<plint, Array<T,3> > > quantities3D;
    std::map<plint, std::map<plint, std::vector<T> > > quantitiesND;
    void calcCellIds();
    std::vector<plint> ccrIds1D, ccrIds3D, ccrIdsND;
    std::vector<Particle3D<T,Descriptor>*> cellReductionParticles;
public:
    std::vector<Particle3D<T,Descriptor>*> & getCellReductionParticles() { return cellReductionParticles ; }
    plint & getNumVerticesPerCell() { return numVerticesPerCell; }
    plint & getNumMeshCells() { return numMeshCells; }
    std::map<plint, Particle3D<T,Descriptor>*>& get_tagToParticle3D() { return tagToParticle3D; }
    plint getNumCells() { return cellIds.size(); }
    T getMaxDiameter() { return maxDiameter; }
    void setMaxDiameter(T maxDiameter_) { maxDiameter = maxDiameter_; }
public:
    void clearQuantities();
    void clearQuantity(plint subscribedQuantity);
    /* Get Class cointainers */
    std::vector<plint> const& getCellIds();
    void setCellIds(std::vector<plint> cellIds_);
    std::map<plint, plint> & getParticlesPerCellId() { return particlesPerCellId ; }


    // 0 -- Sum, 1 -- Mean, 2 -- Max, 3 -- Min, 4 -- Std
    void reduceQuantity1D(plint cellId, plint ccrId, T value, plint numParts=0) ;
    void reduceQuantity3D(plint cellId, plint ccrId, Array<T,3> const& value, plint numParts=0);
    void reduceQuantityND(plint cellId, plint ccrId, std::vector<T> const& value, plint numParts=0) ;
};


//private: // NOT IMPLEMENTED YET;
//    std::map<plint, plint> cellIdToMeshCellId;
//    std::stack<plint> freeMeshCellIds;
//    std::stack<plint> & getFreeMeshCellIds() { return freeMeshCellIds; }
//    std::map<plint, plint> & getCellIdToMeshCellId() { return cellIdToMeshCellId; }
//    void init() {
//        plint nVertices = Cells.getMesh().getNumVertices();
//        numMeshCells = nVertices / numVerticesPerCell;
//        PLB_ASSERT(numMeshCells * numVerticesPerCell == nVertices);
//        for (plint i=0; i < numMeshCells; i++) {
//            freeMeshCellIds.push(i);
//        }
//    }
//    plint getVertexId(plint particleTag) {
//        plint cellId = particleTag/numVerticesPerCell;
//        plint vertexIdOnCell = (particleTag%numVerticesPerCell);
//        return vertexIdOnCell + numVerticesPerCell * cellIdToMeshCellId[cellId];
//    }
//
//};

/*
 * Maps the vertices of the mesh (including envelopes) to an std::map<plint, Particle3D*>
 */
/*
 * Maps the vertices of the mesh (including envelopes) to an std::map<plint, Particle3D*>
 */
template<typename T, template<typename U> class Descriptor>
class MapVertexToParticle3D : public BoxProcessingFunctional3D
{
public:
    MapVertexToParticle3D (
            TriangleBoundary3D<T> const& Cells_,
            std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_);
    virtual ~MapVertexToParticle3D();
    MapVertexToParticle3D(MapVertexToParticle3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MapVertexToParticle3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangleBoundary3D<T> const& Cells;
    std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D;
};



template<typename T, template<typename U> class Descriptor>
class MeshToParticleField3D : public BoxProcessingFunctional3D
{
public:
    MeshToParticleField3D(CellField3D<T,Descriptor> & cellInfoHolder_,
            std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_) :
                tagToParticle3D(&tagToParticle3D_),
                verticesPerCell(cellInfoHolder_.getNumVerticesPerCell()),
                numParticlesPerCellId(),
                cellIdToMeshCellId(&cellInfoHolder_.getCellIdToMeshCellId()),
                freeMeshCellIds(&cellInfoHolder_.getFreeMeshCellIds())	{	};
    MeshToParticleField3D(CellField3D<T,Descriptor> & cellInfoHolder_) :
                tagToParticle3D(),
                verticesPerCell(cellInfoHolder_.getNumVerticesPerCell()),
                numParticlesPerCellId(),
                cellIdToMeshCellId(&cellInfoHolder_.getCellIdToMeshCellId()),
                freeMeshCellIds(&cellInfoHolder_.getFreeMeshCellIds())	{	};
    virtual ~MeshToParticleField3D() { } ;
    MeshToParticleField3D(MeshToParticleField3D<T,Descriptor> const& rhs) :
        tagToParticle3D(rhs.tagToParticle3D),
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
    std::map<plint, Particle3D<T,Descriptor>*> * tagToParticle3D;
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
    virtual ~CopyParticleToMeshVertex3D() { };
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CopyParticleToMeshVertex3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangularSurfaceMesh<T>& mesh;
};


#include "cellField3D.hh"

#endif // CELL_FIELD_3D
