#ifndef MESH_TO_PARTICLE_FIELD_3D_H
#define MESH_TO_PARTICLE_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include "cellReductionTypes.h"
#include <vector>
#include <map>
#include <stack>

template<typename T, template<typename U> class Descriptor>
class CellFieldQuantityHolder
{
public:
    CellFieldQuantityHolder(TriangleBoundary3D<T> const& triangleBoundary_, plint numVerticesPerCell_,
            std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D_) :
        triangleBoundary(triangleBoundary_), numVerticesPerCell(numVerticesPerCell_), tagToParticle3D(tagToParticle3D_) { }
    plint & getNumVerticesPerCell() { return numVerticesPerCell; }
    plint & getNumMeshCells() { return numMeshCells; }
    std::map<plint, Particle3D<T,Descriptor>*>& get_tagToParticle3D() { return tagToParticle3D; }
public:
    plint getVertexId(plint particleTag) {
        plint cellId = particleTag/numVerticesPerCell;
        plint vertexIdOnCell = (particleTag%numVerticesPerCell);
        return vertexIdOnCell + numVerticesPerCell * cellIdToMeshCellId[cellId];
    }

private:
    TriangleBoundary3D<T> const& triangleBoundary;
    plint numMeshCells, numVerticesPerCell;
    std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D;
public:
    void clearQuantities() { cellIds.clear(); quantities1D.clear(); quantities3D.clear(); quantitiesND.clear(); particlesPerCellId.clear(); }
    plint getNumCells() { return cellIds.size(); }
    std::vector<plint> const& getCellIds() { if (cellIds.size() == 0) { calcCellIds(); }; return cellIds; } ;
    void calcCellIds() {
        cellIds.clear();
        typename std::map<plint, std::map<plint, T > >::iterator iter1D;
        for (iter1D  = quantities1D.begin(); iter1D != quantities1D.end(); ++iter1D) {
            cellIds.push_back(iter1D->first);
        }
    } ;
    T getVolume(plint cellId) { return quantities1D[cellId][CCR_VOLUME];} ;
    T getSurface(plint cellId) { return quantities1D[cellId][CCR_SURFACE];} ;
    T getEnergy(plint cellId) { return quantities1D[cellId][CCR_ENERGY];} ;
    T getMeanAngle(plint cellId) { return quantities1D[cellId][CCR_ANGLE_MEAN];} ;
    T getMeanEdgeLength(plint cellId) { return quantities1D[cellId][CCR_EDGE_DISTANCE_MEAN];} ;
    Array<T,3> const& getPosition(plint cellId) { return quantities3D[cellId][CCR_POSITION_MEAN];  } ;
    Array<T,3> const& getVelocity(plint cellId) { return quantities3D[cellId][CCR_VELOCITY_MEAN];} ;
    std::vector<T> const& getInertia(plint cellId) { return quantitiesND[cellId][CCR_INERTIA];} ;
//    Array<Array<T,3>,3> getInertia(plint cellId) ;

    // 0 -- Sum, 1 -- Mean, 2 -- Max, 3 -- Min, 4 -- Std
    void reduceQuantity1D(plint cellId, plint ccrId, T value, plint numParts=0) {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        T prValue = quantities1D[cellId][ccrId];
        if (0 == reductionType)      { quantities1D[cellId][ccrId] += value; }
        else if (1 == reductionType) { quantities1D[cellId][ccrId] = (prValue*particlesPerCellId[cellId] + value*numParts) * 1.0 / (particlesPerCellId[cellId] + numParts); }
        else if (2 == reductionType) { quantities1D[cellId][ccrId] = max(prValue, value); }
        else if (3 == reductionType) { quantities1D[cellId][ccrId] = min(prValue, value);  } // Min can be calculated from the inverse Max
//        else if (4 == reductionType) { false; } // Std not implemented yet
    }
    void reduceQuantity3D(plint cellId, plint ccrId, Array<T,3> const& value, plint numParts=0) {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        Array<T,3> prValue = quantities3D[cellId][ccrId];
        if (0 == reductionType)      { quantities3D[cellId][ccrId] = prValue + value; }
        else if (1 == reductionType) {
            plint prNumParts =  particlesPerCellId[cellId];
            quantities3D[cellId][ccrId] = prNumParts * 1.0 * prValue + numParts * 1.0 * value;
            T factr = 1.0 / (prNumParts + numParts);
            quantities3D[cellId][ccrId] = quantities3D[cellId][ccrId] * factr;
        }
        else if (2 == reductionType) {
            quantities3D[cellId][ccrId] = Array<T,3>( max(prValue[0], value[0]),
                                                      max(prValue[1], value[1]),
                                                      max(prValue[2], value[2]) );
        }
        else if (3 == reductionType) {
            quantities3D[cellId][ccrId] = Array<T,3>( min(prValue[0], value[0]),
                                                      min(prValue[1], value[1]),
                                                      min(prValue[2], value[2]) );
        }
//        else if (4 == reductionType) { false; } // Std not implemented yet
    }
    void reduceQuantityND(plint cellId, plint ccrId, std::vector<T> const& value, plint numParts=0) {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        for (pluint iv = 0; iv < value.size(); ++iv) {
            T prValue = quantitiesND[cellId][ccrId][iv];
            if (0 == reductionType)      { quantitiesND[cellId][ccrId][iv] = prValue + value[iv]; }
            else if (1 == reductionType) { quantitiesND[cellId][ccrId][iv] = (prValue*particlesPerCellId[cellId] + numParts*value[iv]) * 1.0 / (particlesPerCellId[cellId] + numParts); }
            else if (2 == reductionType) { quantitiesND[cellId][ccrId][iv] = max(prValue, value[iv]); }
            else if (3 == reductionType) { quantitiesND[cellId][ccrId][iv] = min(prValue, value[iv]);  } // Min can be calculated from the inverse Max
//            else if (4 == reductionType) { false; } // Std not implemented yet
        }
    }

    std::map<plint, plint> & getParticlesPerCellId() { return particlesPerCellId ; }
private:
    std::map<plint, std::map<plint, T > > quantities1D; // quantities1D[cellId][CCR_EDGE_DISTANCE_STD]= MEAN_EDGE_DISTANCE
    std::map<plint, std::map<plint, Array<T,3> > > quantities3D;
    std::map<plint, std::map<plint, std::vector<T> > > quantitiesND;
    std::vector<plint> cellIds;
    std::map<plint, plint> particlesPerCellId; // particlesPerCellId[cellId]= nParticles;
private: // NOT IMPLEMENTED YET;
    std::map<plint, plint> cellIdToMeshCellId;
    std::stack<plint> freeMeshCellIds;
    std::stack<plint> & getFreeMeshCellIds() { return freeMeshCellIds; }
    std::map<plint, plint> & getCellIdToMeshCellId() { return cellIdToMeshCellId; }
    void init() {
        plint nVertices = triangleBoundary.getMesh().getNumVertices();
        numMeshCells = nVertices / numVerticesPerCell;
        PLB_ASSERT(numMeshCells * numVerticesPerCell == nVertices);
        for (plint i=0; i < numMeshCells; i++) {
            freeMeshCellIds.push(i);
        }
    }
};

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
            TriangleBoundary3D<T> const& triangleBoundary_,
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
    TriangleBoundary3D<T> const& triangleBoundary;
    std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D;
};



template<typename T, template<typename U> class Descriptor>
class MeshToParticleField3D : public BoxProcessingFunctional3D
{
public:
    MeshToParticleField3D(CellFieldQuantityHolder<T,Descriptor> & cellInfoHolder_,
            std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_) :
                tagToParticle3D(&tagToParticle3D_),
                verticesPerCell(cellInfoHolder_.getNumVerticesPerCell()),
                numParticlesPerCellId(),
                cellIdToMeshCellId(&cellInfoHolder_.getCellIdToMeshCellId()),
                freeMeshCellIds(&cellInfoHolder_.getFreeMeshCellIds())	{	};
    MeshToParticleField3D(CellFieldQuantityHolder<T,Descriptor> & cellInfoHolder_) :
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


#include "meshToParticleField3D.hh"

#endif // MESH_TO_PARTICLE_FIELD_3D_H
