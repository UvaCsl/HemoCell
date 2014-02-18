#ifndef CELL_FIELD_3D_HH
#define CELL_FIELD_3D_HH

#include "cellField3D.h"

/* ******** MapVertexToParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CellField3D<T,Descriptor>::CellField3D(TriangleBoundary3D<T> const& Cells_, \
        plint numVerticesPerCell_, plint numTriangles_, T maxDiameter_,
        std::map<plint, Particle3D<T,Descriptor>*>& tagToParticle3D_) :
    Cells(Cells_), numVerticesPerCell(numVerticesPerCell_), numTriangles(numTriangles_), maxDiameter(maxDiameter_), tagToParticle3D(tagToParticle3D_) { }


template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::reduceQuantity1D(plint cellId, plint ccrId, T value, plint numParts)
{
    if (quantities1D[cellId].count(ccrId) == 0) { quantities1D[cellId][ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        T prValue = quantities1D[cellId][ccrId];
        if (0 == reductionType)      { quantities1D[cellId][ccrId] += value; }
        else if (1 == reductionType) { quantities1D[cellId][ccrId] = (prValue*particlesPerCellId[cellId] + value*numParts) * 1.0 / (particlesPerCellId[cellId] + numParts); }
        else if (2 == reductionType) { quantities1D[cellId][ccrId] = max(prValue, value); }
        else if (3 == reductionType) { quantities1D[cellId][ccrId] = min(prValue, value);  } // Min can be calculated from the inverse Max
    //        else if (4 == reductionType) { false; } // Std not implemented yet
    }
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::reduceQuantity3D(plint cellId, plint ccrId, Array<T,3> const& value, plint numParts) {
    plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
    if (quantities3D[cellId].count(ccrId) == 0) { quantities3D[cellId][ccrId] = value; }
    else {
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
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::reduceQuantityND(plint cellId, plint ccrId, std::vector<T> const& value, plint numParts) {
    plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
    if (quantitiesND[cellId].count(ccrId) == 0) { quantitiesND[cellId][ccrId] = value; }
    else {
        for (pluint iv = 0; iv < value.size(); ++iv) {
            T prValue = quantitiesND[cellId][ccrId][iv];
            if (0 == reductionType)      { quantitiesND[cellId][ccrId][iv] = prValue + value[iv]; }
            else if (1 == reductionType) { quantitiesND[cellId][ccrId][iv] = (prValue*particlesPerCellId[cellId] + numParts*value[iv]) * 1.0 / (particlesPerCellId[cellId] + numParts); }
            else if (2 == reductionType) { quantitiesND[cellId][ccrId][iv] = max(prValue, value[iv]); }
            else if (3 == reductionType) { quantitiesND[cellId][ccrId][iv] = min(prValue, value[iv]);  } // Min can be calculated from the inverse Max
    //            else if (4 == reductionType) { false; } // Std not implemented yet
        }
    }
}



template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::setCellIds(std::vector<plint> cellIds_) {
    cellIds.clear();
    cellIds = cellIds_;
} ;


template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::calcCellIds() {
    cellIds.clear();
    typename std::map<plint, std::map<plint, T > >::iterator iter1D;
    for (iter1D  = quantities1D.begin(); iter1D != quantities1D.end(); ++iter1D) {
        cellIds.push_back(iter1D->first);
    }
} ;

template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::clearQuantities() {
        cellIds.clear();
        quantities1D.clear();
        quantities3D.clear();
        quantitiesND.clear();
        particlesPerCellId.clear();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T,Descriptor>::clearQuantity(plint subscribedQuantity) {
        plint ccrId = subscribedQuantity;
        plint dim = ccrId%10; // Dimension: last digit
        for (pluint ci = 0; ci < cellIds.size(); ++ci) {
            plint cellId = cellIds[ci];
            if (dim == 1) {
                if (quantities1D.count(cellId) >0 ) { quantities1D[cellId].erase(ccrId); }
            } else if (dim == 3) {
                if (quantities3D.count(cellId) >0 ) { quantities3D[cellId].erase(ccrId); }
            } else {
                if (quantitiesND.count(cellId) >0 ) { quantitiesND[cellId].erase(ccrId); }
            }
        }
}

template<typename T, template<typename U> class Descriptor>
std::vector<plint> const& CellField3D<T,Descriptor>::getCellIds() {
        if (cellIds.size() == 0) {
            calcCellIds();
        };
        return cellIds;
} ;


/* ******** MapVertexToParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::MapVertexToParticle3D (
            TriangleBoundary3D<T> const& Cells_,
            std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_)
    : Cells(Cells_), tagToParticle3D(tagToParticle3D_)
{
    tagToParticle3D.clear();
}

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::~MapVertexToParticle3D()
{
}

template<typename T, template<typename U> class Descriptor>
MapVertexToParticle3D<T,Descriptor>::MapVertexToParticle3D (
            MapVertexToParticle3D<T,Descriptor> const& rhs)
    : Cells(rhs.Cells),
      tagToParticle3D(rhs.tagToParticle3D)
{ }


template<typename T, template<typename U> class Descriptor>
void MapVertexToParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(particleField.getBoundingBox(), found); // Gets the whole domain.
    tagToParticle3D.clear();
    Array<T,3> pos;
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        plint vertexId = nonTypedParticle->getTag();
        if (tagToParticle3D.count(vertexId) == 0) {
            tagToParticle3D[vertexId] = nonTypedParticle;
        } else if (contained(nonTypedParticle->getPosition(), domain)) {
            tagToParticle3D[vertexId] = nonTypedParticle;
        }
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
    return BlockDomain::bulk; // This changes in processGenericBlocks with getBoundingBox();
}

template<typename T, template<typename U> class Descriptor>
void MapVertexToParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
}





/* ******** MapVertexToParticle3D *********************************** */
//std::map<plint, Particle3D<T,Descriptor>*> * tagToParticle3D;
//std::map<plint, plint> numParticlesPerCellId;
//std::map<plint, plint> * cellIdToMeshCellId;
//std::stack<plint> * freeMeshCellIds;

template<typename T, template<typename U> class Descriptor>
void MeshToParticleField3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);
    numParticlesPerCellId.clear();

    /* Count cellIds of particles in current domain. */
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        numParticlesPerCellId[particle->get_cellId()] += 1;
    }

    /* If there are no particles of cellId in the domain and there is a connection with
     * the meshCellId, erase that connection. */
    std::map<plint, plint>::iterator iter = cellIdToMeshCellId[0].begin();
    for (; iter != cellIdToMeshCellId[0].end();) {
        plint cellId = iter->first;
        plint meshCellId = iter->second;
        if (numParticlesPerCellId.count(cellId)>0) {
            freeMeshCellIds[0].push(meshCellId);
            numParticlesPerCellId.erase(cellId);
            cellIdToMeshCellId[0].erase(iter++);
        } else {
            iter++;
        }
    }
    /* If there are particles of cellId in the domain and there is no connection with
     * the meshCellId, create that connection. */
    for (iter=numParticlesPerCellId.begin(); iter != numParticlesPerCellId.end(); ++iter) {
        plint cellId = iter->first;
        plint numParticles = iter->second;
        if (numParticles > 0) {
            cellIdToMeshCellId[0][cellId] = freeMeshCellIds[0].top();
            freeMeshCellIds[0].pop();
        }
    }
//    pcout << "ncl " << cellIdToMeshCellId[0].size() << std::endl;

    /* Commented because the rest are not yet implemented */
//    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
//        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
//        ImmersedCellParticle3D<T,Descriptor>* particle =
//            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
//        plint vertexId = particle->getTag();
//        tagToParticle3D[0][vertexId] = nonTypedParticle;
//    }

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
        Array<T,3> position(particle->get_pbcPosition());
        plint vertexId = particle->getTag();
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
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void CopyParticleToMeshVertex3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
}



#endif // CELL_FIELD_3D_HH
