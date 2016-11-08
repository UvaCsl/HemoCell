#ifndef CELL_3D_HH
#define CELL_3D_HH
#include "cell3D.h"

/********************* CellQuantityHolder *********************/

template<typename T>
CellQuantityHolder<T>::CellQuantityHolder(CellQuantityHolder<T> const& rhs):
    cellId(rhs.cellId),particlesPerCellId(rhs.particlesPerCellId),
    quantities1D(rhs.quantities1D),quantities3D(rhs.quantities3D),
    quantitiesND(rhs.quantitiesND),scalar_ccrIds(rhs.scalar_ccrIds),
    vector_ccrIds(rhs.vector_ccrIds),tensor_ccrIds(rhs.tensor_ccrIds)
    {    };

template<typename T>
CellQuantityHolder<T>& CellQuantityHolder<T>::operator=(CellQuantityHolder<T> const& rhs){
    cellId=rhs.cellId;
    particlesPerCellId=rhs.particlesPerCellId;
    quantities1D=rhs.quantities1D;
    quantities3D=rhs.quantities3D;
    quantitiesND=rhs.quantitiesND;
    scalar_ccrIds=rhs.scalar_ccrIds;
    vector_ccrIds=rhs.vector_ccrIds;
    tensor_ccrIds=rhs.tensor_ccrIds;
    return *this;
};

template<typename T>
void CellQuantityHolder<T>::clearQuantities(std::vector<plint> ccrIds) {
    for (pluint i=0; i<ccrIds.size(); ++i) {
        plint ccrId = ccrIds[i];
        plint d = getReductionDimension(ccrId);
        if (d==1) { quantities1D.erase(ccrId); }
        else if (d==3) { quantities3D.erase(ccrId); }
        else { quantitiesND.erase(ccrId); }
    }
}

template<typename T>
void CellQuantityHolder<T>::clearQuantities() {
    quantities1D.clear();
    quantities3D.clear();
    quantitiesND.clear();
    scalar_ccrIds.clear(); vector_ccrIds.clear(); tensor_ccrIds.clear();
    numParts_ccrId.clear();
}

template<typename T>
void CellQuantityHolder<T>::updateCQH(CellQuantityHolder<T> const& cqh) {
    this->clearQuantities();
    particlesPerCellId = cqh.getParticlesPerCellId() ;
    quantities1D = cqh.getQuantities1D();
    quantities3D = cqh.getQuantities3D();
    quantitiesND = cqh.getQuantitiesND();
    this->make_ccrId_List();
}

template<typename T>
void CellQuantityHolder<T>::updateCQH(CellQuantityHolder<T> const& cqh, std::vector<plint> const& ccreqIds) {
    this->clearQuantities();
    particlesPerCellId = cqh.getParticlesPerCellId() ;
    for (pluint i=0; i<ccreqIds.size(); ++i) {
        plint ccrId = ccreqIds[i];
        plint d = getReductionDimension(ccrId);
        if (d==1) { quantities1D[ccrId] = cqh.get1D(ccrId); }
        else if (d==3) { quantities3D[ccrId] = cqh.get3D(ccrId); }
        else { quantitiesND[ccrId] = cqh.getND(ccrId); }
    }
    this->make_ccrId_List();
}



template<typename T>
bool CellQuantityHolder<T>::has_ccrId(plint ccrId) {
	plint dim = getReductionDimension(ccrId);
	plint ret;
	if      (dim==1)  { ret = quantities1D.count(ccrId); }
	else if (dim==3)  { ret = quantities3D.count(ccrId); }
	else		      { ret = quantitiesND.count(ccrId); }
	return (ret>0);
}


template<typename T>
void CellQuantityHolder<T>::copyFromBlockStatisticsCCR(BlockStatisticsCCR<T> & reducer) {
    vector<plint> & ccrIds = reducer.get_ccrIds();
    for (pluint i = 0; i < ccrIds.size(); ++i) {
        plint ccrId = ccrIds[i];
        plint dim = getReductionDimension(ccrId);
        if (1==dim) { reducer.get(ccrId, get1D(ccrId)); }
        else if (3==dim) { reducer.get(ccrId, get3D(ccrId)); }
        else { reducer.get(ccrId, getND(ccrId)); }
    }
}


template<typename T>
void CellQuantityHolder<T>::reduceQuantity(plint ccrId, T value, plint numParts) {
    if (quantities1D.count(ccrId) == 0) { quantities1D[ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        T prValue = quantities1D[ccrId];
        if (0 == reductionType)      { quantities1D[ccrId] += value; } // Sum
        else if (1 == reductionType) {
            quantities1D[ccrId] = (prValue*numParts_ccrId[ccrId] + value*numParts) * 1.0 / (numParts_ccrId[ccrId] + numParts);
            numParts_ccrId[ccrId] += numParts;
        }  // Mean
        else if (2 == reductionType) { quantities1D[ccrId] = max(prValue, value); } // Max of
        else if (3 == reductionType) { quantities1D[ccrId] = min(prValue, value);  } // Min
        // STD, not implemented
        // else if (4 == reductionType) { false; } // Std not implemented yet
    }
}


template<typename T>
void CellQuantityHolder<T>::reduceQuantity(plint ccrId, Array<T,3> const& value, plint numParts) {
    plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
    if (quantities3D.count(ccrId) == 0) { quantities3D[ccrId] = value; }
    else {
        Array<T,3> prValue = quantities3D[ccrId];
        if (0 == reductionType)      { quantities3D[ccrId] = prValue + value; }
        else if (1 == reductionType) {
            plint prNumParts =  numParts_ccrId[ccrId];
            quantities3D[ccrId] = prNumParts * 1.0 * prValue + numParts * 1.0 * value;
            T factr = 1.0 / (prNumParts + numParts);
            quantities3D[ccrId] = quantities3D[ccrId] * factr;
            numParts_ccrId[ccrId] += numParts;
        }
        else if (2 == reductionType) {
            quantities3D[ccrId] = Array<T,3>( max(prValue[0], value[0]),
                                                      max(prValue[1], value[1]),
                                                      max(prValue[2], value[2]) );
        }
        else if (3 == reductionType) {
            quantities3D[ccrId] = Array<T,3>( min(prValue[0], value[0]),
                                                      min(prValue[1], value[1]),
                                                      min(prValue[2], value[2]) );
        }
    }
}


template<typename T>
void CellQuantityHolder<T>::reduceQuantity(plint ccrId, std::vector<T> const& value, plint numParts) {
    if (quantitiesND.count(ccrId) == 0) { quantitiesND[ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        for (pluint iv = 0; iv < value.size(); ++iv) {
            T prValue = quantitiesND[ccrId][iv];
            if (0 == reductionType)      { quantitiesND[ccrId][iv] = prValue + value[iv]; }
            else if (1 == reductionType) {
                quantitiesND[ccrId][iv] = (prValue*numParts_ccrId[ccrId] + numParts*value[iv]) * 1.0 / (numParts_ccrId[ccrId] + numParts);
            }
            else if (2 == reductionType) { quantitiesND[ccrId][iv] = max(prValue, value[iv]); }
            else if (3 == reductionType) { quantitiesND[ccrId][iv] = min(prValue, value[iv]);  } // Min can be calculated from the inverse Max
    //            else if (4 == reductionType) { false; } // Std not implemented yet
        }
        if (1 == reductionType) {
            numParts_ccrId[ccrId] += numParts;
        }
    }
}


template<typename T>
void CellQuantityHolder<T>::make_ccrId_List() {
    if (scalar_ccrIds.size() != quantities1D.size()) {
        scalar_ccrIds.clear();
        typename std::map<plint, T >::const_iterator iter1D;
        for (iter1D  = quantities1D.begin(); iter1D != quantities1D.end(); ++iter1D) { scalar_ccrIds.push_back(iter1D->first); }
        std::sort(scalar_ccrIds.begin(), scalar_ccrIds.end());
    }
    if (vector_ccrIds.size() != quantities3D.size()) {
        vector_ccrIds.clear();
        typename std::map<plint, Array<T,3> >::const_iterator iter3D;
        for (iter3D  = quantities3D.begin(); iter3D != quantities3D.end(); ++iter3D) { vector_ccrIds.push_back(iter3D->first); }
        std::sort(vector_ccrIds.begin(), vector_ccrIds.end());
    }
    if (tensor_ccrIds.size() != quantitiesND.size()) {
        tensor_ccrIds.clear();
        typename std::map<plint, std::vector<T> >::const_iterator iterND;
        for (iterND  = quantitiesND.begin(); iterND != quantitiesND.end(); ++iterND) { tensor_ccrIds.push_back(iterND->first); }
        std::sort(tensor_ccrIds.begin(), tensor_ccrIds.end());
    }
}




  /**************************************************/
 /********************* Cell3D *********************/
/**************************************************/


template<typename T, template<typename U> class Descriptor>
Cell3D<T, Descriptor>::Cell3D(TriangularSurfaceMesh<T> & mesh_, plint cellId_) :
    CellQuantityHolder<T>(), mesh(mesh_), cellId(cellId_) {
    setMesh();
    for (plint i = 0; i < allReductions.size(); i++) {
        plint d = getReductionDimension(allReductions[i]);
        if (d == 1) {
            this->quantities1D[allReductions[i]] = *(new T()); 
        }
    }
};

template<typename T, template<typename U> class Descriptor>
Cell3D<T, Descriptor>::Cell3D(Cell3D<T,Descriptor> const& rhs) :
    CellQuantityHolder<T>(rhs), mesh(rhs.mesh), cellId(rhs.cellId) {
    setMesh();
    for (plint i = 0; i < allReductions.size(); i++) {
        plint d = getReductionDimension(allReductions[i]);
        if (d == 1) {
            this->quantities1D[allReductions[i]] = *(new T()); 
        }
    }
};

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::setMesh() {
	cellNumVertices	 = mesh.getNumVertices();
	cellNumTriangles = mesh.getNumTriangles();
};

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::saveMesh(std::string fileName) {
	if ( vertices.size() == mesh.getNumVertices() ) {
		ElementsOfTriangularSurfaceMesh<T> emptyEoTSM;
		TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(mesh, emptyEoTSM);
	    for (int iV = 0; iV< meshCopy->getNumVertices(); ++iV) {
	    	meshCopy->replaceVertex(iV,  getVertex(iV));
	    }
		writeSurfaceMeshAsciiSTL(*meshCopy, fileName);
	}
};


template<typename T, template<typename U> class Descriptor>
plint Cell3D<T, Descriptor>::getMpiProcessor() {
 #ifdef PLB_MPI_PARALLEL
    //  Get the individual process ID.
    int rank = MPI::COMM_WORLD.Get_rank();
#else
    int rank = 0;
#endif
    return rank;
}

template<typename T, template<typename U> class Descriptor>
plint Cell3D<T,Descriptor>::getEdgeId(plint iVertex, plint jVertex) {
    iVertex = iVertex % cellNumVertices;
    jVertex = jVertex % cellNumVertices;
    if (iVertex > jVertex){
        return (iVertex*(iVertex - 1))/2 + jVertex;
    } else if (iVertex < jVertex) {
        return (jVertex*(jVertex - 1))/2 + iVertex;
    }
    return -1;
};

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::push_back(Particle3D<T,Descriptor>* particle3D, bool particleIsInBulk) {
    plint vertexId = castParticleToICP3D(particle3D)->getVertexId();
	if (particleIsInBulk) {
	    iVertexToParticle3D[vertexId] = particle3D;
	    verticesInBulk.insert(vertexId);
	} else {
	    // For some reason Palabos deletes particles that are in the envelope after some functionals;
	    iVertexToParticle3D[vertexId] = particle3D->clone();
	    particlesToDelete.push_back(iVertexToParticle3D[vertexId]);
	}
}

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::close() {
    plint numTrianges = mesh.getNumTriangles();
    triangles.clear();
    vertices.clear();
    edges.clear();
    triangleAreas.clear();
    edgeLengths.clear();
    edgeLengthVectors.clear();
    signedAngles.clear();
    triangleAreas.clear();
    triangleNormals.clear();
    vertexAreas.clear();
    vertexNormals.clear();
    edgeTileSpans.clear();

    std::map<plint, Array<plint,2> > edgeMap;
    std::set<plint> triangleSet;
    for (int iTriangle = 0; iTriangle < numTrianges; ++iTriangle)
    {
        plint vId0 = mesh.getVertexId(iTriangle, 0);
        plint vId1 = mesh.getVertexId(iTriangle, 1);
        plint vId2 = mesh.getVertexId(iTriangle, 2);
        plint v0inBulk = verticesInBulk.count(vId0); // verticesInBulk
        plint v1inBulk = verticesInBulk.count(vId1);
        plint v2inBulk = verticesInBulk.count(vId2);
        plint numVert = v0inBulk + v1inBulk + v2inBulk;
        if (numVert > 0) {
            triangleSet.insert(iTriangle);
            edgeMap[getEdgeId(vId0, vId1)] = Array<plint, 2>(vId0, vId1) ;
            edgeMap[getEdgeId(vId0, vId2)] = Array<plint, 2>(vId0, vId2) ;
            edgeMap[getEdgeId(vId1, vId2)] = Array<plint, 2>(vId1, vId2) ;
            std::vector<plint> adjacentTriangleIds = mesh.getAdjacentTriangleIds(iTriangle);
            // Bending force requires also the adjacent triangles of the edge.
            for (pluint iatId = 0; iatId < adjacentTriangleIds.size(); ++iatId) {
                plint jT = adjacentTriangleIds[iatId];
                triangleSet.insert(jT);
                vId0 = mesh.getVertexId(iTriangle, 0);
                vId1 = mesh.getVertexId(iTriangle, 1);
                vId2 = mesh.getVertexId(iTriangle, 2);
                edgeMap[getEdgeId(vId0, vId1)] = Array<plint, 2>(vId0, vId1) ;
                edgeMap[getEdgeId(vId0, vId2)] = Array<plint, 2>(vId0, vId2) ;
                edgeMap[getEdgeId(vId1, vId2)] = Array<plint, 2>(vId1, vId2) ;
            }
        }
    }

    typename std::map<plint, Array<plint,2> >::iterator iter;
    for (iter = edgeMap.begin(); iter != edgeMap.end(); ++iter)
    {
        edges.push_back(iter->second);
    }

    typename std::set<plint>::iterator iterSet;
    for (iterSet = triangleSet.begin(); iterSet != triangleSet.end(); ++iterSet)
    {
        triangles.push_back(*iterSet);
    }

    typename std::map<plint, Particle3D<T,Descriptor>* >::iterator it;
    for (it  = iVertexToParticle3D.begin(); it != iVertexToParticle3D.end(); ++it) {
        vertices.push_back(it->first);
    }
    this->getParticlesPerCellId() = getNumVertices_Local();
}

  /*________________________________________________*/
 /***** Cell3D -  Calcuate Quantities **************/
/*________________________________________________*/


template<typename T, template<typename U> class Descriptor>
Array<T,3> Cell3D<T, Descriptor>::computeEdgeLengthVector(plint iVertex, plint jVertex) {
	plint edgeId = getEdgeId(iVertex, jVertex);
	if (edgeLengthVectors.count(edgeId) == 0) {
		if (iVertex > jVertex) { edgeLengthVectors[edgeId] = getVertex(iVertex) - getVertex(jVertex) ; }
		else 				   { edgeLengthVectors[edgeId] = getVertex(jVertex) - getVertex(iVertex) ; }
		edgeLengths[edgeId] = norm(edgeLengthVectors[edgeId]);
	}
	if (iVertex > jVertex) { return edgeLengthVectors[edgeId]; }
	else 				   { return edgeLengthVectors[edgeId] * -1.0; }
}

template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeEdgeLength(plint iVertex, plint jVertex) {
	plint edgeId = getEdgeId(iVertex, jVertex);
	if (edgeLengths.count(edgeId) == 0) {
		computeEdgeLengthVector(iVertex, jVertex);
	}
	return edgeLengths[edgeId];
}

template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeSignedAngle(plint iVertex, plint jVertex, bool& found) {
    plint edgeId = getEdgeId(iVertex, jVertex);
    found = true;
    if (signedAngles.count(edgeId) == 0) {
        plint kVertex, lVertex;
        T angle = computeSignedAngle(iVertex, jVertex,kVertex, lVertex,found);
        if (found) { signedAngles[edgeId] = angle;  }
        else       { signedAngles.erase(edgeId); return angle;}
    }
    return signedAngles[edgeId];
}

template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeSignedAngle(plint iVertex, plint jVertex, plint & kVertex, plint & lVertex, bool& found) {
	plint edgeId = getEdgeId(iVertex, jVertex);
	found = true;
    Array<T,3> x1 = getVertex(iVertex), x2(0.,0.,0.), x3(0.,0.,0.), x4(0.,0.,0.);

    std::vector<plint> adjacentTriangles = getAdjacentTriangleIds(iVertex, jVertex);
    plint iTriangle=adjacentTriangles[0], jTriangle=adjacentTriangles[1];
    x3 = getVertex(jVertex);
    T foundVertices=0;
    for (pluint id = 0; id < 3; ++id) {
        kVertex = getVertexId(iTriangle,id);
        if ( (kVertex != iVertex) && (kVertex != jVertex) && (iVertexToParticle3D.count(kVertex)>0)) {
            x2 = getVertex(kVertex);
            foundVertices += 1;
            break;
        }
    }
    for (pluint id = 0; id < 3; ++id) {
        lVertex = getVertexId(jTriangle,id);
        if ( (lVertex != iVertex) && (lVertex != jVertex) && (iVertexToParticle3D.count(lVertex)>0)) {
            x4 = getVertex(lVertex);
            foundVertices += 1;
            break;
        }
    }
    found = (foundVertices == 2); //Assert if some particles are outside of the domain
    if (not found) { signedAngles.erase(edgeId);  return 0.0; }
	if (signedAngles.count(edgeId) == 0) {
	    Array<T,3> V1 = computeTriangleNormal(iTriangle);
	    Array<T,3> V2 = computeTriangleNormal(jTriangle);
	    T angle = angleBetweenVectors(V1, V2);
		plint sign = dot(x2-x1, V2) >= 0?1:-1;

		if (sign <= 0) {
			angle = 2*pi-angle;
		}
		angle = (angle > pi)?angle-2*pi:angle;
		signedAngles[edgeId] = angle;
	}
	return signedAngles[edgeId];
}

template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeTriangleArea(plint iTriangle) {
    PLB_ASSERT(iTriangle >= 0 && iTriangle < plint(cellNumTriangles));
	if (triangleAreas.count(iTriangle) == 0) {
	    Array<T,3> v0 = getVertex(iTriangle, 0);
	    Array<T,3> v1 = getVertex(iTriangle, 1);
	    Array<T,3> v2 = getVertex(iTriangle, 2);

	    Array<T,3> e01 = v1 - v0;
	    Array<T,3> e02 = v2 - v0;

	    Array<T,3> n;
	    crossProduct(e01, e02, n);
	    T normN = norm(n);
        n /= normN;
	    triangleAreas[iTriangle] = normN * 0.5;
	    triangleNormals[iTriangle] = n;
	}
    return triangleAreas[iTriangle];
}

inline void fic_swap(plint& a, plint& b)
{
    plint temp = a;
    a = b;
    b = temp;
}

template<typename T, template<typename U> class Descriptor>
plint Cell3D<T, Descriptor>::findTriangleId(plint iVertex, plint jVertex, plint kVertex) 
{ 
    #if 0
   
    std::vector<plint> ati1 = getAdjacentTriangleIds(iVertex, jVertex); 
    std::vector<plint> ati2 = getAdjacentTriangleIds(iVertex, kVertex); 
    if (ati1[0] == ati2[0]) return ati1[0];
    if (ati1[0] == ati2[1]) return ati1[0];
    if (ati1[1] == ati2[0]) return ati1[1];
    if (ati1[1] == ati2[1]) return ati1[1];
    return -1;
    #else

    plint v[] = {iVertex,jVertex,kVertex};

    if(v[0] > v[1]) fic_swap(v[0], v[1]);
    if(v[1] > v[2]) fic_swap(v[1], v[2]);
    if(v[0] > v[1]) fic_swap(v[0], v[1]);

    plint key = (v[0]<<20) | (v[1]<<10) | v[2];

    if(triangleIDs.count(key) == 0)
    {   
        std::vector<plint> ati1 = getAdjacentTriangleIds(iVertex, jVertex); 
        std::vector<plint> ati2 = getAdjacentTriangleIds(iVertex, kVertex); 
        if (ati1[0] == ati2[0]) triangleIDs[key] = ati1[0];
        else if (ati1[0] == ati2[1]) triangleIDs[key] = ati1[0];
        else if (ati1[1] == ati2[0]) triangleIDs[key] = ati1[1];
        else if (ati1[1] == ati2[1]) triangleIDs[key] = ati1[1];
        else triangleIDs[key] = -1;
    }
    return triangleIDs[key];

    #endif
} ;


template<typename T, template<typename U> class Descriptor>
Array<T,3> Cell3D<T, Descriptor>::computeTriangleNormal(plint iTriangle) {
    PLB_ASSERT(iTriangle >= 0 && iTriangle < plint(cellNumTriangles));
    if (triangleAreas.count(iTriangle) == 0) {
        computeTriangleArea(iTriangle);
    }
    return triangleNormals[iTriangle];
}



template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeVertexArea(plint iVertex) {
	if (vertexAreas.count(iVertex) == 0) {
	    std::vector<plint> neighborTriangleIds = getNeighborTriangleIds(iVertex);
	    Array<T,3> n(0.,0.,0.);
	    T s = 0.0;
	    std::vector<plint>::iterator it = neighborTriangleIds.begin();
	    for (n.resetToZero(); it != neighborTriangleIds.end(); ++it) {
	        s += computeTriangleArea(*it) ;
	    	n += computeTriangleNormal(*it) * computeTriangleArea(*it) ;
	    }
	    T normN = norm(n);
	    n /= normN;
	    vertexAreas[iVertex] = s/3.0;
	    vertexNormals[iVertex] = n;
	}
    return vertexAreas[iVertex];
}

template<typename T, template<typename U> class Descriptor>
Array<T,3>  Cell3D<T, Descriptor>::computeVertexNormal(plint iVertex) {
	if (vertexNormals.count(iVertex) == 0) {
	    computeVertexArea(iVertex);
	}
    return vertexNormals[iVertex];
}



template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeEdgeTileSpan(plint iVertex, plint jVertex) {
	plint edgeId = getEdgeId(iVertex, jVertex);
	if (edgeTileSpans.count(edgeId) == 0) {
	    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex, jVertex);

	    Array<T,3> v0 = getVertex(neighborVertexIds[0]);
	    Array<T,3> v1 = getVertex(iVertex);
	    Array<T,3> v2 = getVertex(jVertex);

	    Array<T,3> v01 = v1 - v0;
	    Array<T,3> v21 = v1 - v2;
	    T angle_012 = angleBetweenVectors(v21, v01);

	    T span = fabs(sin(angle_012)) * norm(v01);

	    if (neighborVertexIds.size() == 1)
	        return span/6.0;

	    Array<T,3> v3 = getVertex(neighborVertexIds[1]);

	    Array<T,3> v32 = v2 - v3;
	    Array<T,3> v12 = -v21;
	    T angle_321 = angleBetweenVectors(v12, v32);
	    
	    span += fabs(sin(angle_321)) * norm(v32);
	    edgeTileSpans[edgeId] = span/6.0;

	}
    return edgeTileSpans[edgeId];
}





template<typename T, template<typename U> class Descriptor>
void calculateCCRQuantities(plint ccrId, BlockStatisticsCCR<T> & reducer, Cell3D<T, Descriptor> * cell, plint iVertex) {
    std::vector<plint> neighbors = cell->getNeighborVertexIds(iVertex);
    plint q = getReductionQuantity(ccrId);
/****** 1D Quantities ******/
    // Calculate ANGLE
    if (q==2) { 
        T edgeAngle = 0.0; bool found;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeAngle += cell->computeSignedAngle(iVertex, neighbors[iB], found);
        }
        reducer.gather(ccrId, edgeAngle*1.0/neighbors.size() );
    // Calculate AREA
    } else if (q==3) {
//        cout << iVertex << "cell->computeVertexArea(iVertex) " << cell->computeVertexArea(iVertex) << std::endl;
        reducer.gather(ccrId, cell->computeVertexArea(iVertex) );
    // Calculate EDGE DISTANCE
    } else if (q==4) { 
        T edgeDistance = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) 
        { edgeDistance +=  cell->computeEdgeLength(iVertex, neighbors[iB]) ; }
        reducer.gather(ccrId, edgeDistance*1.0/neighbors.size());
    // Calculate EDGE TILE SPAN
    } else if (q==5) { 
        T edgeTileSpan = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) 
        { edgeTileSpan +=  cell->computeEdgeTileSpan(iVertex, neighbors[iB]) ; }
        reducer.gather(ccrId, edgeTileSpan = edgeTileSpan*1.0/neighbors.size() );
    // Return Energy of particle
#ifdef PLB_DEBUG // Less Calculations
    } else if (q==9) { reducer.gather(ccrId, cell->get_Energy(iVertex));
#endif
    // Calculate VOLUME
    } else if (q==1) { 
        std::vector<plint> neighborTriangleIds = cell->getNeighborTriangleIds(iVertex);
        T quantity1D = 0.0;
        for (pluint iB = 0; iB < neighborTriangleIds.size(); ++iB) {
            plint iTriangle = neighborTriangleIds[iB];
            Array<T,3> v0 = cell->getVertex(iTriangle, 0);
            Array<T,3> v1 = cell->getVertex(iTriangle, 1);
            Array<T,3> v2 = cell->getVertex(iTriangle, 2);
            Array<T,3> tmp;
            crossProduct(v1, v2, tmp);
            T triangleVolumeT6 =  VectorTemplate<T,Descriptor>::scalarProduct(v0,tmp); // * (1.0/6.0)
            quantity1D += triangleVolumeT6/6.0/3.0; // every volume is evaluated 3 times
        }
        reducer.gather(ccrId, quantity1D);
    // Vertex distance
    } else if (q==18) { reducer.gather(ccrId, norm( cell->getPosition() - cell->getPosition(iVertex)) ); // POSITION
/****** 3D Quantities ******/
    // CCR_NO_PBC_POSITION_MEAN
    } else if (q==0) { reducer.gather(ccrId, cell->getPosition(iVertex)); // POSITION
    // POSITION FROM PERIODIC BOUNDARY CONDITION
    } else if (q==6) { reducer.gather(ccrId, cell->get_pbcPosition(iVertex));
    // VELOCITY
    } else if (q==7) { reducer.gather(ccrId, cell->get_v(iVertex)) ;
    // Force
    } else if (q==16) {
        reducer.gather(ccrId, cell->get_force(iVertex));
	// Normalized force
	} else if (q==26) {
		// Normalize force with area. Force = force * Vi/Vmean;
		Array<T,3> normForce = cell->get_force(iVertex) * cell->computeVertexArea(iVertex) * 1.0/ (cell->getSurface()*1.0/cell->getCellNumVertices()) ;
		reducer.gather(ccrId, normForce );
/****** ND Quantities ******/
    // INERTIA
    } else if (q==8) {
        std::vector<T> quantityND;
        T rx, ry, rz;
        T Ixx=0, Ixy=0, Ixz=0;
        T Iyx=0, Iyy=0, Iyz=0;
        T Izx=0, Izy=0, Izz=0;
        Array<T,3> r0 = cell->getPosition(); // Get the Cell Center;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            T Aj = cell->computeTriangleArea(neighbors[iB]);
            Array<T,3> nj = cell->computeTriangleNormal(neighbors[iB]);
            Array<T,3> rj = (cell->getVertex(neighbors[iB],0) +
                            cell->getVertex(neighbors[iB],1) + cell->getVertex(neighbors[iB],2))/3.0;
            rj = rj - r0;
            T dV = 1.0/5.0 * Aj * dot(nj, rj);
            rx = rj[0]; ry = rj[1]; rz = rj[2];
            Ixx += (rz*rz+ry*ry)*dV;
            Iyy += (rz*rz+rx*rx)*dV;
            Izz += (rx*rx+ry*ry)*dV;
            Iyx += -rx*ry*dV;
            Ixy += -rx*ry*dV;
            Izx += -rx*rz*dV;
            Ixz += -rx*rz*dV;
            Izy += -ry*rz*dV;
            Iyz += -ry*rz*dV;
        }
        quantityND.clear();
        quantityND.push_back(Ixx/3.0); // [0] every element is evaluated 3 times
        quantityND.push_back(Ixy/3.0); // [1] every element is evaluated 3 times
        quantityND.push_back(Ixz/3.0); // [2] every element is evaluated 3 times
        quantityND.push_back(Iyx/3.0); // [3] every element is evaluated 3 times
        quantityND.push_back(Iyy/3.0); // [4] every element is evaluated 3 times
        quantityND.push_back(Iyz/3.0); // [5] every element is evaluated 3 times
        quantityND.push_back(Izx/3.0); // [6] every element is evaluated 3 times
        quantityND.push_back(Izy/3.0); // [7] every element is evaluated 3 times
        quantityND.push_back(Izz/3.0); // [8] every element is evaluated 3 times
        reducer.gather(ccrId, quantityND);
    }

}




#endif  // CELL_3D_HH

