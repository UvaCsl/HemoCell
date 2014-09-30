#ifndef CELL_3D_HH
#define CELL_3D_HH
#include "cell3D.h"

/********************* CellQuantityHolder *********************/

template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::clearQuantities() {
    quantities1D.clear();
    quantities3D.clear();
    quantitiesND.clear();
    scalar_ccrIds.clear(); vector_ccrIds.clear(); tensor_ccrIds.clear();
}

template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::updateCQH(CellQuantityHolder<T> * cqh) {
    this->clearQuantities();
    particlesPerCellId = cqh->getParticlesPerCellId() ;
    quantities1D = cqh->getQuantities1D();
    quantities3D = cqh->getQuantities3D();
    quantitiesND = cqh->getQuantitiesND();
    this->make_ccrId_List();
}

template<typename T, template<typename U> class Descriptor>
bool CellQuantityHolder<T, Descriptor>::count(plint ccrId) {
	plint dim = getReductionDimension(ccrId);
	if      (dim==1)  { ret = quantities1D.count(ccrId); }
	else if (dim==3)  { ret = quantities3D.count(ccrId); }
	else		      { ret = quantitiesND.count(ccrId); }
	return (ret==1)
}


template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::reduceQuantity(plint ccrId, T value, plint numParts) {
    if (quantities1D.count(ccrId) == 0) { quantities1D[ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        T prValue = quantities1D[ccrId];
        if (0 == reductionType)      { quantities1D[ccrId] += value; } // Sum
        else if (1 == reductionType) { quantities1D[ccrId] = 
                (prValue*particlesPerCellId + value*numParts) * 1.0 / (particlesPerCellId + numParts); }  // Mean
        else if (2 == reductionType) { quantities1D[ccrId] = max(prValue, value); } // Max of
        else if (3 == reductionType) { quantities1D[ccrId] = min(prValue, value);  } // Min
        // STD, not implemented
        // else if (4 == reductionType) { false; } // Std not implemented yet
    }
}

template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::reduceQuantity(plint ccrId, T value, plint numParts) {
    if (quantities1D.count(ccrId) == 0) { quantities1D[ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        T prValue = quantities1D[ccrId];
        if (0 == reductionType)      { quantities1D[ccrId] += value; } // Sum
        else if (1 == reductionType) { quantities1D[ccrId] = 
                (prValue*particlesPerCellId + value*numParts) * 1.0 / (particlesPerCellId + numParts); }  // Mean
        else if (2 == reductionType) { quantities1D[ccrId] = max(prValue, value); } // Max of
        else if (3 == reductionType) { quantities1D[ccrId] = min(prValue, value);  } // Min
        // STD, not implemented
        // else if (4 == reductionType) { false; } // Std not implemented yet
    }
}

template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::reduceQuantity(plint ccrId, std::vector<T> const& value, plint numParts) {
    if (quantitiesND.count(ccrId) == 0) { quantitiesND[ccrId] = value; }
    else {
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        for (pluint iv = 0; iv < value.size(); ++iv) {
            T prValue = quantitiesND[ccrId][iv];
            if (0 == reductionType)      { quantitiesND[ccrId][iv] = prValue + value[iv]; }
            else if (1 == reductionType) { quantitiesND[ccrId][iv] = (prValue*particlesPerCellId + numParts*value[iv]) * 1.0 / (particlesPerCellId + numParts); }
            else if (2 == reductionType) { quantitiesND[ccrId][iv] = max(prValue, value[iv]); }
            else if (3 == reductionType) { quantitiesND[ccrId][iv] = min(prValue, value[iv]);  } // Min can be calculated from the inverse Max
    //            else if (4 == reductionType) { false; } // Std not implemented yet
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void CellQuantityHolder<T, Descriptor>::make_ccrId_List() {
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
Cell3D<T, Descriptor>::Cell3D(TriangularSurfaceMesh<T>& mesh_, plint cellId_=-1) {
    cellId = cellId_;
    setMesh(mesh_);
};

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::setMesh(TriangularSurfaceMesh<T>& mesh_) {
	mesh = mesh_;
	cellNumVertices	 = mesh.getNumVertices();
	cellNumTriangles = mesh.getNumTriangles();
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
void Cell3D<T, Descriptor>::push_back(Particle3D<T,Descriptor>* particle3D) {  
	iVertexToParticle3D[castParticleToICP3D(iVertexToParticle3D[iVertex])->getVertexId()] = particle3D; 
}

template<typename T, template<typename U> class Descriptor>
void Cell3D<T, Descriptor>::close() {
    plint numTrianges = mesh.getNumTriangles();
    triangles.clear();
    for (int iTriangle = 0; iTriangle < numTrianges; ++iTriangle)
    {
        vId0 = getVertexId(iTriangle, 0);
        vId1 = getVertexId(iTriangle, 1);
        vId2 = getVertexId(iTriangle, 2);
        numVert = iVertexToParticle3D.count(vId0) 
                + iVertexToParticle3D.count(vId1)
                + iVertexToParticle3D.count(vId2);
        if (numVert == 3) {
            triangles.push_back(iTriangle);
        }
    }
    vertices.clear();
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
T Cell3D<T, Descriptor>::computeSignedAngle(plint iVertex, plint jVertex) {
	plint edgeId = getEdgeId(iVertex, jVertex);
	if (signedAngles.count(edgeId) == 0) {
		plint kVertex, lVertex;
	    Array<T,3> x1 = getVertex(iVertex), x2(0.,0.,0.), x3(0.,0.,0.), x4(0.,0.,0.);

	    std::vector<plint> adjacentTriangles = getAdjacentTriangleIds(iVertex, jVertex);
		plint iTriangle=adjacentTriangles[0], jTriangle=adjacentTriangles[1];
	    x3 = getVertex(jVertex);
	    T foundVertices=0;
	    for (pluint id = 0; id < 3; ++id) {
	        kVertex = getVertexId(iTriangle,id);
	        if ( (kVertex != iVertex) && (kVertex != jVertex) ) {
	            x2 = getVertex(kVertex);
	            foundVertices += 1;
	            break;
	        }
	    }
	    for (pluint id = 0; id < 3; ++id) {
	        lVertex = getVertexId(jTriangle,id);
	        if ( (lVertex != iVertex) && (lVertex != jVertex) ) {
	            x4 = getVertex(lVertex);
	            foundVertices += 1;
	            break;
	        }
	    }
	    PLB_ASSERT(foundVertices == 2); //Assert if some particles are outside of the domain

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
    PLB_ASSERT(iTriangle >= 0 && iTriangle < cellNumTriangles);
	if (triangleAreas.count(iTriangle) == 0) {
	    Array<T,3> v0 = getVertex(iTriangle, 0);
	    Array<T,3> v1 = getVertex(iTriangle, 1);
	    Array<T,3> v2 = getVertex(iTriangle, 2);

	    Array<T,3> e01 = v1 - v0;
	    Array<T,3> e02 = v2 - v0;

	    Array<T,3> n;
	    crossProduct(e01, e02, n);
	    normN = normN;
        n /= normN;
	    triangleAreas[iTriangle] = normN;
	    triangleNormals[iTriangle] = n;
	}
    return triangleAreas[iTriangle];
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> Cell3D<T, Descriptor>::computeTriangleNormal(plint iTriangle) {
    PLB_ASSERT(iTriangle >= 0 && iTriangle < cellNumTriangles);
	if (triangleAreas.count(iTriangle) == 0) {
		computeTriangleArea(iTriangle);
	}
	return triangleNormals[iTriangle];
}

template<typename T, template<typename U> class Descriptor>
T Cell3D<T, Descriptor>::computeVertexArea(plint iVertex) {
	if (vertexAreas.count(iVertex) == 0) {
	    std::vector<plint> neighborTriangleIds = getNeighborTriangleIds(iVertex);
	    Array<T,3> n;
	    std::vector<plint>::iterator it = neighborTriangleIds.begin();
	    for (n.resetToZero(); it != neighborTriangleIds.end(); ++it)
	        n += computeTriangleNormal(*it);
	    T normN = norm(n);
	    n /= normN;
	    vertexAreas[iVertex] = normN;
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




#endif  // CELL_3D_HH
