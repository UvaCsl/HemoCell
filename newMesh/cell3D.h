#ifndef CELL_3D_H
#define CELL_3D_H

/*
 * Heavily depends on triangularSurfaceMesh
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include "cellReductionTypes.h"
#include <vector>
#include <map>

using namespace std;
using namespace plb;



template<typename T>
class CellQuantityHolder
{
public:
    CellQuantityHolder() { } ;
    ~CellQuantityHolder() { };

    void setParticlesPerCellId(plint particlesPerCellId_) { particlesPerCellId = particlesPerCellId_; }

    void clearQuantities() {
        quantities1D.clear();
        quantities3D.clear();
        quantitiesND.clear();
        scalar_ccrIds.clear(); vector_ccrIds.clear(); tensor_ccrIds.clear();
    }


    void insert(plint ccrId, T value) { quantities1D[ccrId] = value; }
    void insert(plint ccrId, Array<T,3> value) { quantities3D[ccrId] = value; }
    void insert(plint ccrId, std::vector<T> value) { quantitiesND[ccrId] = value; }

    std::map<plint, T >& getQuantities1D() { return quantities1D; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> >& getQuantities3D() { return quantities3D; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, std::vector<T> >& getQuantitiesND() { return quantitiesND; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME

    T const& get1D(plint ccrId) { return quantities1D[ccrId]; }
    Array<T,3> const& get3D(plint ccrId) { return quantities3D[ccrId]; }
    std::vector<T> const& getND(plint ccrId) { return quantitiesND[ccrId]; }
    T & get1D(plint ccrId) { return quantities1D[ccrId]; }
    Array<T,3> & get3D(plint ccrId) { return quantities3D[ccrId]; }
    std::vector<T> & getND(plint ccrId) { return quantitiesND[ccrId]; }
    void get(plint ccrId, T& value) { value = get1D(ccrId); return value }
    void get(plint ccrId, Array<T,3>& value) { value = get3D(ccrId); return value }
    void get(plint ccrId, std::vector<T>& value) { value = getND(ccrId); return value }


    void reduceQuantity(plint ccrId, T value, plint numParts) {
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

    void reduceQuantity(plint ccrId, Array<T,3> const& value, plint numParts) {
        if (quantities3D.count(ccrId) == 0) { quantities3D[ccrId] = value; }
        else {
            plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
            Array<T,3> prValue = quantities3D[ccrId];
            if (0 == reductionType)      { quantities3D[ccrId] = prValue + value; }
            else if (1 == reductionType) {
                plint prNumParts =  particlesPerCellId;
                quantities3D[ccrId] = prNumParts * 1.0 * prValue + numParts * 1.0 * value;
                T factr = 1.0 / (prNumParts + numParts);
                quantities3D[ccrId] = quantities3D[ccrId] * factr;
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

    void reduceQuantity(plint ccrId, std::vector<T> const& value, plint numParts) {
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


    void make_ccrId_List() {
        plint containedData = scalar_ccrIds.size() + vector_ccrIds.size() + tensor_ccrIds.size();
        if (containedData == 0) {
            typename std::map<plint, T >::const_iterator iter1D;
            typename std::map<plint, Array<T,3> >::const_iterator iter3D;
            typename std::map<plint, std::vector<T> >::const_iterator iterND;
            scalar_ccrIds.clear();     vector_ccrIds.clear();     tensor_ccrIds.clear();
            for (iter1D  = quantities1D.begin(); iter1D != quantities1D.end(); ++iter1D) { scalar_ccrIds.push_back(iter1D->first); }
            for (iter3D  = quantities3D.begin(); iter3D != quantities3D.end(); ++iter3D) { vector_ccrIds.push_back(iter3D->first); }
            for (iterND  = quantitiesND.begin(); iterND != quantitiesND.end(); ++iterND) { tensor_ccrIds.push_back(iterND->first); }
            std::sort(scalar_ccrIds.begin(), scalar_ccrIds.end());
            std::sort(vector_ccrIds.begin(), vector_ccrIds.end());
            std::sort(tensor_ccrIds.begin(), tensor_ccrIds.end());
        }
    }

    std::vector<plint> const& getScalarCcrIds() { make_ccrId_List(); return scalar_ccrIds; }
    std::vector<plint> const& getVectorCcrIds() { make_ccrId_List(); return vector_ccrIds; }
    std::vector<plint> const& getTensorCcrIds() { make_ccrId_List(); return tensor_ccrIds; }
private:
    plint cellId, particlesPerCellId;

    std::map<plint, T > quantities1D; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> > quantities3D;
    std::map<plint, std::vector<T> > quantitiesND;
    std::vector<plint> scalar_ccrIds, vector_ccrIds, tensor_ccrIds;



public:
    T & getVolume() { return quantities1D[CCR_VOLUME];} ;
    T & getSurface() { return quantities1D[CCR_SURFACE];} ;
    T & getEnergy() { return quantities1D[CCR_ENERGY];} ;
    T & getMeanAngle() { return quantities1D[CCR_ANGLE_MEAN];} ;
    T & getMeanEdgeLength() { return quantities1D[CCR_EDGE_DISTANCE_MEAN];} ;
    T & getMeanEdgeDistance() { return quantities1D[CCR_EDGE_DISTANCE_MEAN];} ;
    T & getMaxEdgeLength() { return quantities1D[CCR_EDGE_DISTANCE_MAX];} ;
    T & getMaxEdgeDistance() { return quantities1D[CCR_EDGE_DISTANCE_MAX];} ;
    T getMeanTriangleArea(plint numTriangles) { return quantities1D[CCR_SURFACE]*1.0/numTriangles;} ;
    T & getMeanTileSpan() { return quantities1D[CCR_TILE_SPAN_MEAN];} ;

    Array<T,3> const& getPosition() { return quantities3D[CCR_POSITION_MEAN];  } ;
    Array<T,3> const& getVelocity() { return quantities3D[CCR_VELOCITY_MEAN];} ;
    Array<T,3> const& getForce() { return quantities3D[CCR_FORCE];} ;
    std::vector<T> & getInertia() { return quantitiesND[CCR_INERTIA];} ;

    Array<T,3> & getTumblingAngles() { return quantities3D[CCR_TUMBLING_ANGLES];  } ;
    Array<T,3> & getTankTreadingAngles() { return quantities3D[CCR_TANK_TREADING_ANGLES];  } ;
    Array<T,3> & getDiameters()  { return quantities3D[CCR_DIAMETERS];  } ;
    T & getSymmetryDeviation()  { return quantities1D[CCR_SYMMETRY_DEVIATION];  } ;
    T & getDeformationIndex()  { return quantities1D[CCR_DEFORMATION_INDEX];  } ;
    T & getTaylorDeformationIndex()  { return quantities1D[CCR_TAYLOR_DEFORMATION_INDEX];  } ;

};




template<typename T, template<typename U> class Descriptor>
class Cell3D : public CellQuantityHolder<T>
{
public:
    Cell3D(TriangularSurfaceMesh<T>& mesh_, plint cellId_=-1) {
        cellId = cellId_;
        this->setParticlesPerCellId(getNumVertices_Global());
    };
    ~Cell3D() {};

    void push_back(Particle3D<T,Descriptor>* particle3D) { 
        ImmersedCellParticle3D<T,Descriptor>* particle = 
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particle3D);
        iVertexToParticle3D[particle->getVertexId()] = particle3D; 
    }

    void setMesh(TriangularSurfaceMesh<T>& mesh_) { mesh = mesh_; }
    plint set_cellId(plint cellId_) { cellId = cellId_; }

    TriangularSurfaceMesh<T>& getMesh() { return mesh; }
    plint get_cellId() { return cellId;  }

    plint getMpiProcessor() {
 #ifdef PLB_MPI_PARALLEL
        //  Get the individual process ID.
        int rank = MPI::COMM_WORLD.Get_rank();
#else
        int rank = 0;
#endif
        return rank;
    }

    plint getNumVertices() const { return vertices.size(); }
    plint getNumTriangles() const { return triangles.size(); }

    plint getNumVertices_Global() const { return mesh.getNumVertices(); }
    plint getNumTriangles_Global() const { return mesh.getNumTriangles(); }


    std::vector<plint> const& getTriangles() const { return triangles; }
    std::vector<plint> & getTriangles() { return triangles; }
    std::vector<plint> const& getVertices() const { return vertices; }
    std::vector<plint> & getVertices() { return vertices; }

   
    plint getVertexId(plint iTriangle, plint id) {   return mesh.getVertexId(iTriangle, id); }
    plint getNeighborVertexIds(plint iVertex) {   return mesh.getNeighborVertexIds(iVertex); }
    std::vector<plint> getNeighborTriangleIds(plint iVertex) { return mesh.getNeighborTriangleIds(iVertex); }

 
    Array<T,3> getVertex(plint iVertex) { 
        ImmersedCellParticle3D<T,Descriptor>* particle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D[iVertex]);
        return particle->get_pbcPosition(); 
    }
    Array<T,3> getVertex(plint iTriangle, plint id) { 
        return getVertex( getVertexId(iTriangle, id) ); 
    }

    Array<T,3> get_v(plint iVertex) {
        ImmersedCellParticle3D<T,Descriptor>* particle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D[iVertex]);
        return particle->get_v(); 
    }

    T get_Energy(plint iVertex) {
        ImmersedCellParticle3D<T,Descriptor>* particle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D[iVertex]);
        return particle->get_Energy(); 
    }
    Array<T,3>  get_force(plint iVertex) {
        ImmersedCellParticle3D<T,Descriptor>* particle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D[iVertex]);
        return particle->get_force(); 
    }
    Array<T,3>  getPosition(plint iVertex) {
        return iVertexToParticle3D[iVertex]->getPosition();
    }

    T computeEdgeLength(plint iVertex, plint jVertex);
    T computeEdgeLengthVector(plint iVertex, plint jVertex);
    T computeSignedAngle(plint iVertex, plint jVertex);
    T computeTriangleArea(plint iTriangle);
    Array<T,3> computeTriangleNormal(plint iTriangle);
    T computeVertexArea(plint iVertex);
    Array<T,3> computeVertexNormal(plint iVertex);
    T computeEdgeTileSpan(plint iVertex, plint jVertex);

    /* data */
    void close() {
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
    }
private:
    plint cellId;
    TriangularSurfaceMesh<T>& mesh;
    std::map<plint, Particle3D<T,Descriptor>*> iVertexToParticle3D;

    /* Quantity containers for less computations */
    std::vector<plint> triangles;
    std::vector<plint> vertices;
    /* Quantity containers for less computations */
    std::map<plint, T> edgeLengths;
    std::map<plint, Array<T,3> > edgeLengthVectors;
    std::map<plint, T> signedAngles;

    std::map<plint, T> triangleAreas;
    std::map<plint, Array<T,3> > triangleNormals;
    std::map<plint, T> vertexAreas;
    std::map<plint, Array<T,3> > vertexNormals;

    std::map<plint, T> edgeTileSpans;
};




#include "cell3D.hh"
#endif  // CELL_3D_HH

