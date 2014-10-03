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
    ~CellQuantityHolder() { } ;
    CellQuantityHolder(CellQuantityHolder<T> const& rhs);
    CellQuantityHolder<T>& operator=(CellQuantityHolder<T> const& rhs);

    plint const& getParticlesPerCellId() const { return particlesPerCellId; };
    plint & getParticlesPerCellId() { return particlesPerCellId; };
    void clearQuantities() ;
    void updateCQH(CellQuantityHolder<T> const& cqh);

    void insert(plint ccrId, T value)              { quantities1D[ccrId] = value; }
    void insert(plint ccrId, Array<T,3> value)     { quantities3D[ccrId] = value; }
    void insert(plint ccrId, std::vector<T> value) { quantitiesND[ccrId] = value; }
    virtual void copyFromBlockStatisticsCCR(BlockStatisticsCCR<T> & reducer);
    bool count(plint ccrId);

    std::map<plint, T >&              getQuantities1D() { return quantities1D; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> >&     getQuantities3D() { return quantities3D; }; 
    std::map<plint, std::vector<T> >& getQuantitiesND() { return quantitiesND; }; 
    std::map<plint, T > const& getQuantities1D() const { return quantities1D; };
    std::map<plint, Array<T,3> > const& getQuantities3D() const { return quantities3D; };
    std::map<plint, std::vector<T> > const& getQuantitiesND() const { return quantitiesND; };

    T const&                get1D(plint ccrId) const  { return quantities1D[ccrId]; }
    Array<T,3> const&       get3D(plint ccrId) const  { return quantities3D[ccrId]; }
    std::vector<T> const&   getND(plint ccrId) const  { return quantitiesND[ccrId]; }
    T &                     get1D(plint ccrId)   { return quantities1D[ccrId]; }
    Array<T,3> &            get3D(plint ccrId)   { return quantities3D[ccrId]; }
    std::vector<T> &        getND(plint ccrId)   { return quantitiesND[ccrId]; }
    void get(plint ccrId, T& value)              { value = get1D(ccrId); return value; }
    void get(plint ccrId, Array<T,3>& value)     { value = get3D(ccrId); return value; }
    void get(plint ccrId, std::vector<T>& value) { value = getND(ccrId); return value; }

    // numParts refer to the number of particles contained in the subdomain
    void reduceQuantity(plint ccrId, T value, plint numParts=0) ;
    void reduceQuantity(plint ccrId, Array<T,3> const& value, plint numParts=0);
    void reduceQuantity(plint ccrId, std::vector<T> const& value, plint numParts=0);

    void make_ccrId_List();

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
    virtual T & getVolume()             { return quantities1D[CCR_VOLUME];} ;
    virtual T & getSurface()            { return quantities1D[CCR_SURFACE];} ;
    virtual T & getEnergy()             { return quantities1D[CCR_ENERGY];} ;
    virtual T & getMeanAngle()          { return quantities1D[CCR_ANGLE_MEAN];} ;
    virtual T & getMeanEdgeLength()     { return quantities1D[CCR_EDGE_DISTANCE_MEAN];} ;
    virtual T & getMeanEdgeDistance()   { return quantities1D[CCR_EDGE_DISTANCE_MEAN];} ;
    virtual T & getMaxEdgeLength()      { return quantities1D[CCR_EDGE_DISTANCE_MAX];} ;
    virtual T & getMaxEdgeDistance()    { return quantities1D[CCR_EDGE_DISTANCE_MAX];} ;
    virtual T getMeanTriangleArea(plint numTriangles) { return quantities1D[CCR_SURFACE]*1.0/numTriangles;} ;
    virtual T & getMeanTileSpan()       { return quantities1D[CCR_TILE_SPAN_MEAN];} ;

    virtual Array<T,3> & getPosition() { return quantities3D[CCR_POSITION_MEAN];  } ;

    virtual Array<T,3> const& getVelocity() { return quantities3D[CCR_VELOCITY_MEAN];} ;
    virtual Array<T,3> getForce(plint numVertices=1)    { return quantities3D[CCR_FORCE] * (1.0*numVertices) / quantities1D[CCR_SURFACE];} ;
    virtual std::vector<T> & getInertia()   { return quantitiesND[CCR_INERTIA];} ;

    virtual Array<T,3> & getTumblingAngles()     { return quantities3D[CCR_TUMBLING_ANGLES];  } ;
    virtual Array<T,3> & getTankTreadingAngles() { return quantities3D[CCR_TANK_TREADING_ANGLES];  } ;
    virtual Array<T,3> & getDiameters()          { return quantities3D[CCR_DIAMETERS];  } ;
    virtual T & getSymmetryDeviation()           { return quantities1D[CCR_SYMMETRY_DEVIATION];  } ;
    virtual T & getDeformationIndex()            { return quantities1D[CCR_DEFORMATION_INDEX];  } ;
    virtual T & getTaylorDeformationIndex()      { return quantities1D[CCR_TAYLOR_DEFORMATION_INDEX];  } ;

};



template<typename T, template<typename U> class Descriptor>
class Cell3D ;


template<typename T, template<typename U> class Descriptor>
void computeCCRQuantities(plint ccrId, BlockStatisticsCCR<T> & reducer, Cell3D<T, Descriptor> * cell, plint iVertex);


template<typename T, template<typename U> class Descriptor>
class Cell3D : public CellQuantityHolder<T>
{
public:
    Cell3D(TriangularSurfaceMesh<T> & mesh_, plint cellId_=-1);
    Cell3D(Cell3D<T,Descriptor> const& rhs);
    ~Cell3D() {};
    Cell3D<T,Descriptor>& operator=(Cell3D<T,Descriptor> const& rhs) {
        CellQuantityHolder<T>::operator=(rhs),
        mesh = rhs.mesh;
        cellId=rhs.cellId;
        setMesh();
        return *this;
     }

    void push_back(Particle3D<T,Descriptor>* particle3D);
    void close();

    void setMesh() ;
    void set_cellId(plint cellId_) { cellId = cellId_; }

    TriangularSurfaceMesh<T> const& getMesh() { return mesh; }
    plint & get_cellId() { return cellId;  }

    plint getMpiProcessor();

    plint getCellNumVertices() const { return cellNumVertices; }
    plint getCellNumTriangles() const { return cellNumTriangles; }
    plint getNumVertices_Local() const { return vertices.size(); }
    plint getNumTriangles_Local() const { return triangles.size(); }
    plint getNumVertices_Global() const { return cellNumVertices; }
    plint getNumTriangles_Global() const { return cellNumTriangles; }


    std::vector<plint> const& getTriangles() const { return triangles; }
    std::vector<plint> & getTriangles() { return triangles; }
    std::vector<plint> const& getVertices() const { return vertices; }
    std::vector<plint> & getVertices() { return vertices; }
    plint getEdgeId(plint iVertex, plint jVertex);

    std::vector<Array<plint,2> > const& getEdges() const { return edges; }
    std::vector<Array<plint,2> > & getEdges() { return edges; }

    std::map<plint, plint> getInvertVertices() {
        std::map<plint, plint> iv;
        for (int var = 0; var < vertices.size(); ++var) {
            iv[ vertices[var] ] = var;
        }
        return iv;
    }
 
    plint getVertexId(plint iTriangle, plint id) {   return mesh.getVertexId(iTriangle, id); }
    std::vector<plint> getNeighborVertexIds(plint iVertex) {   return mesh.getNeighborVertexIds(iVertex); }
    std::vector<plint> getNeighborVertexIds(plint iVertex, plint jVertex) {   return mesh.getNeighborVertexIds(iVertex, jVertex); }
    std::vector<plint> getNeighborTriangleIds(plint iVertex) { return mesh.getNeighborTriangleIds(iVertex); }
    std::vector<plint> getAdjacentTriangleIds(plint iVertex, plint jVertex) { return mesh.getAdjacentTriangleIds(iVertex, jVertex); }
 
    Array<T,3> getVertex(plint iVertex) { return castParticleToICP3D(iVertexToParticle3D[iVertex])->get_pbcPosition(); }
    Array<T,3> getVertex(plint iTriangle, plint id) { return getVertex( getVertexId(iTriangle, id) ); }
    Particle3D<T,Descriptor>* getParticle3D(plint iVertex) { return iVertexToParticle3D[iVertex]; }

    Array<T,3> get_v(plint iVertex) { return castParticleToICP3D(iVertexToParticle3D[iVertex])->get_v(); }

    T get_Energy(plint iVertex) { return castParticleToICP3D(iVertexToParticle3D[iVertex])->get_Energy(); }
    Array<T,3>  get_force(plint iVertex) { return castParticleToICP3D(iVertexToParticle3D[iVertex])->get_force(); }
    Array<T,3>  getPosition(plint iVertex) { return iVertexToParticle3D[iVertex]->getPosition();}
    Array<T,3> get_pbcPosition(plint iVertex) { return castParticleToICP3D(iVertexToParticle3D[iVertex])->get_pbcPosition(); }

    T computeEdgeLength(plint iVertex, plint jVertex);
    Array<T,3> computeEdgeLengthVector(plint iVertex, plint jVertex);
    T computeSignedAngle(plint iVertex, plint jVertex, bool& found);
    T computeSignedAngle(plint iVertex, plint jVertex, plint & kVertex, plint & lVertex, bool& found);
    T computeTriangleArea(plint iTriangle);
    Array<T,3> computeTriangleNormal(plint iTriangle);
    plint findTriangleId(plint iVertex, plint jVertex, plint kVertex) ;
    T computeTriangleArea(plint iVertex, plint jVertex, plint kVertex) { return computeTriangleArea(findTriangleId(iVertex,jVertex,kVertex)); };
    Array<T,3> computeTriangleNormal(plint iVertex, plint jVertex, plint kVertex) { return computeTriangleNormal(findTriangleId(iVertex,jVertex,kVertex)); };
    T computeVertexArea(plint iVertex);
    Array<T,3> computeVertexNormal(plint iVertex);
    T computeEdgeTileSpan(plint iVertex, plint jVertex);
    virtual void clearReducer() { reducer.clear();  }
    virtual void computeCCRQuantities(plint ccrId, plint iVertex) { calculateCCRQuantities(ccrId, reducer, this, iVertex); }
    virtual void computeCCRQuantities(plint ccrId, Particle3D<T,Descriptor> * particle) { calculateCCRQuantities(ccrId, reducer, this, castParticleToICP3D(particle)->getVertexId()); }
    void closeCCRQuantities() { this->copyFromBlockStatisticsCCR(reducer); }
public:
    Array<T,3> getForce()    { return this->getForce(cellNumVertices); } ;
    Array<T,3> & getPosition() { return this->get3D(CCR_POSITION_MEAN);  } ;

private:
    /* data */
    TriangularSurfaceMesh<T> & mesh;
    plint cellId;
    std::map<plint, Particle3D<T,Descriptor>*> iVertexToParticle3D;

    BlockStatisticsCCR<T> reducer;

    pluint cellNumVertices, cellNumTriangles, cellNumEdges;
    /* Quantity containers for less computations */
    std::vector<plint> triangles;
    std::vector<plint> vertices;
    std::vector<Array<plint,2> > edges;
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

