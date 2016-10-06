#ifndef CELL_STRETCHING_3D_H
#define CELL_STRETCHING_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "surfaceParticle3D.h"
#include <math.h>
using namespace std;
using namespace plb;


template<typename T>
bool comparePair (std::pair <plint,T> p1, std::pair <plint,T> p2) {
    T i = p1.second;
    T j = p2.second;
    return (i<j);
}

template<typename T>
std::vector<plint> meshVerticesFromDirection(TriangularSurfaceMesh<T> & mesh, plint direction, plint signDirection, T percentage) {
    /*
     * Returns the vertices of a mesh that are from a specific direction --> (0,1,2) -- (x,y,z).
     * The sign on direction denotes if it's the "left" (lower values of direction) or the "right" (higher values of direction).
     */
    std::vector<plint> verticesToReturn;
    plint numVertices = mesh.getNumVertices() ;
    plint pcNumVertices = ceil(numVertices * percentage);
    std::vector<std::pair <plint,T> > iv2;

    for (plint iV = 0; iV < numVertices; ++iV) {
        Array<T,3> vertex = mesh.getVertex(iV);
        iv2.push_back(  std::make_pair(iV, vertex[direction]) );
    }
    std::sort(iv2.begin(), iv2.end(), comparePair<T>);

    for (plint i = 0; i < pcNumVertices; ++i) {
        if (signDirection >0)  { verticesToReturn.push_back(iv2[i].first); }
        else { verticesToReturn.push_back(iv2[numVertices - i - 1].first); }
    }
    return verticesToReturn;
}


template<typename T, template<typename U> class Descriptor>
void applyForceToCells(CellField3D<T, Descriptor> & cellField_,
              std::vector<plint> const& cellIds_, std::vector<std::vector<plint> > iVertices_,
              std::vector<Array<T,3> > const& forces_);


template<typename T, template<typename U> class Descriptor>
void applyForceToCells(CellField3D<T, Descriptor> & cellField_,
              std::vector<plint> const& cellIds_,
              std::vector<Array<T,3> > const& forces_);

template<typename T, template<typename U> class Descriptor>
class ApplyForce3D ;

template<typename T, template<typename U> class Descriptor>
class CellStretch
{
public:
    CellStretch(CellField3D<T, Descriptor> & cellField_, T force_, T percentage) :
        cellField(cellField_)
    {
        TriangularSurfaceMesh<T> & mesh = cellField.getMesh();
        plint numVertices = mesh.getNumVertices();
        nparPerSide = ceil( percentage * numVertices );
        forcePerSide = force_ / 2.0;   // Totle force / 2 for each side
        std::vector<std::pair <plint,T> > iv2X, iv2Y, iv2Z;
        for (plint iV = 0; iV < numVertices; ++iV) {
            Array<T,3> vertex = mesh.getVertex(iV);
            iv2X.push_back(  std::make_pair(iV, vertex[0]) );
            iv2Y.push_back(  std::make_pair(iV, vertex[1]) );
            iv2Z.push_back(  std::make_pair(iV, vertex[2]) );
        }
        std::sort(iv2X.begin(), iv2X.end(), comparePair<T>);
        std::sort(iv2Y.begin(), iv2Y.end(), comparePair<T>);
        std::sort(iv2Z.begin(), iv2Z.end(), comparePair<T>);

        cellIds.push_back(0); cellIds.push_back(0); // Two cellIds. One for left and one for right
        forces.push_back(Array<T,3>(-forcePerSide, 0, 0) ); forces.push_back(Array<T,3>(forcePerSide, 0, 0) ); // Two forces. One for left and one for right
        xVertices.clear();         yVertices.clear();         zVertices.clear();
        xVertices.resize(2);         yVertices.resize(2);         zVertices.resize(2);
        for (plint i = 0; i < nparPerSide; ++i) {
            xVertices[0].push_back(iv2X[i].first);
            xVertices[1].push_back(iv2X[numVertices - i - 1].first);
            yVertices[0].push_back(iv2Y[i].first);
            yVertices[1].push_back(iv2Y[numVertices - i - 1].first);
            zVertices[0].push_back(iv2Z[i].first);
            zVertices[1].push_back(iv2Z[numVertices - i - 1].first);
        }
    }

    ~CellStretch() {
    } ;

    void stretch() {
        applyProcessingFunctional (
            new ApplyForce3D<T,Descriptor>(cellField, cellIds, xVertices, forces),
            cellField.getBoundingBox(), cellField.getParticleArg() );
    }

    // Measures stretch as the Sx = max{x_i} - min{x_i}
    Array<T,3> measureStretch() {
        SyncRequirements sr(CCR_POSITION_MIN);
        sr.insert(CCR_POSITION_MAX);
        cellField.synchronizeSyncRequirements(sr);
        Array<T,3> stretch = cellField[0]->get3D(CCR_POSITION_MAX) - cellField[0]->get3D(CCR_POSITION_MIN);
        return stretch;
    }

    void setCellIds(std::vector<plint> cellIds_) { cellIds = cellIds_; } ;
    void setVertices(std::vector<std::vector<plint> > iVertices_)  { xVertices = iVertices_;} ;
    void setForces(std::vector<Array<T,3> > forces_)  { forces = forces_; } ;
    void setForce(T force_)  { forcePerSide = force_; } ;
private:
    CellField3D<T, Descriptor> & cellField;
    T forcePerSide;
    plint nparPerSide;
    std::vector<plint> cellIds;
    std::vector<std::vector<plint> > xVertices, yVertices, zVertices;
    std::vector<Array<T,3> > forces;
};




template<typename T, template<typename U> class Descriptor>
class ApplyForce3D : public BoxProcessingFunctional3D
{
public:
    ApplyForce3D (CellField3D<T, Descriptor> & cellField_,
                  std::vector<plint> const& cellIds_, std::vector<std::vector<plint> > iVertices_,
                  std::vector<Array<T,3> > const& forces_);
    ApplyForce3D (CellField3D<T, Descriptor> & cellField_,
                  std::vector<plint> const& cellIds_,
                  std::vector<Array<T,3> > const& forces_);
    ~ApplyForce3D() {
//        std::cout <<" ~ApplyForce3D() " << global::mpi().getRank() << std::endl;
    } ;
    ApplyForce3D(ApplyForce3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ApplyForce3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    CellField3D<T, Descriptor> & cellField;
    std::vector<plint> const& cellIds;
    std::vector<std::vector<plint> > iVertices;
    std::vector<Array<T,3> > const& forces;
    std::vector<bool > allVerticesOfTheCell;
};


#include "cellStretching3D.hh"
#endif  // CELL_STRETCHING_3D_H

