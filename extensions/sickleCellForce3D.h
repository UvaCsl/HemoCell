#ifndef SICKLE_CELL_3D_H
#define SICKLE_CELL_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "surfaceParticle3D.h"
using namespace std;
using namespace plb;

template<typename T, template<typename U> class Descriptor>
class SickleCellForce3D
{
public:
    SickleCellForce3D(HemoCellField & cellField_, Array<T,3> forceLeft, Array<T,3> forceRight, T percentage) :
        cellField(cellField_)
    {
        TriangularSurfaceMesh<T> & mesh = cellField.getMesh();
        plint numVertices = mesh.getNumVertices();
        nparPerSide = ceil( percentage * numVertices);
        Array<T,3> forceForRest = (forceLeft + forceRight) * (-1.0/(numVertices - nparPerSide));
        // The forceForRest is added to the others in order to avoid making a complex list without those
        Array<T,3> forceLeftPerVertex  = forceLeft * (1.0 / nparPerSide)  - forceForRest;
        Array<T,3> forceRightPerVertex = forceRight * (1.0 / nparPerSide) - forceForRest;

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

        cellIds.push_back(0); cellIds.push_back(0); cellIds.push_back(0); // Three cellIds. One for left, one for right, one for the rest (basically, all)
        forces.push_back(forceLeftPerVertex);
        forces.push_back(forceRightPerVertex);
        forces.push_back(forceForRest);  // Three forces. One for left, one for right, one for the rest (basically, all)
        vertices.clear();
        vertices.resize(3);
        for (plint i = 0; i < nparPerSide; ++i) {
        	vertices[0].push_back(iv2X[i].first);
        	vertices[1].push_back(iv2X[numVertices - i - 1].first);
        }
        for (plint i = 0; i < numVertices; ++i) {
        	vertices[2].push_back(i);
        }
    }

    ~SickleCellForce3D() {
    } ;

    void stretch() {
        applyProcessingFunctional (
            new ApplyForce3D<T,Descriptor>(cellField, cellIds, vertices, forces),
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
private:
    HemoCellField & cellField;
    plint nparPerSide;
    std::vector<plint> cellIds;
    std::vector<std::vector<plint> > vertices;
    std::vector<Array<T,3> > forces;
};



#include "sickleCellForce3D.hh"
#endif  // SICKLE_CELL_3D_H

