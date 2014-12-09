#ifndef CELLFIELD_3D_H
#define CELLFIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cell3D.h"
#include "cellModel3D.h"
#include "shapeMemoryModel3D.h"
#include "cellReductionTypes.h"
#include "cellFieldFunctionals3D.h"
#include "cellCellForces3D.h"
#include <set>
#include <string>
#include <algorithm>    // std::sort


using namespace std;
using namespace plb;


template<typename T, template<typename U> class Descriptor>
class CellField3D
{
public:
	/*    
	Declaration:
	-------------
	CellField3D<T, DESCRIPTOR> RBCField(lattice, Cells.getMesh(), npar, cellModel);
	
	Arguments:
	-------------
			MultiBlockLattice3D<T, Descriptor> & lattice. Contents may be changed (like IBM Force spreading)
			ConstitutiveModel<T, Descriptor> * cellModel. Deletion is taken care by CellField3D
			pluint numberOfCells_. Global numer of Cells
	*/
	CellField3D(MultiBlockLattice3D<T, Descriptor> & lattice_, TriangularSurfaceMesh<T> & elementaryMesh_,
			T hematocrit_, ConstitutiveModel<T, Descriptor> * cellModel_, plint ibmKernel_, std::string identifier_);
	~CellField3D();
public:
	/* Set or change parameters */
	void setIBMKernel(plint ibmKernel_) { ibmKernel = ibmKernel_; }
	void setIBMCoupling(bool coupleWithIBM_) { coupleWithIBM = coupleWithIBM_; }
	TriangularSurfaceMesh<T> & getMesh() { return elementaryMesh; }

	Cell3D<T,Descriptor>* operator[](plint cellId) {
	    if(count(cellId) > 0) {
	        return cellIdToCell3D[cellId];
	    } else {
	        return NULL;
	    }
	}
public:
	std::map<plint, Cell3D<T,Descriptor>* > & getCellIdToCell3D() { return cellIdToCell3D; };
	MultiBlockLattice3D<T, Descriptor> & getFluidField3D() { return lattice; };
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & getParticleField3D() { return *immersedParticles; };
    std::string getIdentifier() { return identifier; };
    T getVolumeFraction() { return hematocrit; };
    void setParticleField3D(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * immersedParticles_) {
        delete immersedParticles;
        immersedParticles=immersedParticles_;
        particleArg.clear(); particleLatticeArg.clear(); particleReductionParticleArg.clear();
        particleArg.push_back(immersedParticles);
        particleLatticeArg.push_back(immersedParticles);
        particleLatticeArg.push_back(&lattice);
        particleReductionParticleArg.push_back(immersedParticles);
        particleReductionParticleArg.push_back(reductionParticles);
    };
public:
	virtual void setFluidExternalForce(Array<T,3> force);
	virtual void setFluidExternalForce(T forceScalar);
	virtual void setFluidExternalForce() { return setFluidExternalForce(0.0); }

    virtual void initialize(std::vector<Array<T,3> > & centers);
    virtual void grow(plint growIterations=0);
    virtual void createCellMap();
	virtual void advanceParticles();
	virtual void spreadForceIBM();
    virtual void computeVelocity(T ratio=1.0);
    virtual void interpolateVelocityIBM();
    virtual void applyConstitutiveModel();
    virtual void applyCellCellForce(CellCellForce3D<T> & calcForce_, T cutoffRadius_);
    virtual void applyDifferentCellForce(CellCellForce3D<T> & calcForce, T cutoffRadius,
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * otherCellParticles);
	/* Need implementation */
    virtual void synchronizeSyncRequirements_Local(SyncRequirements ccrRequirements_);
    virtual void synchronizeSyncRequirements_Global(SyncRequirements ccrRequirements_) {};
    virtual void synchronizeSyncRequirements(SyncRequirements ccrRequirements_) { return synchronizeSyncRequirements_Local(ccrRequirements_); } ;
    virtual void synchronizeCellQuantities_Local(SyncRequirements ccrRequirements_=SyncRequirements());
	virtual void synchronizeCellQuantities_Global(SyncRequirements ccrRequirements_=SyncRequirements()) {} ;
	virtual void synchronizeCellQuantities(SyncRequirements ccrRequirements_=SyncRequirements()) { return synchronizeCellQuantities_Local(ccrRequirements_); } ;
	pluint getNumberOfCells_Global() {
	    CountGlobalNumberOfCells<T,Descriptor> nfunctional( elementaryMesh.getNumVertices() );
	    applyProcessingFunctional(nfunctional, immersedParticles->getBoundingBox(), particleArg);
	    pluint gnoc =  nfunctional.getValue();
	    return gnoc;
	};
	pluint getNumberOfCells_Local() { return cellIdToCell3D.size(); } ;
	pluint getNumberOfCells() { return getNumberOfCells_Local(); } ;
    plint count(plint cellId) { return cellIdToCell3D.count(cellId); }
    std::vector<plint> getCellIds() {
        std::vector<plint> cellIds;
        typename std::map<plint, Cell3D<T,Descriptor>* >::iterator itrtr;
        for (itrtr  = cellIdToCell3D.begin(); itrtr != cellIdToCell3D.end(); ++itrtr) {
            cellIds.push_back(itrtr->first);
        }
        std::sort (cellIds.begin(), cellIds.end());
        return cellIds;
    }
    bool has_cellId(plint cellId) { return cellIdToCell3D.count(cellId) > 0; }
private:
	MultiBlockLattice3D<T, Descriptor> & lattice;
	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * immersedParticles;
	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * reductionParticles;
	TriangularSurfaceMesh<T> & elementaryMesh;
    T hematocrit;
	ConstitutiveModel<T, Descriptor> * cellModel;
	plint ibmKernel;
	bool coupleWithIBM;
    SyncRequirements ccrRequirements;
    std::string identifier;

public:
    std::vector<MultiBlock3D*> & getReductionParticleArg()                 { return reductionParticleArg; }
    std::vector<MultiBlock3D*> & getParticleArg()                 { return particleArg; }
    std::vector<MultiBlock3D*> & getParticleLatticeArg()          { return particleLatticeArg; }
    std::vector<MultiBlock3D*> & getparticleReductionParticleArg() { return particleReductionParticleArg; }
    Box3D getBoundingBox() { return immersedParticles->getBoundingBox(); }
private:
    std::vector<MultiBlock3D*> reductionParticleArg;
    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;
    std::vector<MultiBlock3D*> particleReductionParticleArg;
    std::map<plint, Cell3D<T,Descriptor>* > cellIdToCell3D;
    std::vector<plint> cellIds;
	/* data */
};


#include "cellField3D.hh"
#endif

