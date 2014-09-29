#ifndef CELLFIELD_3D_H
#define CELLFIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cell3D.h"
#include "cellReductionTypes.h"

#include <set>

using namespace std;
using namespace plb;

particleEnvelopeWidth = 2;
ibmKernel = 2;


class SyncRequirements
{
public:
	SyncRequirements();
	~SyncRequirements();

	virtual std::set<plint> getSyncRequirements() { return ccrRequirements; }
	virtual std::vector<plint> getSyncRequirements() { 
		std::vector<plint> ccrRequirementsVector(ccrRequirements.begin(), ccrRequirements.end()); 
		return ccrRequirementsVector; 
	}


	virtual void insert(std::set<plint> const& ccrReq) { 
		for (std::set<plint>::const_iterator it=ccrReq.begin(); it!=ccrReq.end(); ++it) 
	    	{	ccrRequirements.insert(*it); 	}
	}

	virtual void insert(std::set<plint> const& ccrReq) { 
		for (std::set<plint>::const_iterator it=ccrReq.begin(); it!=ccrReq.end(); ++it) 
			{ ccrRequirements.insert(*it); }
	}

	virtual void insert(std::vector<plint> const& ccrReq) { 
		for (int iV = 0; iV < ccrReq.size(); ++iV)
		{ ccrRequirements.insert( ccrReq[iV] ); }
	}

private:
	std::set<plint> ccrRequirements;
};



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
	CellField3D(MultiBlockLattice3D<T, Descriptor> & lattice_, TriangularSurfaceMesh<T> const& elementaryMesh_, 
			pluint numberOfCells_, ConstitutiveModel<T, Descriptor> * cellModel_);
	~CellField3D();
public:
	/* Set or change parameters */
	void setIBMKernel(plint ibmKernel_) { ibmKernel = ibmKernel_; }
	void setIBMCoupling(bool coupleWithIBM_) { coupleWithIBM = coupleWithIBM_; }

    Cell3D<T,Descriptor> & operator[](plint cellId) { return cellIdToCell3D[cellId]; }
public:
	std::map<plint, Cell3D<T,Descriptor> > & getCellIdToCell3D() { return cellIdToCell3D; };
	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * getParticleField3D() { return immersedParticles; };
public:
	virtual void setFluidExternalForce(Array<T,3> force);
	virtual void setFluidExternalForce(T forceScalar);
	virtual void setFluidExternalForce() { return setFluidExternalForce(0.0); }

	virtual void advanceParticles();
	virtual void spreadForceIBM();
	virtual void interpolateVelocityIBM();
	virtual void applyConstitutiveModel();
	/* Need implementation */
	virtual void synchronizeCellQuantities_Local();
	virtual void synchronizeCellQuantities_Global() =0;
	virtual void synchronizeCellQuantities() { return synchronizeCellQuantities_Local(); } ;
	pluint getNumberOfCells_Global();
	pluint getNumberOfCells_Local() { return cellIdToCell3D.size(); } ;
	pluint getNumberOfCells() { return getNumberOfCells_Local(); } ;
	bool has_cellId(plint cellId) { return cellIdToCell3D.count(cellId) > 0; }
private:
	MultiBlockLattice3D<T, Descriptor> & lattice
	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * immersedParticles;
	TriangularSurfaceMesh<T> & elementaryMesh;
	ConstitutiveModel<T, Descriptor> * cellModel;
	plint ibmKernel;
	bool coupleWithIBM;

    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;

    std::map<plint, Cell3D<T,Descriptor> > cellIdToCell3D;
    std::vector<plint> cellIds;
	/* data */
	std::set<plint> ccrRequirements;
};


#include "CellField3D.hh"
#endif