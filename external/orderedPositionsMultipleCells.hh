#ifndef ORDERED_POSISIONS_OF_MULTIPLE_CELLS_HH
#define ORDERED_POSISIONS_OF_MULTIPLE_CELLS_HH
#include "orderedPositionsMultipleCells.h"
#include <algorithm>    // std::random_shuffle


template<typename T>
Array<T,3> randomMeshRotationAngleAndBoundingBox(TriangularSurfaceMesh<T>* mesh, T & dx, T & dy, T & dz) {
    Array<T,3> randomAngles;
    randomAngles[0] = guessRandomNumber()* 2 * pi;
    randomAngles[1] = guessRandomNumber()* pi;
    randomAngles[2] = guessRandomNumber()* 2 * pi;

    ElementsOfTriangularSurfaceMesh<T> emptyEoTSM;
    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*mesh, emptyEoTSM);

    meshRandomRotation(meshCopy, randomAngles);
    Array<T,2> xRange, yRange, zRange;
    meshCopy->computeBoundingBox (xRange, yRange, zRange);
    dx = (xRange[1] - xRange[0]);
    dy = (yRange[1] - yRange[0]);
    dz = (zRange[1] - zRange[0]);

    delete meshCopy;
    return randomAngles;
}

template<typename T>
void getOrderedPositionsMultipleCellsVector(Box3D realDomain, 
    std::vector<TriangularSurfaceMesh<T>* > & meshes,
    std::vector<plint> & Np,
    std::vector<std::vector<Array<T,3> > > & positions, 
    std::vector<std::vector<plint> > & cellIds,
    std::vector<std::vector<Array<T,3> > > & randomAngles) {
    PLB_PRECONDITION( meshes.size()==Np.size() );


	int fac=10; 	// We scale everything by a factor of 10,
					// because the solver only works with integers.
	float dl=0.1;
	int W=realDomain.getNx() * fac,
		H=realDomain.getNy() * fac,
		D=realDomain.getNz() * fac;
	int nodelimit=0, iterlimit=0, timelimit=10;
	int packingtype=0;
	int ub, lb;
	int nodeused, iterused, timeused;

	int n=0;
	for (pluint i = 0; i < Np.size(); ++i)	{ n+=Np[i];	}
	std::vector<plint> indices(n);
	for (pluint i = 0; i < n; ++i)	{ indices[i]=i;	}
	std::random_shuffle ( indices.begin(), indices.end() );

	Array<T,3> randomAngle(0.,0.,0.);
	int w[n], h[n], d[n];
	int x[n], y[n], z[n], bno[n];
    int xt, yt,zt;
    Array<T,2> xRange, yRange, zRange;
    T dim[3];
	plint indx, ni=0;
	T dx, dy, dz;

	randomAngles.clear();  randomAngles.resize(Np.size());
	for (pluint i = 0; i < Np.size(); ++i)
	{
        meshes[i]->computeBoundingBox (xRange, yRange, zRange);
        dx = (xRange[1] - xRange[0]);
        dy = (yRange[1] - yRange[0]);
        dz = (zRange[1] - zRange[0]);
        randomAngles[i].resize(Np[i]);
		for (pluint j = 0; j < Np[i]; ++j)
		{
//	        randomAngle = randomMeshRotationAngleAndBoundingBox(meshes[i], dx, dy, dz);
	        randomAngles[i][j] = randomAngle;
	        dim[0] = dx + 2*dl; // Give some air to breath with dl LU
	        dim[1] = dy + 2*dl;
	        dim[2] = dz + 2*dl;
	        xt=int(ceil(dim[0] * fac));
	        yt=int(ceil(dim[1] * fac));
	        zt=int(ceil(dim[2] * fac));
			indx=indices[ni];
			w[indx]=xt;
			h[indx]=yt;
			d[indx]=zt;
			ni++;
		}
		printf("%d) %d with size (%d,%d,%d) \n", i, Np[i], xt, yt,zt);

	}

	printf("Total N=%d WHD (%d,%d,%d) \n", n, W, H, D);

	binpack3d(n, W, H, D,
               w, h, d, 
               x, y, z, bno,
               &lb, &ub, 
              nodelimit, iterlimit, timelimit, 
              &nodeused, &iterused, &timeused, 
              packingtype);
	pcout << "Packing Done." << std::endl;

	positions.clear();	positions.resize(Np.size());
	cellIds.clear();	cellIds.resize(Np.size());
	ni=0;
	int mpiRank = global::mpi().getRank();
	Array<T,3> rdomain(realDomain.x0, realDomain.y0, realDomain.z0);
	for (pluint i = 0; i < Np.size(); ++i)
	{
		positions[i].resize(Np[i]);
		cellIds[i].resize(Np[i]);
		for (pluint j = 0; j < Np[i]; ++j)
		{
			plint indx=indices[ni];
			positions[i][j] = rdomain + Array<T,3>(x[indx]*1.0/fac + dl, y[indx]*1.0/fac + dl, z[indx]*1.0/fac + dl);
			cellIds[i][j] = ni + 1000*mpiRank;
			ni++;
//			printf("pack (%f,%f,%f) \n", positions[i][j][0], positions[i][j][1], positions[i][j][2]);
		}
	}

}







/* ******** OrderedPositionMultipleCellField3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    pluint numberOfCellFields=cellFields.size();
    PLB_PRECONDITION( blocks.size()==(numberOfCellFields+1) );
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();
    srand (mpiSize * mpiRank);
    T ratio;
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);

    std::vector<ParticleField3D<T,Descriptor>* > particleFields(numberOfCellFields);
    std::vector<T> volumes(numberOfCellFields);
    std::vector<T> volumeFractions(numberOfCellFields);
    std::vector<TriangularSurfaceMesh<T>* > meshes(numberOfCellFields);
    std::vector<ElementsOfTriangularSurfaceMesh<T> > emptyEoTSM(numberOfCellFields);
    std::vector<Box3D> relativeDomains(numberOfCellFields);
	std::vector<plint> Np(numberOfCellFields);

    T totalVolumeFraction=0;

    T Vdomain = domain.getNx() * domain.getNy() * domain.getNz();

    Dot3D fLocation(fluid.getLocation());
    Box3D realDomain(
                        domain.x0 + fLocation.x, domain.x1 + fLocation.x,
                        domain.y0 + fLocation.y, domain.y1 + fLocation.y,
                        domain.z0 + fLocation.z, domain.z1 + fLocation.z );
    Array<T,2> xRange, yRange, zRange;

    for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
        T vf = cellFields[iCF]->getVolumeFraction();
        if (vf > 1.0) { vf /= 100; }
        volumeFractions[iCF] = vf;

        TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(cellFields[iCF]->getMesh(), emptyEoTSM[iCF]);
        mesh->computeBoundingBox (xRange, yRange, zRange);
        mesh->translate(Array<T,3>(-xRange[0], -yRange[0], -zRange[0]));
        meshes[iCF] = mesh;
        volumes[iCF] = MeshMetrics<T>(*mesh).getVolume();

        Np[iCF] = (1 + 0.05*(guessRandomNumber() -0.5)) * (Vdomain/volumes[iCF] * volumeFractions[iCF]) ;
        totalVolumeFraction += volumes[iCF];
        particleFields[iCF] = ( dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[iCF+1]) );
        particleFields[iCF]->removeParticles(particleFields[iCF]->getBoundingBox());
        Dot3D relativeDisplacement = computeRelativeDisplacement(fluid, *(particleFields[iCF]));
        relativeDomains[iCF] = ( Box3D(
                            domain.x0 + relativeDisplacement.x, domain.x1 + relativeDisplacement.x,
                            domain.y0 + relativeDisplacement.y, domain.y1 + relativeDisplacement.y,
                            domain.z0 + relativeDisplacement.z, domain.z1 + relativeDisplacement.z ) );
    }

    std::vector<std::vector<Array<T,3> > > positions;
    std::vector<std::vector<plint> > cellIds;
    std::vector<std::vector<Array<T,3> > > randomAngles;
    getOrderedPositionsMultipleCellsVector(realDomain, meshes, Np, positions, cellIds, randomAngles);

    for (pluint iCF = 0; iCF < positions.size(); ++iCF)
    {
        cout << "(OrderedPositionMultipleCellField3D) ";
    	for (pluint c = 0; c < positions[iCF].size(); ++c)
    	{
//    	    ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
//    	    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);
//    	    meshRandomRotation(meshCopy, randomAngles[iCF][c]);
			positionCellInParticleField(*(particleFields[iCF]), fluid, 
			        meshes[iCF], positions[iCF][c]-0.5, cellIds[iCF][c]);
//			delete meshCopy;
    	}

	    // DELETE CELLS THAT ARE NOT WHOLE
    	plint nVertices=meshes[iCF]->getNumVertices();
        cout << mpiRank;
	    plint cellsDeleted = deleteIncompleteCells(*(particleFields[iCF]), fluid, relativeDomains[iCF], nVertices);
	    std::vector<Particle3D<T,Descriptor>*> particles;
	    particleFields[iCF]->findParticles(particleFields[iCF]->getBoundingBox(),   particles);
		cout    << ", number of particles/nVertices " << particles.size()*1.0/nVertices
		        << " (deleted:" << cellsDeleted << ") for particleId:" << iCF << std::endl;
		delete meshes[iCF];
	}
}


template<typename T, template<typename U> class Descriptor>
OrderedPositionMultipleCellField3D<T,Descriptor>* OrderedPositionMultipleCellField3D<T,Descriptor>::clone() const {
    return new OrderedPositionMultipleCellField3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT OrderedPositionMultipleCellField3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}





template<typename T, template<typename U> class Descriptor>
void orderedPositionMultipleCellField3D(std::vector<CellField3D<T, Descriptor>* > & cellFields) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D*> fluidAndParticleFieldsArg;
    fluidAndParticleFieldsArg.push_back( &(cellFields[0]->getFluidField3D()) );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back( &(cellFields[icf]->getParticleField3D()) );
    }

    std::vector<MultiBlock3D*> particleFieldsArg;
    particleFieldsArg.push_back( &(cellFields[0]->getParticleField3D()) );

    applyProcessingFunctional (
        new OrderedPositionMultipleCellField3D<T,Descriptor>(cellFields),
        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );
    pcout << "Ready to start.." << std::endl;

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        cellFields[icf]->advanceParticles();
        cellFields[icf]->synchronizeCellQuantities();
    }

    global::timer("CellInit").stop();
}

#endif // ORDERED_POSISIONS_OF_MULTIPLE_CELLS_HH


