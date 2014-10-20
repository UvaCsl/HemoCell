#ifndef FICSION_HELPER_FUNCTIONALS_HH
#define FICSION_HELPER_FUNCTIONALS_HH

#include "helperFunctionals.h"



/* ******** PositionBoundaryParticles *********************************** */
template<typename T, template<typename U> class Descriptor>
PositionBoundaryParticles<T,Descriptor>::PositionBoundaryParticles () { };

template<typename T, template<typename U> class Descriptor>
void PositionBoundaryParticles<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
#ifdef PLB_MPI_PARALLEL
    //  Get the number of processes.
    int  p = MPI::COMM_WORLD.Get_size();
    //  Get the individual process ID.
     int id = MPI::COMM_WORLD.Get_rank();
#else
     int p=1;
     int id = 0;
#endif
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& boundaryParticleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);

    Dot3D latticeLocation= fluid.getLocation();
    Array<T,3> relativePosition(latticeLocation.x, latticeLocation.y, latticeLocation.z);

    plint vertexId=0, cellId=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            	if (fluid.get(iX, iY, iZ).getDynamics().isBoundary()) {
            		Array<T,3> vertex = Array<T,3>(iX,iY,iZ) + relativePosition;
//            		std::cout << "(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ") " << id << std::endl;
            		std::vector<Array<T,3> > vertices;
					vertices.push_back(Array<T,3>(0.25, 0.25, 0.25) + vertex);
					vertices.push_back(Array<T,3>(0.25, 0.25, 0.50) + vertex);
					vertices.push_back(Array<T,3>(0.25, 0.50, 0.50) + vertex);
					vertices.push_back(Array<T,3>(0.25, 0.50, 0.25) + vertex);
					vertices.push_back(Array<T,3>(0.50, 0.50, 0.25) + vertex);
					vertices.push_back(Array<T,3>(0.50, 0.25, 0.25) + vertex);
					vertices.push_back(Array<T,3>(0.50, 0.25, 0.50) + vertex);
					vertices.push_back(Array<T,3>(0.50, 0.50, 0.50) + vertex);
            	    for (plint indx=0; indx<8; ++indx) {
            	    	boundaryParticleField.addParticle(
            	    			boundaryParticleField.getBoundingBox(),
            	    			new ImmersedCellParticle3D<T,Descriptor>(vertexId, vertices[indx], cellId));
            	    }
            	}
            }
        }
    }
    std::vector<Particle3D<T,Descriptor>*> particles;
    boundaryParticleField.findParticles(domain, particles);
	std::cout << id << " Number of Particles " << particles.size() << std::endl;

}

template<typename T, template<typename U> class Descriptor>
PositionBoundaryParticles<T,Descriptor>* PositionBoundaryParticles<T,Descriptor>::clone() const {
    return new PositionBoundaryParticles<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void PositionBoundaryParticles<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field
    modified[1] = modif::nothing; // Fluid field
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT PositionBoundaryParticles<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}




#endif  // FICSION_HELPER_FUNCTIONALS_HH
