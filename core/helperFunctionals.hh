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

                Array<T,3> vertex = Array<T,3>(iX,iY,iZ) + relativePosition;

                if(fluid.get(iX, iY, iZ).getDynamics().isBoundary())
                {
                    for (int px = -1; px <= +1; ++px){ for (int py = -1; py <= +1; ++py) { for (int pz = -1; pz <= +1; ++pz) {
                        
                        if(abs(px)+abs(py)+abs(pz) != 1) 
                            continue;

                        try{
                            if(!fluid.get(px+iX, py+iY, pz+iZ).getDynamics().isBoundary()) {
                                cellId +=1;
                                boundaryParticleField.addParticle(
                                    boundaryParticleField.getBoundingBox(),
                                    new SurfaceParticle3D(Array<T,3>(0.5+(px*0.5), 0.5+(py*0.5), 0.5+(pz*0.5)) + vertex, cellId, 0));
                            }
                        }
                        catch (int e) { }

                    } } }

                }

                // -------------------
                /*
                // Put boundary particles in the first fluid layer
                //  rational: deny particles from the outer shear layer -> glycocalyx

                if(!fluid.get(iX, iY, iZ).getDynamics().isBoundary())
                {
                    bool neighboringBoundariesAnywhere=false;

                    for (int px = iX-1; px <= iX+1; ++px) {  for (int py = iY-1; py <= iY+1; ++py) { for (int pz = iZ-1; pz <= iZ+1; ++pz) {
                        try {
                            neighboringBoundariesAnywhere = neighboringBoundariesAnywhere or fluid.get(px, py, pz).getDynamics().isBoundary();
                        } catch (int e) { neighboringBoundariesAnywhere=true; }
                        
                    }  }  }


                    if(neighboringBoundariesAnywhere)
                    {
                        Array<T,3> vertex = Array<T,3>(iX,iY,iZ) + relativePosition;

                        std::vector<Array<T,3> > vertices;
                        cellId +=1;
                        boundaryParticleField.addParticle(
                                boundaryParticleField.getBoundingBox(),
                                new SurfaceParticle3D(Array<T,3>(0.5, 0.5, 0.5) + vertex, cellId, 0));
                    }
                }
                */

                // ----------------------

                /*
                if (fluid.get(iX, iY, iZ).getDynamics().isBoundary() and (not neighboringBoundariesEverywhere) ) {
            		
            		
                    // Position 8 particles inside the fluid cell
                    
                    T step= 0.25;
                    vertices.push_back(Array<T,3>(step, step, step) + vertex);
                    vertices.push_back(Array<T,3>(step, step, step*3) + vertex);
                    vertices.push_back(Array<T,3>(step, step*3, step*3) + vertex);
                    vertices.push_back(Array<T,3>(step, step*3, step) + vertex);
                    vertices.push_back(Array<T,3>(step*3, step*3, step) + vertex);
                    vertices.push_back(Array<T,3>(step*3, step, step) + vertex);
                    vertices.push_back(Array<T,3>(step*3, step, step*3) + vertex);
                    vertices.push_back(Array<T,3>(step*3, step*3, step*3) + vertex);
                    */

                    // Position 6 particles face centered, if they face a fluid cell
                    /*
                    T step = 0.5;
                    if(!fluid.get(iX, iY, iZ-1).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(step, step, 0) + vertex);

                    if(!fluid.get(iX, iY-1, iZ).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(step, 0, step) + vertex);

                    if(!fluid.get(iX-1, iY, iZ).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(0, step, step) + vertex);

                    if(!fluid.get(iX+1, iY, iZ).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(2*step, step, step) + vertex);

                    if(!fluid.get(iX, iY+1, iZ).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(step, 2*step, step) + vertex);

                    if(!fluid.get(iX, iY, iZ+1).getDynamics().isBoundary())
                        vertices.push_back(Array<T,3>(step, step, 2*step) + vertex);

                    cellId +=1;
            	    for (plint indx=0; indx<vertices.size(); ++indx) {
            	    	boundaryParticleField.addParticle(
            	    			boundaryParticleField.getBoundingBox(),
            	    			new SurfaceParticle3D(vertices[indx], cellId, indx));
            	    }
                    
            	} */
            }
        }
    }
//    std::vector<Particle3D<T,Descriptor>*> particles;
//    boundaryParticleField.findParticles(domain, particles);
//	std::cout << "(PositionBoundaryParticles) pid:" << id << " Number of Particles " << particles.size() << std::endl;

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




/* ******** DeleteParticles3D *********************************** */
template<typename T, template<typename U> class Descriptor>
DeleteParticles3D<T,Descriptor>::DeleteParticles3D () { };

template<typename T, template<typename U> class Descriptor>
void DeleteParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()>0 );
    ParticleField3D<T,Descriptor>& boundaryParticleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    boundaryParticleField.removeParticles(domain); // CAUTION: Domains might not overlap-- palabos deletes by what is contained in the cell-list and not by coordinate.
}

template<typename T, template<typename U> class Descriptor>
DeleteParticles3D<T,Descriptor>* DeleteParticles3D<T,Descriptor>::clone() const {
    return new DeleteParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void DeleteParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field
    modified[1] = modif::nothing; // Fluid field
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DeleteParticles3D<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}


#endif  // FICSION_HELPER_FUNCTIONALS_HH

