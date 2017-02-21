#ifndef SURFACE_PARTICLE_3D_HH
#define SURFACE_PARTICLE_3D_HH

#include "surfaceParticle3D.h"
#include <limits>       // std::numeric_limits



namespace plb {

/* *************** class ImmersedCellParticle3D ************************************ */


SurfaceParticle3D::SurfaceParticle3D (Array<double,3> const& position, plint cellId_, plint vertexId_,pluint celltype_)
    : Particle3D<double,DESCRIPTOR>(cellId_, position), // The cellId initializor does nothing
      pbcPosition(position),
      v(),
      force(),
      vPrevious(),
      cellId(cellId_), 
      vertexId(vertexId_),
      celltype(celltype_),
      rank(getMpiProcessor())
{

}


void SurfaceParticle3D::advance(){

    /* scheme:
     *  1: Euler
     *  2: Adams-Bashforth
     */
    #if HEMOCELL_MATERIAL_INTEGRATION == 1
            Array<double,3> dxyz = v;
        	this->getPosition() += dxyz;         
        	pbcPosition += dxyz;                 

    #elif HEMOCELL_MATERIAL_INTEGRATION == 2
            Array<double,3> dxyz = (1.5*v - 0.5*vPrevious);
        	this->getPosition() +=  dxyz;
        	pbcPosition += dxyz;
        	
        	vPrevious = v;  // Store velocity
    #endif
}


void SurfaceParticle3D::reset(Array<double,3> const& position_, Array<double,3> const& velocity_) {
        Particle3D<double,DESCRIPTOR>::reset(position_);
        pbcPosition = position_;
        v = velocity_;
        vPrevious = velocity_;
        resetForces();
        rank = this->getMpiProcessor();
}


void SurfaceParticle3D::resetForces() {
        force.resetToZero();
}


void SurfaceParticle3D::reset(Array<double,3> const& position_)
{
        reset(position_,Array<double,3>(0,0,0));
}

void SurfaceParticle3D::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<double,DESCRIPTOR>::serialize(serializer);
    serializer.addValues<double,3>(v);
    serializer.addValues<double,3>(pbcPosition);
    serializer.addValues<double,3>(force);
    serializer.addValues<double,3>(vPrevious);
    serializer.addValue<int>(rank);
    serializer.addValue<plint>(cellId);
    serializer.addValue<plint>(vertexId);
}

void SurfaceParticle3D::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<double,DESCRIPTOR>::unserialize(unserializer);
    unserializer.readValues<double,3>(v);
    unserializer.readValues<double,3>(pbcPosition);
    unserializer.readValues<double,3>(force);
    unserializer.readValues<double,3>(vPrevious);
    unserializer.readValue<int>(rank);
    unserializer.readValue<plint>(cellId);
    unserializer.readValue<plint>(vertexId);
}


SurfaceParticle3D* SurfaceParticle3D::clone() const {
    return new SurfaceParticle3D(*this);
}

}  // namespace plb

#endif  // SURFACE_PARTICLE_3D_HH
