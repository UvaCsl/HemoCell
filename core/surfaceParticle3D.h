#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "particleHelperFunctions.h"
#include <vector>


namespace plb {

class SurfaceParticle3D : public Particle3D<double,DESCRIPTOR> {
public:
    SurfaceParticle3D(Array<double,3> const& position, plint cellId_ = -1, plint vertexId_ = 0, pluint celltype_=0);
    SurfaceParticle3D* clone() const override;

    void velocityToParticle(TensorField3D<double,3>& velocityField, double scaling=1.) override { }
    void velocityToParticle(NTensorField3D<double>& velocityField, double scaling=1.) override { }
    void rhoBarJtoParticle(NTensorField3D<double>& rhoBarJfield, bool velIsJ, double scaling=1.) override { }
    void fluidToParticle(BlockLattice3D<double,DESCRIPTOR>& fluid, double scaling=1.) override { }

    /// Implements Euler integration with velocity alone.
    void advance() override;
    void serialize(HierarchicSerializer& serializer) const override;
    void unserialize(HierarchicUnserializer& unserializer) override;

    void reset(Array<double,3> const& position, Array<double,3> const& velocity_);
    void reset(Array<double,3> const& position) override;
    void resetForces();
    
    int getId() const {return 0;}

    Array<double,3> const& get_v() const { return v; }
    Array<double,3> const& getVelocity() const { return get_v(); }
    Array<double,3> const& get_pbcPosition() const { return pbcPosition; }
    Array<double,3> const& get_vPrevious() const { return vPrevious; }
    Array<double,3> const& get_force() const { return force; }
    plint const& get_cellId() const { return cellId; }
    plint const& getVertexId() const { return vertexId; }
    // Difference between getVelocity and get_v:
    // get_v holds the actual interpolated velocity, while
    // getVelocity holds the velocity according to the
    // update of the position particle "method advance()".
    Array<double,3>& get_v() { return v; }
    //Array<T,3>& get_pbcPosition() { return pbcPosition; }
    //Array<T,3>& get_vPrevious() { return vPrevious; }
    //Array<T,3>& get_a() { return a;double}
    Array<double,3>& get_force() { return force; }
    //plint& get_cellId() { return cellId; }
    //plint& getVertexId() { return vertexId; }
    //plint& get_processor() { return processor; }
    int getMpiProcessor() { MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank;}
private:
    Array<double,3> pbcPosition, v, force, vPrevious;
    plint cellId;
    plint vertexId;
    pluint celltype;
    int rank;
public:
    std::vector<Dot3D> & getIBMcoordinates() { return cellPos; }
    std::vector<double> & getIBMweights() { return weights; }
private:
    std::vector<Dot3D> cellPos;
    std::vector<double> weights;
};

}  // namespace plb

#include "surfaceParticle3D.hh"

#endif  // SURFACE_PARTICLE_3D_H

