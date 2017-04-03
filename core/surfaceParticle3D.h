#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

class SurfaceParticle3D;
#include "hemocell_internal.h"

class SurfaceParticle3D : public Particle3D<double,DESCRIPTOR> {
public:
    SurfaceParticle3D();
    SurfaceParticle3D(Array<double,3> const& position, plint cellId_ = -1, plint vertexId_ = 0, pluint celltype_=0);
    SurfaceParticle3D* clone() const override;

    void velocityToParticle(TensorField3D<double,3>& velocityField, double scaling=1.) override;
    void velocityToParticle(NTensorField3D<double>& velocityField, double scaling=1.) override;
    void rhoBarJtoParticle(NTensorField3D<double>& rhoBarJfield, bool velIsJ, double scaling=1.) override;
    void fluidToParticle(BlockLattice3D<double,DESCRIPTOR>& fluid, double scaling=1.) override;

    /// Implements Euler integration with velocity alone.
    void advance() override;
    void serialize(HierarchicSerializer& serializer) const override;
    void unserialize(HierarchicUnserializer& unserializer) override;

    void reset(Array<double,3> const& position, Array<double,3> const& velocity_);
    void reset(Array<double,3> const& position) override;
    void resetForces();

    static int id;
    
    int getId() const;

    Array<double,3> const& get_v() const;
    Array<double,3> const& getVelocity() const; 
    Array<double,3> const& get_vPrevious() const;
    Array<double,3> const& get_force() const;
    plint const& get_cellId() const;
    pluint const& get_celltype() const;
    plint const& getVertexId() const;
    // Difference between getVelocity and get_v:
    // get_v holds the actual interpolated velocity, while
    // getVelocity holds the velocity according to the
    // update of the position particle "method advance()".
    Array<double,3>& get_v();
    //Array<T,3>& get_pbcPosition() { return pbcPosition; }
    //Array<T,3>& get_vPrevious() { return vPrevious; }
    //Array<T,3>& get_a() { return a;double}
    Array<double,3>& get_force();
    //plint& get_cellId() { return cellId; }
    //plint& getVertexId() { return vertexId; }
    //plint& get_processor() { return processor; }
    int getMpiProcessor();

    //Is vector, optimize with array possible
    vector<Cell<double,DESCRIPTOR>*> kernelLocations;
    vector<double>         kernelWeights;
    Array<plint,3> grid_pos;
    Array<double,3> v;
    Array<double,3> force, force_total, vPrevious;
    Array<double,3> *force_volume = &force;
    Array<double,3> *force_bending = &force;
    Array<double,3> *force_inplane = &force;
    Array<double,3> *force_area = &force;; //Default to pointing to force, if output is desired, it can be stored seperately
public:
    plint cellId;
    plint vertexId;
private:
    pluint celltype;
    int rank;

//TODO, remove before production, now to let legacy code compile
public:
   // std::vector<Dot3D> & getIBMcoordinates() { return cellPos; }
  //  std::vector<double> & getIBMweights() { return weights; }
private:
    std::vector<Dot3D> cellPos;
    std::vector<double> weights;
};


#endif  // SURFACE_PARTICLE_3D_H

