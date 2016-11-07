#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "particleHelperFunctions.h"
#include <vector>


namespace plb {

template<typename T, template<typename U> class Descriptor>
class SurfaceParticle3D : public Particle3D<T,Descriptor> {
public:
    SurfaceParticle3D();
    ~SurfaceParticle3D() {
//        std::cout <<" ~ImmersedCellParticle3D() " << global::mpi().getRank() << " vertexId " << this->getVertexId() ;
//        std::cout << " pos (" << this->getPosition()[0] << ", "<< this->getPosition()[1] << ", "<< this->getPosition()[2] << ") " << std::endl;
    };
    /* scheme:
     *  0: Euler
     *  1: Adams-Bashforth
     *  2: Velocity Verlet
     */
    SurfaceParticle3D(SurfaceParticle3D<T,Descriptor> const& rhs);
    SurfaceParticle3D(Array<T,3> const& position, plint cellId_ = -1, plint vertexId_ = 0, T dt_= 1.0);
    SurfaceParticle3D(Array<T,3> const& position,
                          Array<T,3> const& v_, Array<T,3> const& pbcPosition_,
                            Array<T,3> const& a_, Array<T,3> const& force_, Array<T,3> const& vPrevious_,
                            plint cellId_ = -1, plint vertexId_=0, plint scheme_=0, T dt_ = 1);
    virtual SurfaceParticle3D<T,Descriptor>* clone() const;
    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
    virtual void velocityToParticle(NTensorField3D<T>& velocityField, T scaling=1.) { }
    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.) { }
    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
    /// Implements Euler integration with velocity alone.
    virtual void advance();
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    virtual int getId() const;
    virtual void reset(Array<T,3> const& position, Array<T,3> const& velocity_, bool allVariables=false);
    virtual void reset(Array<T,3> const& position);
    virtual void resetForces();
    /// Return the cellId through a generic interface (vector id=0).
    virtual bool getScalar(plint whichScalar, T& scalar) const;
    std::string getScalarName(plint whichScalar) const;
    plint getScalarsNumber() const ;
    /// Return the velocity, acceleration or pbcPosition through a generic interface (vector id=0,1,2).
    virtual bool getVector(plint whichVector, Array<T,3>& vector) const;
    std::string getVectorName(plint whichVector) const;
    plint getVectorsNumber() const ;
    Array<T,3> const& get_v() const { return v; }
    Array<T,3> const& get_pbcPosition() const { return pbcPosition; }
    Array<T,3> const& get_vPrevious() const { return vPrevious; }
    Array<T,3> const& get_a() const { return a; }
    Array<T,3> const& get_force() const { return force; }
    // Difference between getVelocity and get_v:
    // get_v holds the actual interpolated velocity, while
    // getVelocity holds the velocity according to the
    // update of the position particle "method advance()".
    Array<T,3> getVelocity() { return vProgressed; }
    Array<T,3>& get_v() { return v; }
    Array<T,3>& get_pbcPosition() { return pbcPosition; }
    Array<T,3>& get_vPrevious() { return vPrevious; }
    Array<T,3>& get_a() { return a; }
    Array<T,3>& get_force() { return force; }
    T& get_dt() { return dt; }
    void set_dt( T dt_ ) { dt = dt_; }
    plint const& get_cellId() const { return cellId; }
    plint& get_cellId() { return cellId; }
    plint const& getVertexId() const { return vertexId; }
    plint& getVertexId() { return vertexId; }
    plint const& get_processor() const { return processor; }
    plint& get_processor() { return processor; }
    plint getMpiProcessor() {
    	int myrank = 0;
#ifdef PLB_MPI_PARALLEL
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    	return plint(myrank); }
private:
    Array<T,3> v, pbcPosition, a, force, vPrevious, vProgressed;
    T dt=1.0;
    static int id;
private:
    plint processor;
    plint cellId;
    plint vertexId;
public:
    std::vector<Dot3D> & getIBMcoordinates() { return cellPos; }
    std::vector<T> & getIBMweights() { return weights; }
private:
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
private:
    plint scheme;
public:
    plint & get_scheme() { return scheme; }
public:
#ifdef PLB_DEBUG // Less Calculations

    Array<T,3> const& get_f_wlc() const { return f_wlc; }
    Array<T,3> const& get_f_bending() const { return f_bending; }
    Array<T,3> const& get_f_volume() const { return f_volume; }
    Array<T,3> const& get_f_surface() const { return f_surface; }
    Array<T,3> const& get_f_shear() const { return f_shear; }
    Array<T,3> const& get_f_viscosity() const { return f_viscosity; }
    Array<T,3> const& get_f_repulsive() const { return f_repulsive; }
    Array<T,3> const& get_stress() const { return stress; }

    T const& get_E_other() const { return E_other; }
    T const& get_E_inPlane() const { return E_inPlane; }
    T const& get_E_bending() const { return E_bending; }
    T const& get_E_area() const { return E_area; }
    T const& get_E_volume() const { return E_volume; }
    T const& get_E_repulsive() const { return E_repulsive; }

    Array<T,3>& get_f_wlc() { return f_wlc; }
    Array<T,3>& get_f_bending() { return f_bending; }
    Array<T,3>& get_f_volume() { return f_volume; }
    Array<T,3>& get_f_surface() { return f_surface; }
    Array<T,3>& get_f_shear() { return f_shear; }
    Array<T,3>& get_f_viscosity() { return f_viscosity; }
    Array<T,3>& get_f_repulsive() { return f_repulsive; }
    Array<T,3>& get_stress() { return stress; }

    T& get_E_other() { return E_other; }
    T& get_E_inPlane() { return E_inPlane; }
    T& get_E_bending() { return E_bending; }
    T& get_E_area() { return E_area; }
    T& get_E_volume() { return E_volume; }
    T& get_E_repulsive() { return E_repulsive; }

    T const get_E_total() const { return (E_other + E_inPlane + E_bending + E_area + E_volume + E_repulsive);}
    T const get_Energy() const { return (E_other + E_inPlane + E_bending + E_area + E_volume + E_repulsive);}
private:
    Array<T,3> f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity, f_repulsive;
    Array<T,3> stress;
    T E_other, E_inPlane, E_bending, E_area,  E_volume, E_repulsive;
#endif
// TROMBOSIT related variables
public:
    T & getBondTypeSaturation(plint bondType) { return (bondTypeSaturation[bondType]); }
    std::map<plint, T> & getBondTypeSaturation() { return bondTypeSaturation; };
private:
    std::map<plint, T> bondTypeSaturation;
};

template<typename T, template<typename U> class Descriptor>
int SurfaceParticle3D<T,Descriptor>::id = meta::registerGenericParticle3D<T,Descriptor,SurfaceParticle3D<T,Descriptor> >("ImmersedCellParticle3D");

template<typename T, template<typename U> class Descriptor>
SurfaceParticle3D<T,Descriptor>* castParticleToICP3D(Particle3D<T,Descriptor>* particle) {
    return  dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particle);
}

template<typename T, template<typename U> class Descriptor>
SurfaceParticle3D<T,Descriptor>* castParticleToSurfaceParticle3D(Particle3D<T,Descriptor>* particle) {
    return  dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particle);
}


}  // namespace plb

#include "surfaceParticle3D.hh"

#endif  // SURFACE_PARTICLE_3D_H

