#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

class SurfaceParticle3D;
#include "hemocell_internal.h"

class SurfaceParticle3D : public Particle3D<double,DESCRIPTOR> {
public:
    SurfaceParticle3D(){
      force_volume = &force; //These pointers are only changed for nice outputs
      force_area = &force; //These pointers are only changed for nice outputs
      force_inplane = &force; //These pointers are only changed for nice outputs
      force_bending = &force; //These pointers are only changed for nice outputs
    }
    SurfaceParticle3D (Array<double,3> const& position, plint cellId_, plint vertexId_,pluint celltype_)
        : Particle3D<double,DESCRIPTOR>(-1, position), // The cellId initializor does nothing
          v(),
          force(),
          vPrevious(),
          cellId(cellId_), 
          vertexId(vertexId_),
          celltype(celltype_)
    {
      force_volume = &force; //These pointers are only changed for nice outputs
      force_area = &force; //These pointers are only changed for nice outputs
      force_inplane = &force; //These pointers are only changed for nice outputs
      force_bending = &force; //These pointers are only changed for nice outputs
    }
    SurfaceParticle3D* clone() const override {
        SurfaceParticle3D* sparticle = new SurfaceParticle3D(*this);
        sparticle->force_volume = &sparticle->force;
        sparticle->force_bending = &sparticle->force;
        sparticle->force_inplane = &sparticle->force;
        sparticle->force_area = &sparticle->force;
        return sparticle;
    }

    void velocityToParticle(TensorField3D<double,3>& velocityField, double scaling=1.) override {}
    void velocityToParticle(NTensorField3D<double>& velocityField, double scaling=1.) override {}
    void rhoBarJtoParticle(NTensorField3D<double>& rhoBarJfield, bool velIsJ, double scaling=1.) override {}
    void fluidToParticle(BlockLattice3D<double,DESCRIPTOR>& fluid, double scaling=1.) override {}

    /// Implements Euler integration with velocity alone.
    void advance() override{

        /* scheme:
         *  1: Euler
         *  2: Adams-Bashforth
         */
        #if HEMOCELL_MATERIAL_INTEGRATION == 1
              this->getPosition() += v;         

        #elif HEMOCELL_MATERIAL_INTEGRATION == 2
                Array<double,3> dxyz = (1.5*v - 0.5*vPrevious);
              this->getPosition() +=  dxyz;
              vPrevious = v;  // Store velocity
        #endif
        v = {0.0,0.0,0.0};
    }
    void serialize(HierarchicSerializer& serializer) const override
    {
        Particle3D<double,DESCRIPTOR>::serialize(serializer);
        serializer.addValues<double,3>(v);
        serializer.addValues<double,3>(force);
        serializer.addValues<double,3>(vPrevious);
        serializer.addValue<plint>(cellId);
        serializer.addValue<plint>(vertexId);
        serializer.addValue<pluint>(celltype);
    }
    void unserialize(HierarchicUnserializer& unserializer) override
    {
        Particle3D<double,DESCRIPTOR>::unserialize(unserializer);
        unserializer.readValues<double,3>(v);
        unserializer.readValues<double,3>(force);
        unserializer.readValues<double,3>(vPrevious);
        unserializer.readValue<plint>(cellId);
        unserializer.readValue<plint>(vertexId);
        unserializer.readValue<pluint>(celltype);
        //These pointers are only changed for nice outputs
        force_volume = &force; 
        force_area = &force; 
        force_inplane = &force;
        force_bending = &force;
    }

    static int id;
    inline int getId() const override {return id;}

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
    pluint celltype;

private:
    std::vector<Dot3D> cellPos;
    std::vector<double> weights;
};


#endif  // SURFACE_PARTICLE_3D_H

