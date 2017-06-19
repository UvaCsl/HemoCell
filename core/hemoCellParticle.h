#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

class HemoCellParticle;
#include "hemocell_internal.h"


class HemoCellParticle {
  //VARIABLES
  //Store variables in struct for fast serialization
  struct serializeValues_t {
    Array<double,3> v;
    Array<double,3> position;
    Array<double,3> force;
    Array<double,3> force_repulsion;
    plint cellId;
    plint vertexId;
    pluint celltype;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    Array<double,3> vPrevious;
#endif
  };
  
  public:

    serializeValues_t serializeValues = {};
    //Is vector, optimize with array possible
    vector<Cell<double,DESCRIPTOR>*> kernelLocations;
    vector<double>         kernelWeights;
    Array<double,3> & position;// = serializeValues.sv.v;
    Array<double,3> & v;// = serializeValues.sv.v;
    Array<double,3> & force;// = serializeValues.sv.force;
    Array<double,3> force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    Array<double,3> vPrevious;
#endif
    Array<double,3> *force_volume = &force;
    Array<double,3> *force_bending = &force;
    Array<double,3> *force_link = &force;
    Array<double,3> & force_repulsion;// = serializeValues.sv.force_repulsion;
    Array<double,3> *force_area = &force; //Default to pointing to force, if output is desired, it can be stored seperately
    Array<double,3> *force_visc = &force;
    plint tag;
    plint & cellId;// = serializeValues.sv.cellId;
    plint & vertexId;// = serializeValues.sv.vertexId;
    pluint & celltype;// = serializeValues.sv.celltype;
private:
    std::vector<Dot3D> cellPos;
    std::vector<double> weights;

public:
  ~HemoCellParticle(){};
  HemoCellParticle(const HemoCellParticle & copy) :
    position(serializeValues.position), v(serializeValues.v),
    force(serializeValues.force), force_repulsion(serializeValues.force_repulsion),
    cellId(serializeValues.cellId), vertexId(serializeValues.vertexId),
    celltype(serializeValues.celltype)
  {
    kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
    position = copy.position;
    v = copy.v;
    force = copy.force;
    force_total = copy.force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    vPrevious = copy.vPrevious;
#endif
    if (&copy.force == copy.force_volume) {
      force_volume = &force;
      force_bending = &force;
      force_link = &force;
      force_area = &force;
      force_visc = &force;
    } else {
      force_volume = copy.force_volume;
      force_bending = copy.force_bending;
      force_link = copy.force_link;
      force_area = copy.force_area;
      force_visc = copy.force_visc;
    }
    force_repulsion = copy.force_repulsion;
    tag = copy.tag;
    cellId = copy.cellId;
    vertexId = copy.vertexId;
    celltype = copy.celltype;
    cellPos = copy.cellPos;
    weights = copy.weights;
  }
    HemoCellParticle() :
    position(serializeValues.position), v(serializeValues.v),
    force(serializeValues.force), force_repulsion(serializeValues.force_repulsion),
    cellId(serializeValues.cellId), vertexId(serializeValues.vertexId),
    celltype(serializeValues.celltype)
    {
      position = {0.0,0.0,0.0};
      v = {0.0,0.0,0.0};
      force = {0.0,0.0,0.0};
#if HEMOCELL_MATERIAL_INTEGRATION == 2
      vPrevious = {0.0,0.0,0.0};
#endif
      tag = -1;
      cellId = 0;
      vertexId = 0;
      celltype = 0;
      force_volume = &force; //These pointers are only changed for nice outputs
      force_area = &force; //These pointers are only changed for nice outputs
      force_link = &force; //These pointers are only changed for nice outputs
      force_bending = &force; //These pointers are only changed for nice outputs
      force_visc = &force;
    }
    HemoCellParticle (Array<double,3> position_, plint cellId_, plint vertexId_,pluint celltype_) :
    position(serializeValues.position), v(serializeValues.v),
    force(serializeValues.force), force_repulsion(serializeValues.force_repulsion),
    cellId(serializeValues.cellId), vertexId(serializeValues.vertexId),
    celltype(serializeValues.celltype)
{
      force = {0.,0.,0.};
      force_repulsion = {0.,0.,0.};
      cellId = cellId_;
      vertexId = vertexId_;
      celltype=celltype_;
      position = position_;
      v = {0.0,0.0,0.0};
#if HEMOCELL_MATERIAL_INTEGRATION == 2
      vPrevious = {0.0,0.0,0.0};
#endif
      tag = -1;
      force_volume = &force; //These pointers are only changed for nice outputs
      force_area = &force; //These pointers are only changed for nice outputs
      force_link = &force; //These pointers are only changed for nice outputs
      force_bending = &force; //These pointers are only changed for nice outputs
      force_visc = &force;
    }

    HemoCellParticle & operator =(const HemoCellParticle & copy) {
       kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
    position = copy.position;
    v = copy.v;
    force = copy.force;
    force_total = copy.force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    vPrevious = copy.vPrevious;
#endif
    if (&copy.force == copy.force_volume) {
      force_volume = &force;
      force_bending = &force;
      force_link = &force;
      force_area = &force;
      force_visc = &force;
    } else {
      force_volume = copy.force_volume;
      force_bending = copy.force_bending;
      force_link = copy.force_link;
      force_area = copy.force_area;
      force_visc = copy.force_visc;
    }
    force_repulsion = copy.force_repulsion;
    tag = copy.tag;
    cellId = copy.cellId;
    vertexId = copy.vertexId;
    celltype = copy.celltype;
    cellPos = copy.cellPos;
    weights = copy.weights;
      return *this;
    }
    HemoCellParticle* clone() const {
        HemoCellParticle* sparticle = new HemoCellParticle(*this);
        sparticle->force_volume = &sparticle->force;
        sparticle->force_bending = &sparticle->force;
        sparticle->force_link = &sparticle->force;
        sparticle->force_area = &sparticle->force;
        sparticle->force_visc = &sparticle->force;
        return sparticle;
    }

    inline void repoint_force_vectors() {
      force_volume = &force;
      force_bending = &force;
      force_link = &force;
      force_area = &force;
      force_visc = &force;
    }

    /// Implements Euler integration with velocity alone.
    void advance() {

        /* scheme:
         *  1: Euler 
         *  2: Adams-Bashforth
         */
        #if HEMOCELL_MATERIAL_INTEGRATION == 1
              position += v;    
             //vPrevious = v;    // Store previous velocity for viscosity terms     

        #elif HEMOCELL_MATERIAL_INTEGRATION == 2
              Array<double,3> dxyz = (1.5*v - 0.5*vPrevious);
              position +=  dxyz;
              vPrevious = v;  // Store velocity
        #endif
        //v = {0.0,0.0,0.0};
    }
    void serialize(HierarchicSerializer& serializer) const 
    {
        serializer.addValue<serializeValues_t>(serializeValues);
    }
    void unserialize(HierarchicUnserializer& unserializer) 
    {
        unserializer.readValue<serializeValues_t>(serializeValues);
        //These pointers are only changed for nice outputs
        force_volume = &force; 
        force_area = &force; 
        force_link = &force;
        force_bending = &force;
        force_visc = &force;
        
    }

    static int id;
    inline int getId() const {return id;}

    //inline Array<double,3> & getPosition() {return position;} 
    inline plint getTag() { return tag; } //TODO remove for direct access
    inline void setTag(plint tag_) { tag = tag_; }

};

//TODO better way to override this function
inline void serialize(HemoCellParticle& particle, vector<char>& data) {
  HierarchicSerializer serializer(data, HemoCellParticle::id);
  particle.serialize(serializer);
}


#endif  // SURFACE_PARTICLE_3D_H

