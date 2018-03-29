/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SURFACE_PARTICLE_3D_H
#define SURFACE_PARTICLE_3D_H

class HemoCellParticle;
#include "hemocell_internal.h"

#ifndef PARTICLE_ID
#define PARTICLE_ID 0
#endif

class HemoCellParticle {
  //VARIABLES
  //Store variables in struct for fast serialization
  struct serializeValues_t {
    hemo::Array<T,3> v;
    hemo::Array<T,3> position;
    hemo::Array<T,3> force;
    hemo::Array<T,3> force_repulsion;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    hemo::Array<T,3> vPrevious;
#endif
    plint cellId;
    plint vertexId;
    pluint celltype;

    bool fromPreInlet;

  };
  
  public:

    serializeValues_t sv;
    //Is vector, optimize with hemo::Array possible
    vector<Cell<T,DESCRIPTOR>*> kernelLocations;
    vector<T>         kernelWeights;
    //hemo::Array<T,3> & position = sv.position;
    //hemo::Array<T,3> & v = sv.v;
    //hemo::Array<T,3> & force = sv.force;
    hemo::Array<T,3> force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    hemo::Array<T,3> vPrevious;
#endif
    hemo::Array<T,3> *force_volume;// = &sv.force;
    hemo::Array<T,3> *force_bending;// = &sv.force;
    hemo::Array<T,3> *force_link;// = &sv.force;
    //hemo::Array<T,3> & force_repulsion = sv.force_repulsion;
    hemo::Array<T,3> *force_area;// = &sv.force; //Default to pointing to force, if output is desired, it can be stored seperately
    hemo::Array<T,3> *force_visc;// = &sv.force;
    hemo::Array<T,3> *force_inner_link;// = &sv.force;
    bool & fromPreInlet = sv.fromPreInlet;
    plint tag;
    //plint & cellId = sv.cellId;
    //plint & vertexId = sv.vertexId;
    //pluint & celltype = sv.celltype;
private:
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;

public:
  ~HemoCellParticle(){};
  HemoCellParticle(const HemoCellParticle & copy) /*:
    position(sv.position), v(sv.v),
    force(sv.force), force_repulsion(sv.force_repulsion),
    fromPreInlet(sv.fromPreInlet),
    cellId(sv.cellId), vertexId(sv.vertexId),
    celltype(sv.celltype)*/
  {
    //sv = copy.sv;
    memcpy(this,&copy,sizeof(HemoCellParticle));
    //kernelLocations = copy.kernelLocations;
    //kernelWeights = copy.kernelWeights;
    //force_total = copy.force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    vPrevious = copy.vPrevious;
#endif
    if (&copy.sv.force == copy.force_volume) {
      force_volume = &sv.force;
      force_bending = &sv.force;
      force_link = &sv.force;
      force_area = &sv.force;
      force_visc = &sv.force;
      force_inner_link = &sv.force;
    } else {
      force_volume = copy.force_volume;
      force_bending = copy.force_bending;
      force_link = copy.force_link;
      force_area = copy.force_area;
      force_visc = copy.force_visc;
      force_inner_link = copy.force_inner_link;
    }
    //force_repulsion = copy.force_repulsion;
    //fromPreInlet = copy.fromPreInlet;
    tag = copy.tag;
    //cellId = copy.cellId;
    //vertexId = copy.vertexId;
    //celltype = copy.celltype;
    cellPos = copy.cellPos;
    weights = copy.weights;
  }
    HemoCellParticle() /*:
    position(sv.position), v(sv.v),
    force(sv.force), force_repulsion(sv.force_repulsion),
    fromPreInlet(sv.fromPreInlet),
    cellId(sv.cellId), vertexId(sv.vertexId),
    celltype(sv.celltype)*/
    {
      memset(&sv,0,sizeof(serializeValues_t));
      //sv.position = {0.0,0.0,0.0};
      //sv.v = {0.0,0.0,0.0};
      //force = {0.0,0.0,0.0};
#if HEMOCELL_MATERIAL_INTEGRATION == 2
      vPrevious = {0.0,0.0,0.0};
#endif
      //fromPreInlet = false;
      tag = -1;
      //cellId = 0;
      //vertexId = 0;
      //celltype = 0;
      force_volume = &sv.force; //These pointers are only changed for nice outputs
      force_area = &sv.force; //These pointers are only changed for nice outputs
      force_link = &sv.force; //These pointers are only changed for nice outputs
      force_bending = &sv.force; //These pointers are only changed for nice outputs
      force_visc = &sv.force;
      force_inner_link = &sv.force;
    }
    HemoCellParticle (hemo::Array<T,3> position_, plint cellId_, plint vertexId_,pluint celltype_)/* :
    position(sv.position), v(sv.v),
    force(sv.force), force_repulsion(sv.force_repulsion),
    fromPreInlet(sv.fromPreInlet),
    cellId(sv.cellId), vertexId(sv.vertexId),
    celltype(sv.celltype)*/
{
      memset(&sv,0,sizeof(serializeValues_t));
      //force = {0.,0.,0.};
      //force_repulsion = {0.,0.,0.};
      sv.cellId = cellId_;
      sv.vertexId = vertexId_;
      sv.celltype=celltype_;
      sv.position = position_;
      //v = {0.0,0.0,0.0};
#if HEMOCELL_MATERIAL_INTEGRATION == 2
      vPrevious = {0.0,0.0,0.0};
#endif
      //fromPreInlet = false;
      tag = -1;
      force_volume = &sv.force; //These pointers are only changed for nice outputs
      force_area = &sv.force; //These pointers are only changed for nice outputs
      force_link = &sv.force; //These pointers are only changed for nice outputs
      force_bending = &sv.force; //These pointers are only changed for nice outputs
      force_visc = &sv.force;
      force_inner_link = &sv.force;
    }

    HemoCellParticle & operator =(const HemoCellParticle & copy) {
      //sv = copy.sv;
      memcpy(this,&copy,sizeof(HemoCellParticle));
       kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
    //position = copy.position;
    //v = copy.v;
    //force = copy.force;
    force_total = copy.force_total;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    vPrevious = copy.vPrevious;
#endif
    if (&copy.sv.force == copy.force_volume) {
      force_volume = &sv.force;
      force_bending = &sv.force;
      force_link = &sv.force;
      force_area = &sv.force;
      force_visc = &sv.force;
      force_inner_link = &sv.force;
    } else {
      force_volume = copy.force_volume;
      force_bending = copy.force_bending;
      force_link = copy.force_link;
      force_area = copy.force_area;
      force_visc = copy.force_visc;
      force_inner_link = copy.force_inner_link;
    }
    //force_repulsion = copy.force_repulsion;
    //fromPreInlet = copy.fromPreInlet;
    tag = copy.tag;
    //cellId = copy.cellId;
    //vertexId = copy.vertexId;
    //celltype = copy.celltype;
    cellPos = copy.cellPos;
    weights = copy.weights;
      return *this;
    }
    HemoCellParticle* clone() const {
        HemoCellParticle* sparticle = new HemoCellParticle(*this);
        sparticle->force_volume = &sparticle->sv.force;
        sparticle->force_bending = &sparticle->sv.force;
        sparticle->force_link = &sparticle->sv.force;
        sparticle->force_area = &sparticle->sv.force;
        sparticle->force_visc = &sparticle->sv.force;
        sparticle->force_inner_link = &sparticle->sv.force;
        return sparticle;
    }
    
    inline void repoint_force_vectors() {
      force_volume = &sv.force;
      force_bending = &sv.force;
      force_link = &sv.force;
      force_area = &sv.force;
      force_visc = &sv.force;
      force_inner_link = &sv.force;
    }

    /// Implements Euler integration with velocity alone.
    void advance() {

        /* scheme:
         *  1: Euler 
         *  2: Adams-Bashforth
         */
        #if HEMOCELL_MATERIAL_INTEGRATION == 1
              sv.position += sv.v;    
             //vPrevious = v;    // Store previous velocity for viscosity terms     

        #elif HEMOCELL_MATERIAL_INTEGRATION == 2
              hemo::Array<T,3> dxyz = (1.5*v - 0.5*vPrevious);
              position +=  dxyz;
              vPrevious = v;  // Store velocity
        #endif
        //v = {0.0,0.0,0.0};
    }
    void serialize(HierarchicSerializer& serializer) const 
    {
        serializer.addValue<serializeValues_t>(sv);
    }
    void unserialize(HierarchicUnserializer& unserializer) 
    {
        unserializer.readValue<serializeValues_t>(sv);
        //These pointers are only changed for nice outputs
        force_volume = &sv.force; 
        force_area = &sv.force; 
        force_link = &sv.force;
        force_bending = &sv.force;
        force_visc = &sv.force;
        force_inner_link = &sv.force;
        
    }

    inline int getId() const {return PARTICLE_ID;}
    //inline hemo::Array<T,3> & getPosition() {return position;} 
    inline plint getTag() { return tag; } //TODO remove for direct access
    inline void setTag(plint tag_) { tag = tag_; }

};

//TODO better way to override this function
inline void serialize(HemoCellParticle& particle, vector<char>& data) {
  HierarchicSerializer serializer(data,PARTICLE_ID);
  particle.serialize(serializer);
}

#endif  // SURFACE_PARTICLE_3D_H

