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
#include "helper/array.h"
#include "core/cell.hh"

#include <cstdint> 

#ifndef PARTICLE_ID
#define PARTICLE_ID 0
#endif

namespace hemo {

class HemoCellParticle {
public:

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
    
    uint16_t vertexId;
    
    unsigned char celltype;

    bool fromPreInlet;

#ifdef SOLIDIFY_MECHANICS
    bool solidify;
#endif
  };
  
  serializeValues_t sv;

  hemo::Array<T,3> force_total;
  plint tag;
  //Is vector, optimize with hemo::Array possible
  std::vector<plb::Cell<T,DESCRIPTOR>*> kernelLocations;
  std::vector<T>         kernelWeights;

  hemo::Array<T,3> *force_volume = &sv.force;
  hemo::Array<T,3> *force_bending = &sv.force;
  hemo::Array<T,3> *force_link = &sv.force;
  hemo::Array<T,3> *force_area = &sv.force; //Default to pointing to force, if output is desired, it can be stored seperately
  hemo::Array<T,3> *force_visc = &sv.force;
  hemo::Array<T,3> *force_inner_link = &sv.force;


public:
  ~HemoCellParticle(){};
  HemoCellParticle(const HemoCellParticle & copy) {
    sv = copy.sv;
    force_total = copy.force_total;
    tag = copy.tag;
    kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
    
    if (!(&copy.sv.force == copy.force_volume)) {
      force_volume = copy.force_volume;
      force_bending = copy.force_bending;
      force_link = copy.force_link;
      force_area = copy.force_area;
      force_visc = copy.force_visc;
      force_inner_link = copy.force_inner_link;
    }
  }
  HemoCellParticle() {
    sv = {};
    force_total = {0.,0.,0.};
    tag = -1;
  }
  
  HemoCellParticle (hemo::Array<T,3> position_, plint cellId_, plint vertexId_,pluint celltype_) {
    sv.v = {0.,0.,0.};
    sv.position = position_;
    sv.force = {0.,0.,0.};
    sv.force_repulsion = {0.,0.,0.};
    sv.cellId = cellId_;
    sv.vertexId = vertexId_;
    sv.celltype=celltype_;
#ifdef PREINLET_MECHNICS
    sv.fromPreInlet = false;
#endif
#ifdef SOLIDIFY_MECHANICS
    sv.solidify = false;
#endif
    force_total = {0.,0.,0.};
    tag = -1;
    
    if (vertexId_ > UINT16_MAX) {
      std::cerr << "(HemoCellParticle) Trying to add more vertexes to a single cell than UINT16_MAX, consider converting vertexid to a long int" << std::endl;
      exit(1);
    }
  }
  
  HemoCellParticle (const serializeValues_t & sv_) {
    sv = sv_;
    force_total = {0.,0.,0.};
    tag = -1;
  }

  HemoCellParticle & operator =(const HemoCellParticle & copy) {
    sv = copy.sv;
    kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
    
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
   
    tag = copy.tag;
    sv.cellId = copy.sv.cellId;
    sv.vertexId = copy.sv.vertexId;
    sv.celltype = copy.sv.celltype;

    return *this;
  }
  
  HemoCellParticle* clone() const {
    HemoCellParticle* sparticle = new HemoCellParticle(*this);
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

        #elif HEMOCELL_MATERIAL_INTEGRATION == 2
              hemo::Array<T,3> dxyz = (1.5*v - 0.5*vPrevious);
              position +=  dxyz;
              vPrevious = v;  // Store velocity
        #endif
        //v = {0.0,0.0,0.0};
    }
    
    inline int getId() const {return PARTICLE_ID;}
    inline plint getTag() { return tag; } //TODO remove for direct access
    inline void setTag(plint tag_) { tag = tag_; }

};

}
#endif  // SURFACE_PARTICLE_3D_H

