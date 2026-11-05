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
namespace hemo {
class HemoCellParticle;
}
#include "helper/array.h"
#include "core/cell.hh"

#include <cstdint>
#include <cstring>
#include <type_traits>
#include <utility>

#ifndef PARTICLE_ID
#define PARTICLE_ID 0
#endif

namespace hemo {

class HemoCellParticle {
public:
  //VARIABLES
  //Store variables in struct for fast serialization
  struct serializeValues_t {
    hemo::Array<T, 3> v;
    hemo::Array<T, 3> position;
    hemo::Array<T, 3> force;
    hemo::Array<T, 3> force_repulsion;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    hemo::Array<T, 3> vPrevious;
#endif
    plint cellId;

    uint16_t vertexId;
    unsigned int restime;

    unsigned char celltype;

#ifdef SOLIDIFY_MECHANICS
    bool solidify;
#endif
  };

  serializeValues_t sv;

  hemo::Array<T, 3> force_total;
  plint tag;
#ifdef INTERIOR_VISCOSITY
  hemo::Array<T, 3> normalDirection;
  //Is vector, optimize with hemo::Array possible
  std::vector<hemo::Array<plint, 3>> kernelCoordinates;
#endif

  std::vector<plb::Cell<T, DESCRIPTOR> *> kernelLocations;
  std::vector<T> kernelWeights;

  hemo::Array<T, 3> * force_volume = &sv.force;
  hemo::Array<T, 3> * force_bending = &sv.force;
  hemo::Array<T, 3> * force_link = &sv.force;
  hemo::Array<T, 3> * force_area =
      &sv.force; //Default to pointing to force, if output is desired, it can be stored seperately
  hemo::Array<T, 3> * force_visc = &sv.force;
  hemo::Array<T, 3> * force_inner_link = &sv.force;


private:
  /// The force-component pointers either alias this particle's own sv.force
  /// (normal operation) or point to separately allocated output arrays owned
  /// by the particle field (between separateForceVectors() and
  /// unifyForceVectors()). Copies keep that state: aliased pointers are
  /// re-aimed at the copy's own sv.force, separate arrays are shared.
  inline void adoptForcePointers(const HemoCellParticle & from) {
    if (from.force_volume == &from.sv.force) {
      repoint_force_vectors();
    } else {
      force_volume = from.force_volume;
      force_bending = from.force_bending;
      force_link = from.force_link;
      force_area = from.force_area;
      force_visc = from.force_visc;
      force_inner_link = from.force_inner_link;
    }
  }

public:
  ~HemoCellParticle(){};
  HemoCellParticle(const HemoCellParticle & copy)
      : sv(copy.sv), force_total(copy.force_total), tag(copy.tag),
#ifdef INTERIOR_VISCOSITY
        normalDirection(copy.normalDirection), kernelCoordinates(copy.kernelCoordinates),
#endif
        kernelLocations(copy.kernelLocations), kernelWeights(copy.kernelWeights) {
    adoptForcePointers(copy);
  }

  /// Move: steals the kernel vectors instead of copying them (used by
  /// std::vector reallocation and swap-removal of particles).
  HemoCellParticle(HemoCellParticle && other) noexcept
      : sv(other.sv), force_total(other.force_total), tag(other.tag),
#ifdef INTERIOR_VISCOSITY
        normalDirection(other.normalDirection), kernelCoordinates(std::move(other.kernelCoordinates)),
#endif
        kernelLocations(std::move(other.kernelLocations)), kernelWeights(std::move(other.kernelWeights)) {
    adoptForcePointers(other);
  }

  HemoCellParticle(hemo::Array<T, 3> position_, plint cellId_, plint vertexId_, pluint celltype_) {
    // Zero the padding bytes as well: sv is written as raw bytes to
    // checkpoints, which are then the same for identical runs.
    static_assert(std::is_trivially_copyable<serializeValues_t>::value, "sv is transferred as raw bytes");
    std::memset(static_cast<void *>(&sv), 0, sizeof(sv));
    sv.v = {0., 0., 0.};
    sv.position = position_;
    sv.force = {0., 0., 0.};
    sv.force_repulsion = {0., 0., 0.};
    sv.cellId = cellId_;
    sv.vertexId = vertexId_;
    sv.celltype = celltype_;
    sv.restime = 0;
#if HEMOCELL_MATERIAL_INTEGRATION == 2
    sv.vPrevious = {0., 0., 0.};
#endif
#ifdef SOLIDIFY_MECHANICS
    sv.solidify = false;
#endif
    force_total = {0., 0., 0.};
#ifdef INTERIOR_VISCOSITY
    normalDirection = {0., 0., 0.};
#endif
    tag = -1;

    if (vertexId_ > UINT16_MAX) {
      std::cerr
          << "(HemoCellParticle) Trying to add more vertexes to a single cell than UINT16_MAX, consider converting vertexid to a long int"
          << std::endl;
      exit(1);
    }
  }

  HemoCellParticle(const serializeValues_t & sv_) {
    sv = sv_;
    force_total = {0., 0., 0.};
    tag = -1;
  }

  HemoCellParticle & operator=(const HemoCellParticle & copy) {
    if (this == &copy) { return *this; }
    sv = copy.sv;
    force_total = copy.force_total;
    tag = copy.tag;
    kernelLocations = copy.kernelLocations;
    kernelWeights = copy.kernelWeights;
#ifdef INTERIOR_VISCOSITY
    normalDirection = copy.normalDirection;
    kernelCoordinates = copy.kernelCoordinates;
#endif
    adoptForcePointers(copy);
    return *this;
  }

  HemoCellParticle & operator=(HemoCellParticle && other) noexcept {
    if (this == &other) { return *this; }
    sv = other.sv;
    force_total = other.force_total;
    tag = other.tag;
    kernelLocations = std::move(other.kernelLocations);
    kernelWeights = std::move(other.kernelWeights);
#ifdef INTERIOR_VISCOSITY
    normalDirection = other.normalDirection;
    kernelCoordinates = std::move(other.kernelCoordinates);
#endif
    adoptForcePointers(other);
    return *this;
  }

  HemoCellParticle * clone() const {
    HemoCellParticle * sparticle = new HemoCellParticle(*this);
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
    // vPrevious starts at zero, so the first step uses 1.5*v.
    const hemo::Array<T, 3> dxyz = (1.5 * sv.v - 0.5 * sv.vPrevious);
    sv.position += dxyz;
    sv.vPrevious = sv.v; // Store velocity
#endif
    //v = {0.0,0.0,0.0};
  }

  inline int getId() const { return PARTICLE_ID; }
  inline plint getTag() { return tag; } //TODO remove for direct access
  inline void setTag(plint tag_) { tag = tag_; }
};

}
#endif // SURFACE_PARTICLE_3D_H
