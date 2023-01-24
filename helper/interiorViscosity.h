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
#ifndef HEMO_INTERIOR_VISCOSITY_H
#define HEMO_INTERIOR_VISCOSITY_H

#include "hemoCellParticleField.h"
#include "multiBlock/multiDataField3D.h"

namespace hemo {
  class InteriorViscosityHelper {
  public:
     static InteriorViscosityHelper& get(HemoCellFields & cellFields) {

      static InteriorViscosityHelper instance(cellFields);
      return instance;
    }
      
    void checkpoint();
    static void restore(HemoCellFields & cellFields);
    
    //Called within functional, from particlefield
    void add(HemoCellParticleField & pf, const Dot3D & bindingSite, T tau);
    void add(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites, T tau);
    void remove(HemoCellParticleField & pf, const Dot3D & bindingSite);
    void remove(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites);
    void empty(HemoCellParticleField & pf);
    
  private:
    HemoCellFields & cellFields;

    plb::MultiScalarField3D<T> *multiInteriorViscosityField = nullptr,
     *preinlet_multiInteriorViscosityField = nullptr, *domain_multiInteriorViscosityField = nullptr;
    
    InteriorViscosityHelper(HemoCellFields & cellFields);
    ~InteriorViscosityHelper();
    
    void refillBindingSites();
    
    //Singleton Behaviour
  public:
    InteriorViscosityHelper(InteriorViscosityHelper const&) = delete;
    void operator=(InteriorViscosityHelper const&) = delete;  
  };
}
#endif