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
#ifndef HEMO_BINDING_FIELD_H
#define HEMO_BINDING_FIELD_H

#include "hemoCellParticleField.h"
#include "multiBlock/multiDataField3D.h"

namespace hemo {
  class bindingFieldHelper {
  public:
    static bindingFieldHelper& get(HemoCellFields & cellFields) {
      static bindingFieldHelper instance(&cellFields);
      return instance;
    }
      
    void checkpoint();
    static void restore(HemoCellFields & cellFields);
    
    //Called within functional, from particlefield
    void add(HemoCellParticleField & pf, const Dot3D & bindingSite);
    void add(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites);
    void remove(HemoCellParticleField & pf, const Dot3D & bindingSite);
    void remove(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites);
    
  private:
    HemoCellFields & cellFields;
    plb::MultiScalarField3D<bool> * multiBindingField = 0;
    
    bindingFieldHelper(HemoCellFields * cellFields);
    ~bindingFieldHelper();
    
    void refillBindingSites();
    
    //Singleton Behaviour
  public:
    bindingFieldHelper(bindingFieldHelper const&) = delete;
    void operator=(bindingFieldHelper const&) = delete;  
  };
}
#endif