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
#ifndef HEMO_MOLLER_H
#define HEMO_MOLLER_H

#include "array.h" // Need to make pointers to particle object

namespace hemo {
  int MollerTrumbore(const hemo::Array<double,3> v0, const hemo::Array<double,3> v1,
	  const hemo::Array<double,3> v2, hemo::Array<double, 3> rayVector, 
	  hemo::Array<plint, 3> latticeSite, const double EPSILON); 
}

#endif
