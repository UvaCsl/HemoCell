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
#ifndef HEMOCELL_GEOMETRY_UTILS_H
#define HEMOCELL_GEOMETRY_UTILS_H

#include "array.h"
#include "constant_defaults.h"

/*
- returns atan2((Va x Vb) . Vn, Va . Vb)
- change Va x Vb to Vb x Va for left handed coordinate system
Explanation:
Va . Vb == |Va| * |Vb| * cos(alpha)    (by definition) 
        == |Va| * |Vb| * cos(beta)     (cos(alpha) == cos(-alpha) == cos(360° - alpha)


Va x Vb == |Va| * |Vb| * sin(alpha) * n1  
    (by definition; n1 is a unit vector perpendicular to Va and Vb with 
     orientation matching the right-hand rule)

Therefore (again assuming Vn is normalized):
   n1 . Vn == 1 when beta < 180
   n1 . Vn == -1 when beta > 180

==>  (Va x Vb) . Vn == |Va| * |Vb| * sin(beta)
==>  tan(beta) = sin(beta) / cos(beta) == ((Va x Vb) . Vn) / (Va . Vb)
*/
inline T getAngleBetweenFaces(const hemo::Array<T,3> n1, const hemo::Array<T,3> n2, const hemo::Array<T,3> edge) {
	hemo::Array<T,3> cross = crossProduct (n1, n2);
	return std::atan2(dot(cross, edge), dot(n1, n2));
};

#endif