#ifndef HEMOCELL_GEOMETRY_UTILS_H
#define HEMOCELL_GEOMETRY_UTILS_H

#include "hemocell_internal.h"

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
//double getAngleBetweenFaces(const Array<double,3> n1, const Array<double,3> n2, const Array<double,3> edge);
//inline since performance critical
inline double getAngleBetweenFaces(const Array<double,3> n1, const Array<double,3> n2, const Array<double,3> edge) {
	Array<T,3> cross; 
	crossProduct (n1, n2, cross);
	return std::atan2(dot(cross, edge), dot(n1, n2));
}

#endif