#include "geometryUtils.h"

double getAngleBetweenFaces(const Array<double,3> n1, const Array<double,3> n2, const Array<double,3> edge) {
	Array<T,3> cross; 
	crossProduct (n1, n2, cross);
	return std::atan2(dot(cross, edge), dot(n1, n2));
}
