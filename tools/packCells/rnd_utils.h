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

#ifndef RND_UTILS_H
#define RND_UTILS_H

#ifndef PI
	#define PI 3.141592653589793
#endif

// Misc utils ---------------------------

inline double max (double a, double b) {
	return (a >= b) ? a : b;
}

inline double min (double a, double b) {
	return (a <= b) ? a : b;
}

inline double sign (double a, double b) {
	return (b >= 0) ? fabs(a) : -fabs(a);
}

// Random number generator -----------------

class Random {
public:
	static void init () {
		srand48(time(NULL));
	}
	static double getRand() {
		return drand48();
	}
	static double getRandLimits(double xFrom, double xTo) {
		return getRand() * (xTo-xFrom) + xFrom;
	}

};

#endif
