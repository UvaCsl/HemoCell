#ifndef RND_UTILS_H
#define RND_UTILS_H

const double PI = 3.141592653589793;

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
	static bool flag;
	static double r1, r2;
public:
	static void init () {
		long seed = 0;
		time(&seed);
		srand48(seed);
		flag = false;
	}
	static double ran() {
		return drand48();
	}
	static double ran1(double xl, double xu) {
		return ran() * (xu-xl) + xl;
	}
	static double rand_normal() {
		if (flag) {
			flag = false;
			return r2;
		}
		double rbase = sqrt(-2.0 * log(ran()));
		double rt = 2 * PI * ran();
		r1 = rbase * cos (rt);
		r2 = rbase * sin (rt);
		flag = true;
		return r1;
	}
	static double randiam () {
		double dmin = 0.5, dmax = 1;
		double x, y, c2 = dmax * dmin / (dmax - dmin);
		do {
			x = ran1(dmin, dmax);
			y = ran1(0.0, c2 / (dmin * dmin));
		} while (y > c2 / (x * x));
		return  x;
	}
/*
	static double ran3 (double z) {
		return (Diam_max - Diam_min)*z*z / (Box_z * Box_z) + Diam_min;
	}
*/
};

bool Random::flag;
double Random::r1;
double Random::r2;

#endif
