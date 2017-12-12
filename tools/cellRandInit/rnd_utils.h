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
