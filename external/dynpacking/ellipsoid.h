#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

using namespace std;

class vector3 {
	double vec[3];
public:
	vector3() {
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = 0;
	}
	vector3(double* v) {
		vec[0] = v[0];
		vec[1] = v[1];
		vec[2] = v[2];
	}
	vector3(double v0, double v1, double v2) {
		vec[0] = v0;
		vec[1] = v1;
		vec[2] = v2;
	}
	vector3& random() {
		vec[0] = Random::ran();
		vec[1] = Random::ran();
		vec[2] = Random::ran();
		return *this;
	}
	vector3& scale(const vector3& v) {
		vec[0] *= v[0];
		vec[1] *= v[1];
		vec[2] *= v[2];
		return *this;
	}
	vector3& operator+= (const vector3& v) {
		vec[0] += v[0];
		vec[1] += v[1];
		vec[2] += v[2];
		return *this;
	}
	vector3& operator-= (const vector3& v) {
		vec[0] -= v[0];
		vec[1] -= v[1];
		vec[2] -= v[2];
		return *this;
	}
	vector3& operator*= (double d) {
		vec[0] *= d;
		vec[1] *= d;
		vec[2] *= d;
		return *this;
	}
	double& operator[] (int i) {
		return vec[i];
	}
	double operator[] (int i) const {
		return vec[i];
	}
	vector3 pbc_x(double x) {
		vec[0] += x;
		return *this;
	}
	vector3 pbc_y(double y) {
		vec[1] += y;
		return *this;
	}
	vector3 pbc_z(double z) {
		vec[2] += z;
		return *this;
	}
	vector3& pbc(const vector3& v) {
		if (vec[0] < 0) vec[0] += v[0];
		else if (vec[0] > v[0]) vec[0] -= v[0];
		if (vec[1] < 0) vec[1] += v[1];
		else if (vec[1] > v[1]) vec[1] -= v[1];
		if (vec[2] < 0) vec[2] += v[2];
		else if (vec[2] > v[2]) vec[2] -= v[2];
		return *this;
	}
	vector3& pbc_diff(const vector3& v1, const vector3& v2) {
		if (vec[0] < -v1[0]) vec[0] += v2[0];
		else if (vec[0] > v1[0]) vec[0] -= v2[0];
		if (vec[1] < -v1[1]) vec[1] += v2[1];
		else if (vec[1] > v1[1]) vec[1] -= v2[1];
		if (vec[2] < -v1[2]) vec[2] += v2[2];
		else if (vec[2] > v1[2]) vec[2] -= v2[2];
		return *this;
	}
	friend ostream& operator<< (ostream& out, const vector3& v) {
		out << setprecision(4) << "<" << 
			setw(8) << v[0] << ", " <<
			setw(8) << v[1] << ", " <<
			setw(8) << v[2] << ">";
		return out;
	}
};

inline vector3 operator+(const vector3& v1, const vector3& v2) {
	vector3 v(v1);
	v += v2;
	return v;
}
inline vector3 operator-(const vector3& v1, const vector3& v2) {
	vector3 v(v1);
	v -= v2;
	return v;
}
inline double operator*(const vector3& v1, const vector3& v2) {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}
inline vector3 operator*(const vector3& v1, double d) {
	vector3 v(v1);
	v *= d;
	return v;
}
inline vector3 operator*(double d, const vector3& v1) {
	vector3 v(v1);
	v *= d;
	return v;
}
inline vector3 prod(const vector3& v1, const vector3& v2) {
	return vector3(v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]);
}
inline vector3 norm(const vector3& vec) {
	double sp = vec*vec;
	if (fabs(sp - 1) < 1e-10) return vec;
	return vec * (1 / sqrt(sp));
}

class matrix33 {
	double mat[3][3];
public:
	matrix33() {
		memset((void*) mat, 0, 9*sizeof(double));
	}
	matrix33(double val[3][3]) {
		memcpy((void*) mat, (void*) val, 9*sizeof(double));
	}
	matrix33& random() {
		mat[0][0] = Random::ran();
		mat[0][1] = Random::ran();
		mat[0][2] = Random::ran();
		mat[1][0] = Random::ran();
		mat[1][1] = Random::ran();
		mat[1][2] = Random::ran();
		mat[2][0] = Random::ran();
		mat[2][1] = Random::ran();
		mat[2][2] = Random::ran();
		return *this;
	}
	matrix33& operator*=(const matrix33& m) {
		matrix33 t(*this);
		mat[0][0] = t(0,0)*m(0,0) + t(0,1)*m(1,0) + t(0,2)*m(2,0);
		mat[0][1] = t(0,0)*m(0,1) + t(0,1)*m(1,1) + t(0,2)*m(2,1);
		mat[0][2] = t(0,0)*m(0,2) + t(0,1)*m(1,2) + t(0,2)*m(2,2);
		mat[1][0] = t(1,0)*m(0,0) + t(1,1)*m(1,0) + t(1,2)*m(2,0);
		mat[1][1] = t(1,0)*m(0,1) + t(1,1)*m(1,1) + t(1,2)*m(2,1);
		mat[1][2] = t(1,0)*m(0,2) + t(1,1)*m(1,2) + t(1,2)*m(2,2);
		mat[2][0] = t(2,0)*m(0,0) + t(2,1)*m(1,0) + t(2,2)*m(2,0);
		mat[2][1] = t(2,0)*m(0,1) + t(2,1)*m(1,1) + t(2,2)*m(2,1);
		mat[2][2] = t(2,0)*m(0,2) + t(2,1)*m(1,2) + t(2,2)*m(2,2);
		return *this;
	}
	matrix33& operator*=(double d) {
		mat[0][0] *= d;
		mat[0][1] *= d;
		mat[0][2] *= d;
		mat[1][0] *= d;
		mat[1][1] *= d;
		mat[1][2] *= d;
		mat[2][0] *= d;
		mat[2][1] *= d;
		mat[2][2] *= d;
		return *this;
	}
	matrix33& operator+= (const matrix33& m) {
		mat[0][0] += m(0,0);
		mat[0][1] += m(0,1);
		mat[0][2] += m(0,2);
		mat[1][0] += m(1,0);
		mat[1][1] += m(1,1);
		mat[1][2] += m(1,2);
		mat[2][0] += m(2,0);
		mat[2][1] += m(2,1);
		mat[2][2] += m(2,2);
		return *this;
	}
	matrix33& operator-= (const matrix33& m) {
		mat[0][0] -= m(0,0);
		mat[0][1] -= m(0,1);
		mat[0][2] -= m(0,2);
		mat[1][0] -= m(1,0);
		mat[1][1] -= m(1,1);
		mat[1][2] -= m(1,2);
		mat[2][0] -= m(2,0);
		mat[2][1] -= m(2,1);
		mat[2][2] -= m(2,2);
		return *this;
	}
	double& operator()(int i, int j) {
		return mat[i][j];
	}
	double operator()(int i, int j) const {
		return mat[i][j];
	}
	bool not_sym() const {
		if (fabs(mat[0][1] - mat[1][0]) > 1e-5) return true;
		if (fabs(mat[0][2] - mat[2][0]) > 1e-5) return true;
		if (fabs(mat[1][2] - mat[2][1]) > 1e-5) return true;
		return false;
	}
	friend ostream& operator<< (ostream& out, const matrix33& m) {
		out << setprecision(4);
		out << setw(10) << m(0,0);
		out << setw(10) << m(0,1);
		out << setw(10) << m(0,2);
		out << endl;
		out << setw(10) << m(1,0);
		out << setw(10) << m(1,1);
		out << setw(10) << m(1,2);
		out << endl;
		out << setw(10) << m(2,0);
		out << setw(10) << m(2,1);
		out << setw(10) << m(2,2);
		return out << endl;
	}
};

inline matrix33 vvt(const vector3& v1, const vector3& v2) {
	matrix33 m;
	m(0,0) = v1[0]*v2[0];
	m(0,1) = v1[0]*v2[1];
	m(0,2) = v1[0]*v2[2];
	m(1,0) = v1[1]*v2[0];
	m(1,1) = v1[1]*v2[1];
	m(1,2) = v1[1]*v2[2];
	m(2,0) = v1[2]*v2[0];
	m(2,1) = v1[2]*v2[1];
	m(2,2) = v1[2]*v2[2];
	return m;
}
inline matrix33 operator+ (const matrix33& m1, const matrix33& m2) {
	matrix33 sum = m1;
	sum += m2;
	return sum;
}
inline matrix33 operator- (const matrix33& m1, const matrix33& m2) {
	matrix33 sum = m1;
	sum -= m2;
	return sum;
}
inline matrix33 operator* (const matrix33& m1, const matrix33& m2) {
	matrix33 m = m1;
	m *= m2;
	return m;
}
inline matrix33 operator* (const matrix33& m1, double d) {
	matrix33 m = m1;
	m *= d;
	return m;
}
inline matrix33 operator* (double d, const matrix33& m1) {
	matrix33 m = m1;
	m *= d;
	return m;
}
inline vector3 operator* (const matrix33& m, const vector3& v) {
	vector3 vv;
	vv[0] = m(0,0)*v[0] + m(0,1)*v[1] + m(0,2)*v[2];
	vv[1] = m(1,0)*v[0] + m(1,1)*v[1] + m(1,2)*v[2];
	vv[2] = m(2,0)*v[0] + m(2,1)*v[1] + m(2,2)*v[2];
	return vv;
}
inline vector3 operator* (const vector3& v, const matrix33& m) {
	vector3 vv;
	vv[0] = m(0,0)*v[0] + m(1,0)*v[1] + m(2,0)*v[2];
	vv[1] = m(0,1)*v[0] + m(1,1)*v[1] + m(2,1)*v[2];
	vv[2] = m(0,2)*v[0] + m(1,2)*v[1] + m(2,2)*v[2];
	return vv;
}
inline matrix33 eye() {
	matrix33 m;
	m(0,0) = m(1,1) = m(2,2) = 1;
	return m;
}
inline matrix33 adj(const matrix33& m) {
	matrix33 m1;
	m1(0,0) = m(1,1)*m(2,2) - m(1,2)*m(2,1);
	m1(0,1) = m(1,2)*m(2,0) - m(1,0)*m(2,2);
	m1(0,2) = m(1,0)*m(2,1) - m(1,1)*m(2,0);
	m1(1,0) = m(0,2)*m(2,1) - m(2,2)*m(0,1);
	m1(1,1) = m(0,0)*m(2,2) - m(0,2)*m(2,0);
	m1(1,2) = m(0,1)*m(2,0) - m(0,0)*m(2,1);
	m1(2,0) = m(0,1)*m(1,2) - m(0,2)*m(1,1);
	m1(2,1) = m(1,0)*m(0,2) - m(0,0)*m(1,2);
	m1(2,2) = m(1,1)*m(0,0) - m(1,0)*m(0,1);
	return m1;
}
inline double det(const matrix33& m) {
	return m(0,0)*m(1,1)*m(2,2) + m(1,0)*m(2,1)*m(0,2) + m(2,0)*m(0,1)*m(1,2) -
		   m(0,2)*m(1,1)*m(2,0) - m(0,0)*m(1,2)*m(2,1) - m(0,1)*m(1,0)*m(2,2);
}
inline matrix33 inverse(const matrix33& m) {
	return (1. / det(m)) * adj(m);
}
inline matrix33 skew(const vector3& v) {
	matrix33 m;
	m(0,1) = -v[2];
	m(0,2) =  v[1];
	m(1,0) =  v[2];
	m(1,2) = -v[0];
	m(2,0) = -v[1];
	m(2,1) =  v[0];
	return m;
}
inline double vmv(const matrix33& m, const vector3& v) {
	double w = m(0,0)*v[0]*v[0] + m(1,1)*v[1]*v[1] + m(2,2)*v[2]*v[2];
	w += m(0,1)*v[0]*v[1];
	w += m(0,2)*v[0]*v[2];
	w += m(1,0)*v[1]*v[0];
	w += m(1,2)*v[1]*v[2];
	w += m(2,0)*v[2]*v[0];
	w += m(2,1)*v[2]*v[1];
	return w;
}

inline matrix33 transp(const matrix33& m) {
	matrix33 mt;
	mt(0,0) = m(0,0);
	mt(0,1) = m(1,0);
	mt(0,2) = m(2,0);
	mt(1,0) = m(0,1);
	mt(1,1) = m(1,1);
	mt(1,2) = m(2,1);
	mt(2,0) = m(0,2);
	mt(2,1) = m(1,2);
	mt(2,2) = m(2,2);
	return mt;
}

class Quaternion;
Quaternion operator* (const Quaternion&, const Quaternion&);

class Quaternion {
	double s;
	vector3 p;
public:
	Quaternion() {
		double phi = Random::ran1(0, M_PI);
		s = cos(phi);
		p.random();
		p *= sin(phi);
	}
	Quaternion(double ss, vector3 pp) : s(ss), p(pp) { norm(); }
	Quaternion(vector3 euler_angles) {
		vector3 px(sin(0.5*euler_angles[0]), 0, 0);
		vector3 py(0, sin(0.5*euler_angles[1]), 0);
		vector3 pz(0, 0, sin( 0.5*euler_angles[2]));
		double sx = cos(0.5*euler_angles[0]);
		double sy = cos(0.5*euler_angles[1]);
		double sz = cos(0.5*euler_angles[2]);
		Quaternion qx(sx, px);
		Quaternion qy(sy, py);
		Quaternion qz(sz, pz);
		Quaternion q1 = qx * qy;
		*this = (q1 * qz).norm();
	}
	Quaternion& operator*= (const Quaternion& q) {
		Quaternion t(*this);
		s = t.s*q.s - t.p*q.p;
		p = t.s*q.p + q.s*t.p - prod(t.p, q.p);
		return *this;
	}
	Quaternion& norm() {
		double sp = s*s + p*p;
		if (fabs(sp - 1) > 1e-10) {
			sp = 1 / sqrt(sp);
			s *= sp;
			p *= sp;
		}
		return *this;
	}
	matrix33 countQ() const {
		matrix33 Q;
		Q = vvt(p, p) - s * skew(p) + (s*s - 0.5) * eye();
		return 2*Q;
	}
	friend class Ellipsoid_basic;
	friend ostream& operator<< (ostream& o, const Quaternion& q) {
		o << "s: " << q.s << " p: " << q.p << endl;
		return o;
	}
};
		
Quaternion operator* (const Quaternion& q1, const Quaternion& q2) {
	Quaternion q(q1);
	q *= q2;
	return q;
}

class Species {
	int number;
	vector3 rad;
	bool sph;
public:
	Species () : number(0), rad(vector3()), sph(0) {}
	Species (int n, vector3 r) : number(n), rad(r) {
		sph = ((fabs(rad[0] - rad[2]) < 1e-5) && (fabs(rad[0] - rad[1]) < 1e-5) ) ? 1 : 0;
	}
	int getn() { return number; }
	vector3 getr() { return rad; }
	bool gets() { return sph; }
	friend ostream& operator<< (ostream& o, const Species& s) {
		o << s.number << endl;
		o << s.rad << endl;
		o << s.sph << endl;
		return o;
	}
};

class Ellipsoid_basic {
protected:
	Species* kind;
	vector3 pos;
	Quaternion q;
public:
	Ellipsoid_basic() {}
	Ellipsoid_basic(Species* k, vector3 box) : kind(k), q(vector3(0.0,0.0,0.0)) { pos.random().scale(box); }
	Ellipsoid_basic(Species* k, const vector3& ps, const Quaternion& qq) : kind(k), pos(ps), q(qq) {}
	matrix33 countX() const {
		vector3 r(kind->getr());
		matrix33 Q, O;
		O(0,0) = 1 / (r[0]*r[0]);
		O(1,1) = 1 / (r[1]*r[1]);
		O(2,2) = 1 / (r[2]*r[2]);
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	matrix33 countX_12() const {
		vector3 r(kind->getr());
		matrix33 Q, O;
		O(0,0) = 1 / r[0];
		O(1,1) = 1 / r[1];
		O(2,2) = 1 / r[2];
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	matrix33 countX_m1() const {
		vector3 r(kind->getr());
		matrix33 Q, O;
		O(0,0) = r[0]*r[0];
		O(1,1) = r[1]*r[1];
		O(2,2) = r[2]*r[2];
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	Ellipsoid_basic& rotate(const Quaternion& qr) {
		q *= qr;
		return *this;
	}
	vector3& get_pos() {
		return pos;
	}
	const vector3& get_pos() const {
		return pos;
	}
	Species* get_k() const { return kind; }
	void set_k(Species * k) { kind = k; }
	const Quaternion& get_q() const { return q; }
	Quaternion& get_q() { return q; }
	friend ostream &operator<< (ostream &o, Ellipsoid_basic &e) {
		o << *e.kind << endl;
		o << e.pos << endl;
		o << e.q << endl;
		return o;
	}
};

class Poly2 {
	double p[3];
public:
	Poly2(double p0, double p1, double p2) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
	}
	Poly2(double* pp = 0) {
		if (pp == 0) return;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
	}
	void set(double p0, double p1, double p2) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
	}
	double operator[] (int i) const { return p[i]; }
	double val(double x) {
		return p[0] + x*(p[1] + x*p[2]);
	}
	double vald(double x) {
		return p[1] + 2*x*p[2];
	}
	friend ostream& operator<<(ostream& out, const Poly2& p) {
		out << "p[0] = " << p.p[0] << endl;
		out << "p[1] = " << p.p[1] << endl;
		out << "p[2] = " << p.p[2] << endl;
		return out << endl;
	}
};

class Poly3 {
	double p[4];
public:
	Poly3(double p0, double p1, double p2, double p3) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
	}
	Poly3(double* pp = 0) {
		if (pp == 0) return;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
		p[3] = pp[3];
	}
	void set(double p0, double p1, double p2, double p3) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
	}
	double operator[] (int i) const { return p[i]; }
	double val(double x) {
		return p[0] + x*(p[1] + x*(p[2] + x*p[3]));
	}
	double vald(double x) {
		return p[1] + x*(2*p[2] + x*3*p[3]);
	}
	friend ostream& operator<<(ostream& out, const Poly3& p) {
		out << "p[0] = " << p.p[0] << endl;
		out << "p[1] = " << p.p[1] << endl;
		out << "p[2] = " << p.p[2] << endl;
		out << "p[3] = " << p.p[3] << endl;
		return out << endl;
	}
};

class Poly4 {
	double p[5];
public:
	Poly4(double p0, double p1, double p2, double p3, double p4) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
	}
	Poly4(double* pp = 0) {
		if (pp == 0) return;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
		p[3] = pp[3];
		p[4] = pp[4];
	}
	void set(double p0, double p1, double p2, double p3, double p4) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
	}
	double operator[] (int i) const { return p[i]; }
	double val(double x) {
		return p[0] + x*(p[1] + x*(p[2] + x*(p[3] + x*p[4])));
	}
	double vald(double x) {
		return p[1] + x*(2*p[2] + x*(3*p[3] + x*4*p[4]));
	}
	friend ostream& operator<<(ostream& out, const Poly4& p) {
		out << "p[0] = " << p.p[0] << endl;
		out << "p[1] = " << p.p[1] << endl;
		out << "p[2] = " << p.p[2] << endl;
		out << "p[3] = " << p.p[3] << endl;
		out << "p[4] = " << p.p[4] << endl;
		return out << endl;
	}
};

class Poly5 {
	double p[6];
public:
	Poly5(double p0, double p1, double p2, double p3, double p4, double p5) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
		p[5] = p5;
	}
	Poly5(double* pp = 0) {
		if (pp == 0) return;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
		p[3] = pp[3];
		p[4] = pp[4];
		p[5] = pp[5];
	}
	void set(double p0, double p1, double p2, double p3, double p4, double p5) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
		p[5] = p5;
	}
	double operator[] (int i) const { return p[i]; }
	double val(double x) {
		return p[0] + x*(p[1] + x*(p[2] + x*(p[3] + x*(p[4] + x*p[5]))));
	}
	double vald(double x) {
		return p[1] + x*(2*p[2] + x*(3*p[3] + x*(4*p[4] + x*5*p[5])));
	}
	friend ostream& operator<<(ostream& out, const Poly5& p) {
		out << "p[0] = " << p.p[0] << endl;
		out << "p[1] = " << p.p[1] << endl;
		out << "p[2] = " << p.p[2] << endl;
		out << "p[3] = " << p.p[3] << endl;
		out << "p[4] = " << p.p[4] << endl;
		out << "p[5] = " << p.p[5] << endl;
		return out << endl;
	}
};

class Poly6 {
	double p[7];
public:
	Poly6(double p0, double p1, double p2, double p3, double p4, double p5, double p6) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
		p[5] = p5;
		p[6] = p6;
	}
	Poly6(double* pp = 0) {
		if (pp == 0) return;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
		p[3] = pp[3];
		p[4] = pp[4];
		p[5] = pp[5];
		p[6] = pp[6];
	}
	void set(double p0, double p1, double p2, double p3, double p4, double p5, double p6) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		p[3] = p3;
		p[4] = p4;
		p[5] = p5;
		p[6] = p6;
	}
	double operator[] (int i) const { return p[i]; }
	double val(double x) {
		return p[0] + x*(p[1] + x*(p[2] + x*(p[3] + x*(p[4] + x*(p[5] + x*p[6])))));
	}
	double vald(double x) {
		return p[1] + x*(2*p[2] + x*(3*p[3] + x*(4*p[4] + x*(5*p[5] + x*6*p[6]))));
	}
	friend ostream& operator<<(ostream& out, const Poly6& p) {
		out << "p[0] = " << p.p[0] << endl;
		out << "p[1] = " << p.p[1] << endl;
		out << "p[2] = " << p.p[2] << endl;
		out << "p[3] = " << p.p[3] << endl;
		out << "p[4] = " << p.p[4] << endl;
		out << "p[5] = " << p.p[5] << endl;
		out << "p[6] = " << p.p[6] << endl;
		return out << endl;
	}
};

class Ellipsoid_2 {
	matrix33 XA_m1, XB_m1, XB_12;
	vector3 r_AB, n, r_AC, r_BC;
	matrix33 A_AB;
	vector3 a_AB;
	double lambda, lam0;
	Poly4 p;
	Poly3 q;
	Poly6 h;
	matrix33 countY(double lambda) const {
		return lambda * XB_m1 + (1-lambda) * XA_m1;
	}
public:
	Ellipsoid_2(const Ellipsoid_basic& a, const Ellipsoid_basic& b, const vector3 rij) :
		XA_m1(a.countX_m1()), XB_m1(b.countX_m1()),
		XB_12(b.countX_12()),
		r_AB(rij), n(0, 0, 0),
		r_AC(0, 0, 0), r_BC(0, 0, 0),
		A_AB(XB_12 * XA_m1 * XB_12), a_AB(XB_12 * r_AB),
		lambda(0), lam0((a.get_k()->getr())[0] / ((a.get_k()->getr())[0] + (b.get_k()->getr())[0])) {
/*
		if (A_AB.not_sym()) {
			cout << "A_AB" << endl << A_AB << endl;
			exit(1);
		}
*/
		double detA = det(A_AB);
		matrix33 Z = adj(A_AB);
		matrix33 C0(Z), C1(Z), C2(Z);
		C1 *= -2;
		C1(0,0) += A_AB(1,1) + A_AB(2,2);
		C1(1,1) += A_AB(0,0) + A_AB(2,2);
		C1(2,2) += A_AB(1,1) + A_AB(0,0);
		C1(0,1) -= A_AB(1,0);
		C1(0,2) -= A_AB(2,0);
		C1(1,0) -= A_AB(0,1);
		C1(1,2) -= A_AB(2,1);
		C1(2,0) -= A_AB(0,2);
		C1(2,1) -= A_AB(1,2);
		C2(0,0) -= A_AB(1,1) + A_AB(2,2) - 1;
		C2(1,1) -= A_AB(0,0) + A_AB(2,2) - 1;
		C2(2,2) -= A_AB(1,1) + A_AB(0,0) - 1;
		C2(0,1) += A_AB(1,0);
		C2(0,2) += A_AB(2,0);
		C2(1,0) += A_AB(0,1);
		C2(1,2) += A_AB(2,1);
		C2(2,0) += A_AB(0,2);
		C2(2,1) += A_AB(1,2);
		double w0 = vmv(C0, a_AB);
		double w1 = vmv(C1, a_AB);
		double w2 = vmv(C2, a_AB);
		double sa = A_AB(0,0) + A_AB(1,1) + A_AB(2,2);
		double sz = Z(0,0) + Z(1,1) + Z(2,2);
		if (fabs(detA) > 1e-5) {
			p.set(0, w0, w1 - w0, w2 - w1, -w2);
			q.set(detA, -3*detA + sz, 3*detA + sa - 2*sz, 1 - detA - sa + sz);
			h.set(q[0]*p[1], 2*q[0]*p[2], 3*q[0]*p[3] + q[1]*p[2] - q[2]*p[1],
				4*q[0]*p[4] + 2*q[1]*p[3] - 2*q[3]*p[1],
				3*q[1]*p[4] + q[2]*p[3] - q[3]*p[2],
				2*q[2]*p[4], q[3]*p[4]);
		} else {
			p.set(w0, w1 - w0, w2 - w1, -w2, 0);
			q.set(sz, sa - 2*sz, 1 - sa + sz, 0);
			h.set(q[0]*p[1], 2*q[0]*p[2], 3*q[0]*p[3] + q[1]*p[2] - q[2]*p[1],
				2*q[1]*p[3], q[2]*p[3], 0, 0);
		}
	}
	double f_AB(double l) { 
//		if (fabs(q.val(l)) < 1e-5)
//			cout << "******************* " << l << " " << p.val(l) << " " << q.val(l) << endl;
		return p.val(l) / q.val(l);
	}
	double f_AB_d(double l) { return h.val(l); }
	double f_AB_d2(double l) { return h.vald(l); }
	double newton() {
		double x = lam0, fv = f_AB_d(x), fvd = f_AB_d2(x);
		while (fabs(fv) > 1e-10) {
			x -= (fv = f_AB_d(x)) / (fvd = f_AB_d2(x));
/*
			cout << x << " " << fv << " " << fvd << endl;
			if (x < 0 || x > 1) {
				cout << "------------------------------------------" << endl;
				cout << "A_AB " << endl << A_AB << endl;
				for (double ll = 0; ll <= 1; ll += 0.01) 
					cout << ll << " " << f_AB(ll) << " " << f_AB_d(ll) << " " <<
					f_AB_d2(ll) << " " << p.val(ll) << " " << q.val(ll) << " " << det(A_AB) << endl;
				cout << "------------------------------------------" << endl;
			}
*/
		}
		return lambda = x;
	}
	vector3 get_r() { return r_AB; }
	vector3 count_n() {
		return n = inverse(countY(lambda)) * r_AB;
	}
	vector3 count_rac() {
		return r_AC = (1-lambda) * XA_m1 * n;
	}
	vector3 count_rbc() {
		return r_BC = -lambda * XB_m1 * n;
	}
	vector3 get_n() { return n; }
	vector3& get_rab() { return r_AB; }
	double& get_lambda() { return lambda; }
	friend ostream& operator <<(ostream& o, const Ellipsoid_2& e) {
		o << "rab:" << endl;
		o << e.r_AB << endl;
		o << "XA_m1:" << endl;
		o << e.XA_m1 << endl;
		o << "XB_m1:" << endl;
		o << e.XB_m1 << endl;
		o << "lambda = " << e.lambda << endl;
		o << "Y:" << endl;
		o << e.countY(e.lambda) << endl;
		o << "Ym1:" << endl;
		o << inverse(e.countY(e.lambda)) << endl;
		return o;
	}
};

#endif
