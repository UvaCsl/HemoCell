#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include "geometry.h"

using namespace std;

struct CellType {
  string name;
  double dx,dy,dz;
  int number;
};

class Species {
	int number;
	vector3 rad;
	bool sph;
public:
        string name;

	Species (string & name_, int n, vector3 r) :  number(n), rad(r),name(name_) {
		sph = ((fabs(rad[0] - rad[2]) < 1e-5) && (fabs(rad[0] - rad[1]) < 1e-5) ) ? 1 : 0;
	}
        Species (int n, vector3 r) :  number(n), rad(r), name("") {
		sph = ((fabs(rad[0] - rad[2]) < 1e-5) && (fabs(rad[0] - rad[1]) < 1e-5) ) ? 1 : 0;
	}
	int getNum() { return number; }
	vector3 getRot() { return rad; }
	bool getIsSphere() { return sph; }
	friend ostream& operator<< (ostream& o, const Species& s) {
		o << s.number << endl;
		o << s.rad << endl;
		o << s.sph << endl;
		return o;
	}
};

class Ellipsoid {
	vector3 f, fu;
protected:
	Species* kind;
	vector3 pos;
	Quaternion q;
public:
	Ellipsoid() {}
	Ellipsoid(Species* k, vector3 box) : kind(k), q(vector3(0.0,0.0,0.0)) { pos.random().scale(box); }
	Ellipsoid(Species* k, const vector3& ps, const Quaternion& qq) : kind(k), pos(ps), q(qq) {}
	matrix33 countX() const {
		vector3 r(kind->getRot());
		matrix33 Q, O;
		O(0,0) = 1 / (r[0]*r[0]);
		O(1,1) = 1 / (r[1]*r[1]);
		O(2,2) = 1 / (r[2]*r[2]);
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	matrix33 countX_12() const {
		vector3 r(kind->getRot());
		matrix33 Q, O;
		O(0,0) = 1 / r[0];
		O(1,1) = 1 / r[1];
		O(2,2) = 1 / r[2];
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	matrix33 countX_m1() const {
		vector3 r(kind->getRot());
		matrix33 Q, O;
		O(0,0) = r[0]*r[0];
		O(1,1) = r[1]*r[1];
		O(2,2) = r[2]*r[2];
		Q = q.countQ();
		return transp(Q) * O * Q;
	}
	Ellipsoid& rotate(const Quaternion& qr) {
		q *= qr;
		return *this;
	}
	vector3& get_pos() {
		return pos;
	}
	const vector3& get_pos() const {
		return pos;
	}
	void set_force(vector3 ff = vector3()) { f = ff; }
    void set_forceu(vector3 ffu = vector3()) { fu = ffu; }
    vector3& get_f() { return f; }
    vector3& get_fu() { return fu; }
	Species* getSpecies() const { return kind; }
	void setSpecies(Species * k) { kind = k; }
	const Quaternion& get_q() const { return q; }
	Quaternion& get_q() { return q; }
	
	friend ostream &operator<< (ostream &o, Ellipsoid &e) {
		o << *e.kind << endl;
		o << e.pos << endl;
		o << e.q << endl;
		return o;
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

class BinaryEllipsoidSystem {
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
	BinaryEllipsoidSystem(const Ellipsoid& a, const Ellipsoid& b, const vector3 rij) :
		XA_m1(a.countX_m1()), XB_m1(b.countX_m1()),
		XB_12(b.countX_12()),
		r_AB(rij), n(0, 0, 0),
		r_AC(0, 0, 0), r_BC(0, 0, 0),
		A_AB(XB_12 * XA_m1 * XB_12), a_AB(XB_12 * r_AB),
		lambda(0), lam0((a.getSpecies()->getRot())[0] / ((a.getSpecies()->getRot())[0] + (b.getSpecies()->getRot())[0])) {

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
		return p.val(l) / q.val(l);
	}
	double f_AB_d(double l) { return h.val(l); }
	double f_AB_d2(double l) { return h.vald(l); }
	double newton() {
		double x = lam0, fv = f_AB_d(x), fvd = f_AB_d2(x);
		while (fabs(fv) > 1e-10) {
			x -= (fv = f_AB_d(x)) / (fvd = f_AB_d2(x));

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
	friend ostream& operator <<(ostream& o, const BinaryEllipsoidSystem& e) {
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
