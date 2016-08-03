#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

#include "rnd_utils.h"


// Vector class --------------------

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

// Matrix 3x3 class -------------------------------------

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
        double phi = Random::ran1(0, PI);
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

#endif //GEOMETRY_H
