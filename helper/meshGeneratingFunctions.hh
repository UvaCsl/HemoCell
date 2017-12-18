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
#ifndef MESH_GENERATING_FUNCTIONS_HH
#define MESH_GENERATING_FUNCTIONS_HH

#include "meshGeneratingFunctions.h"

namespace plb {

template<typename T>
TriangleSet<T> constructSphereIcosahedron(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles)
{
    std::vector<typename TriangleSet<T>::Triangle> triangles;
#ifdef PLB_DEBUG
    static const T eps = std::numeric_limits<T>::epsilon();
#endif
    PLB_ASSERT(radius > (T) 0.0 && !util::fpequal(radius, (T) 0.0, eps) &&
               minNumOfTriangles >= 8);

    // Create a triangularized unit sphere
    // Twelve vertices of icosahedron on unit sphere
    T tau = -0.8506508084; // t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
    T one = -0.5257311121; // one=1/sqrt(1+t^2) , unit sphere

    // Initial 6 vertices
    plb::Array<T,3> v1( tau, one, 0); //ZA
    plb::Array<T,3> v2(-tau, one, 0); //ZB
    plb::Array<T,3> v3(-tau, -one, 0); //ZC
    plb::Array<T,3> v4( tau, -one, 0); //ZD
    plb::Array<T,3> v5( one, 0 , tau); //YA
    plb::Array<T,3> v6( one, 0 , -tau); //YB
    plb::Array<T,3> v7(-one, 0 , -tau); //YC
    plb::Array<T,3> v8(-one, 0 , tau); //YD
    plb::Array<T,3> v9( 0 , tau, one); //XA
    plb::Array<T,3> v10( 0 , -tau, one); //XB
    plb::Array<T,3> v11( 0 , -tau, -one); //XC
    plb::Array<T,3> v12( 0 , tau, -one); //XD

    // Structure for unit icosahedron
    typename TriangleSet<T>::Triangle tmp;

    tmp[0] = v5; tmp[1] = v8; tmp[2] = v9;
    triangles.push_back(tmp);
    tmp[0] = v5; tmp[1] = v10; tmp[2] = v8;
    triangles.push_back(tmp);
    tmp[0] = v6; tmp[1] = v12; tmp[2] = v7;
    triangles.push_back(tmp);
    tmp[0] = v6; tmp[1] = v7; tmp[2] = v11;
    triangles.push_back(tmp);
    tmp[0] = v1; tmp[1] = v4; tmp[2] = v5;
    triangles.push_back(tmp);
    tmp[0] = v1; tmp[1] = v6; tmp[2] = v4;
    triangles.push_back(tmp);
    tmp[0] = v3; tmp[1] = v2; tmp[2] = v8;
    triangles.push_back(tmp);
    tmp[0] = v3; tmp[1] = v7; tmp[2] = v2;
    triangles.push_back(tmp);
    tmp[0] = v9; tmp[1] = v12; tmp[2] = v1;
    triangles.push_back(tmp);
    tmp[0] = v9; tmp[1] = v2; tmp[2] = v12;
    triangles.push_back(tmp);
    tmp[0] = v10; tmp[1] = v4; tmp[2] = v11;
    triangles.push_back(tmp);
    tmp[0] = v10; tmp[1] = v11; tmp[2] = v3;
    triangles.push_back(tmp);
    tmp[0] = v9; tmp[1] = v1; tmp[2] = v5;
    triangles.push_back(tmp);
    tmp[0] = v12; tmp[1] = v6; tmp[2] = v1;
    triangles.push_back(tmp);
    tmp[0] = v5; tmp[1] = v4; tmp[2] = v10;
    triangles.push_back(tmp);
    tmp[0] = v6; tmp[1] = v11; tmp[2] = v4;
    triangles.push_back(tmp);
    tmp[0] = v8; tmp[1] = v2; tmp[2] = v9;
    triangles.push_back(tmp);
    tmp[0] = v7; tmp[1] = v12; tmp[2] = v2;
    triangles.push_back(tmp);
    tmp[0] = v8; tmp[1] = v10; tmp[2] = v3;
    triangles.push_back(tmp);
    tmp[0] = v7; tmp[1] = v3; tmp[2] = v11;
    triangles.push_back(tmp);

    // Perform refinement iterations

    plint size;
    plb::Array<T,3>  va,vb,vc,vd,ve,vf;
    while ((size = triangles.size()) < minNumOfTriangles) {
        for (plint i = 0; i < size; i++) {
            va = triangles[i][0];
            vb = triangles[i][1];
            vc = triangles[i][2];

            vd = (T) 0.5 * (va + vb);
            ve = (T) 0.5 * (vb + vc);
            vf = (T) 0.5 * (vc + va);

            vd /= norm(vd);
            ve /= norm(ve);
            vf /= norm(vf);

            triangles[i][0] = vd;
            triangles[i][1] = ve;
            triangles[i][2] = vf;

            tmp[0] = va;
            tmp[1] = vd;
            tmp[2] = vf;
            triangles.push_back(tmp);

            tmp[0] = vd;
            tmp[1] = vb;
            tmp[2] = ve;
            triangles.push_back(tmp);

            tmp[0] = vf;
            tmp[1] = ve;
            tmp[2] = vc;
            triangles.push_back(tmp);
        }
    }

    // Scale and translate the mesh

    TriangleSet<T> triangleSet(triangles);

    triangleSet.scale(radius);
    triangleSet.translate(center);

    return triangleSet;
}

template<typename T>
plb::Array<T,3> spherePointToRBCPoint(const plb::Array<T,3> point, T R) {
    plb::Array<T,3> rbcPoint(point);
    T r2 = rbcPoint[0]*rbcPoint[0] + rbcPoint[1]*rbcPoint[1];
    T val = rbcPoint[2];
    plint sign = (T(0) < val) - (val < T(0));
    rbcPoint[0] *= R;
    rbcPoint[1] *= R;
    if (1-r2 <0) {
        r2 =1 ;
    }
    //T C0 = 0.0135805, C2 = 1.001279, C4 = -0.561381;
    T C0 = 0.054322, C2 = 1.001279, C4 = -0.561381;
    rbcPoint[2] = sign * R * sqrt(1-r2) * (C0 + C2*r2 + C4*r2*r2);
    return rbcPoint;
}

template<typename T>
plb::Array<T,3> spherePointToEllipsoidPoint(const plb::Array<T,3> point, T R, T aspectRatio) {
    plb::Array<T,3> elPoint(point);
    T r2 = elPoint[0]*elPoint[0] + elPoint[1]*elPoint[1];
    T val = elPoint[2];
    plint sign = (T(0) < val) - (val < T(0));
    if (1-r2 <0) {
        r2 =1 ;
    }
    elPoint[0] *= R;
    elPoint[1] *= R;
    elPoint[2] = sign * aspectRatio * R * sqrt(1-r2) ;
    return elPoint;
}

template<typename T>
plb::Array<T,3> mapMeshAsRBC(const plb::Array<T,3> point, const plb::Array<T,3> center, T R) {
    plb::Array<T,3> rbcPoint(point - center);
    rbcPoint[0] = rbcPoint[0] > R ? R : rbcPoint[0];
    rbcPoint[1] = rbcPoint[1] > R ? R : rbcPoint[1];
    T r2 = rbcPoint[0]*rbcPoint[0] + rbcPoint[1]*rbcPoint[1];
    //T C0 = 0.02716, C2 = 2.0024, C4 = -1.123;
    T C0 = 0.02716, C2 = 2.0024, C4 = -1.123;
    T val = rbcPoint[2];
    plint sign = (T(0) < val) - (val < T(0));
//    rbcPoint[0] *= R;
//    rbcPoint[1] *= R;
    r2 = r2 / (R*R);
    if (1-r2 <0) {
        r2 =1 ;
    }
    rbcPoint[2] = sign * 0.5 * R * sqrt(1-r2) * (C0 + C2*r2 + C4*r2*r2);
    rbcPoint = (rbcPoint + center);
    return rbcPoint;
}


template<typename T>
TriangleSet<T> constructRBC(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles, plb::Array<T,3> const& eulerAngles) {
    return constructCell(center, radius, "./lib/RBC.stl", eulerAngles);
}


template<typename T>
TriangleSet<T> constructRBCFromSphere(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles,
        plb::Array<T,3> const& eulerAngles, pluint initialSphereShape)
{
    TriangleSet<T> sphere;
    if (initialSphereShape == 1) {
        sphere = constructSphereIcosahedron<T>(plb::Array<T,3>(0,0,0), 1.0, minNumOfTriangles);
    } else if (initialSphereShape == 0) {
        sphere = constructSphere<T>(plb::Array<T,3>(0,0,0), 1.0, minNumOfTriangles);
    }
    sphere.rotate(
            PI/2.0 + eulerAngles[0],
            PI/2.0 + eulerAngles[1],
            0. + eulerAngles[2]);
    std::vector<typename TriangleSet<T>::Triangle> rbcTriangles = sphere.getTriangles();
    for (pluint var = 0; var < rbcTriangles.size(); ++var) {
        rbcTriangles[var][0] = spherePointToRBCPoint(rbcTriangles[var][0]);
        rbcTriangles[var][1] = spherePointToRBCPoint(rbcTriangles[var][1]);
        rbcTriangles[var][2] = spherePointToRBCPoint(rbcTriangles[var][2]);
    }
    TriangleSet<T> rbc(rbcTriangles);
    rbc.scale(radius);
    rbc.rotate(
            PI/2.0 + eulerAngles[0],
            PI/2.0 + eulerAngles[1],
            0. + eulerAngles[2]);
    rbc.translate(center);
    return rbc;
}


template<typename T>
TriangleSet<T> constructEllipsoidFromSphere(plb::Array<T,3> const& center, T radius, T aspectRatio, plint minNumOfTriangles,
        plb::Array<T,3> const& eulerAngles, pluint initialSphereShape)
{
    TriangleSet<T> sphere;
    if (initialSphereShape == 1) {
        sphere = constructSphereIcosahedron<T>(plb::Array<T,3>(0,0,0), 1.0, minNumOfTriangles);
    } else if (initialSphereShape == 0) {
        sphere = constructSphere<T>(plb::Array<T,3>(0,0,0), 1.0, minNumOfTriangles);
    }
    sphere.rotate(
            PI/2.0 + eulerAngles[0],
            PI/2.0 + eulerAngles[1],
            0. + eulerAngles[2]);
    std::vector<typename TriangleSet<T>::Triangle> ellipsoidTriangles = sphere.getTriangles();
    for (pluint var = 0; var < ellipsoidTriangles.size(); ++var) {
        ellipsoidTriangles[var][0] = spherePointToEllipsoidPoint(ellipsoidTriangles[var][0], radius, aspectRatio);
        ellipsoidTriangles[var][1] = spherePointToEllipsoidPoint(ellipsoidTriangles[var][1], radius, aspectRatio);
        ellipsoidTriangles[var][2] = spherePointToEllipsoidPoint(ellipsoidTriangles[var][2], radius, aspectRatio);
    }
    TriangleSet<T> ellipsoid(ellipsoidTriangles);
    ellipsoid.rotate(
            PI/2.0 + eulerAngles[0],
            PI/2.0 + eulerAngles[1],
            0. + eulerAngles[2]);
    ellipsoid.translate(center);
    return ellipsoid;
}


template<typename T>
TriangleSet<T> constructCell(plb::Array<T,3> const& center, T radius, std::string cellFilename, plb::Array<T,3> const& eulerAngles) {
//    Cuboid<T> boundingCuboid;
    TriangleSet<T> Cell(cellFilename);
    Cuboid<T> cb = Cell.getBoundingCuboid();
    plb::Array<T,3> dr = (cb.upperRightCorner - cb.lowerLeftCorner);
    T scaleFactor = std::max(dr[0],std::max(dr[1],dr[2]));
    Cell.scale(radius*2.0/scaleFactor);
    Cell.rotate(
            PI/2.0 + eulerAngles[0],
            PI/2.0 + eulerAngles[1],
            0. + eulerAngles[2]);
    Cell.translate(center);
    return Cell;
}

} // namespace plb

#endif  // MESH_GENERATING_FUNCTIONS_HH
