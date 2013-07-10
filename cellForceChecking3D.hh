#ifndef CELL_FORCE_CHECKING_3D_HH
#define CELL_FORCE_CHECKING_3D_HH

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellForceChecking3D.h"

void writeDataVTK(
        Array<T,3> const& x1, Array<T,3> const& x2,
        Array<T,3> const& x3, Array<T,3> const& x4,
        Array<T,3> const& fx1, Array<T,3> const& fx2,
        Array<T,3> const& fx3, Array<T,3> const& fx4,
        std::string fName) {
    std::ofstream ofile(fName.c_str());
    ofile << "# vtk DataFile Version 3.0\n";
    ofile << "Surface point created with Palabos and ficsion\n";
    ofile << "ASCII\n";
    ofile << "DATASET UNSTRUCTURED_GRID\n";

    ofile << "POINTS 4 double\n";
    ofile << x1[0] << " " <<x1[1] << " " <<x1[2] << "\n";
    ofile << x2[0] << " " <<x2[1] << " " <<x2[2] << "\n";
    ofile << x3[0] << " " <<x3[1] << " " <<x3[2] << "\n";
    ofile << x4[0] << " " <<x4[1] << " " <<x4[2] << "\n";
    ofile << "CELLS 2 8\n";
    ofile << "3 0 1 2\n";
    ofile << "3 2 3 0\n";
    ofile << "CELL_TYPES 4\n";
    ofile << "5\n";
    ofile << "5\n";
    ofile << "5\n";
    ofile << "5\n";


    ofile << "POINT_DATA 4\nVECTORS force double\n";
//    ofile << "LOOKUP_TABLE default \n";
    ofile << fx1[0] << " " <<fx1[1] << " " <<fx1[2] << "\n";
    ofile << fx2[0] << " " <<fx2[1] << " " <<fx2[2] << "\n";
    ofile << fx3[0] << " " <<fx3[1] << " " <<fx3[2] << "\n";
    ofile << fx4[0] << " " <<fx4[1] << " " <<fx4[2] << "\n";

    ofile.close();
}



void testBending() {
    Array<T,3> x1(1.,0.,0.);
    Array<T,3> x3(-1.,0.,0.);
    Array<T,3> x4(0.,-1.,0.);
    Array<T,3> x2, f1,f2,f3,f4,f;
    int eqAngle=-90;
//    for (eqAngle=0; eqAngle<360; eqAngle+=1) {
        for (int th=0; th<360; th+=1) {
            x2 = Array<T,3>(0., cos(th*pi/180), sin(th*pi/180));

        	Array<T,3> ni(0.,0.,0.),  nj(0.,0.,0.);
        	T Ai, Aj;
            crossProduct(x2-x1,  x3-x1, ni);
            crossProduct(x4-x3,  x1-x3, nj);
            Ai = norm(ni)/2; ni = ni*1.0/(2*Ai);
            Aj = norm(nj)/2; nj = nj*1.0/(2*Aj);
//            f1 = computeBendingForceFromPotential (x1, x2, x3, x4, 0., 0.,eqAngle*pi/180, 1.0, f2, f3, f4);
            f1 = computeBendingForce(x1, x2, x3, x4, ni, nj, Ai, Aj, 0.0, 0.0, eqAngle*pi/180, 1.0, f2, f3, f4);

            f = f1+f2+f3+f4;
            if (norm(f) < 1.e-10) {
                pcout << "BENDOK;"  ;
            } else {
                pcout << "BENDNOTOK;" ;
            }
            writeDataVTK(x1,x2,x3,x4,f1,f2,f3,f4, createFileName("./BTEST/DANGLE_",th,10)+".vtk");
            pcout << eqAngle * pi/180
                  << ";" << th* pi/180
                  << ";" << (th-eqAngle)* pi/180
                  << ";" << norm(f) << std::endl;
//            pcout << "eqAngle:" << eqAngle * 180/pi
//                  << "\tAngle:" << th* 180/pi
//                  << "\tdiff:" << (th-eqAngle)* 180/pi
//                  << "\tforce " << norm(f) << std::endl;
//        }
    }
}



#endif  // CELL_FORCE_CHECKING_3D_HH
