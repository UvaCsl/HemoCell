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
#ifndef FCN_GENERIC_FUNCTIONS_H
#define FCN_GENERIC_FUNCTIONS_H

#include <sys/stat.h>
#include <vector>
#include <string>

#include <errno.h>
#include <math.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */

namespace hemo {
  
typedef struct stat Stat;


using namespace std;


/*
 * Returns the new dimensions in case of weak scaling
 * eg.
 *       dim -> (32, 32, 32), numberOfProcesses-> 4
 *    newdim -> (64, 64, 32)
 *
 */
void weakScaling(int Nx, int Ny, int Nz, int numberOfProcesses, vector<int> & newNxNyNz);

inline bool file_exists (const std::string& name) {
    /* Checks if a file exists */
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int renameFileToDotOld(std::string fName);



/*
@(#)File:           $RCSfile: mkpath.c,v $
@(#)Version:        $Revision: 1.13 $
@(#)Last changed:   $Date: 2012/07/15 00:40:37 $
@(#)Purpose:        Create all directories in path
@(#)Author:         J Leffler
@(#)Copyright:      (C) JLSS 1990-91,1997-98,2001,2005,2008,2012
*/
// http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux

int do_mkdir(const char *path, mode_t mode);

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
// mkpath(argv[i], 0777);
int mkpath(const char *path, mode_t mode);

std::string zeroPadNumber(int num,int w = 12);

void printHeader();

}
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "constant_defaults.h"
namespace hemo {
void boundaryFromFlagMatrix(plb::MultiBlockLattice3D<T,DESCRIPTOR> * fluid, plb::MultiScalarField3D<int> * flagMatrix, bool); 
inline std::ostream& operator<<(std::ostream& stream, const plb::Box3D& box) {
    return stream << "Box3D: " << box.x0 << " "<<box.x1<<" "<<box.y0<<" "<<box.y1<< " "<<box.z0<<" "<<box.z1<<endl;
}
inline std::ostream& operator<<(std::ostream& stream, const plb::Dot3D& dot) {
    return stream << "Dot3D: " << dot.x << " "<<dot.y<<" "<<dot.z<<endl;
}
}
#endif // FCN_GENERIC_FUNCTIONS_HH
