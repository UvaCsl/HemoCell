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
#ifndef LOGFILE_H
#define LOGFILE_H

#include <fstream>
#include <iostream>
namespace hemo {
class Logfile
{
public:
    Logfile()
    : logfile(*(new std::fstream()))
    {}
    Logfile(std::fstream& logfile_) 
    : logfile(logfile)
    {}
    std::fstream & logfile;
   
};
template<typename Val>
Logfile & operator << (Logfile & lf, Val const & rhs) {
    if (lf.logfile.good()) {
        std::cout << rhs;
        lf.logfile << rhs;
    }
    return lf;
}
inline Logfile & operator<< (Logfile & lf, std::ostream & (*op) (std::ostream&)) {
    if (lf.logfile.good()){
        std::cout << op << std::flush;
        lf.logfile << op;
    }
    return lf;
}
extern Logfile hlog;

class Logfile_only
{};
template<typename Val>
Logfile_only & operator << (Logfile_only & lf, Val const & rhs) {
    if (hlog.logfile.good()) {
        hlog.logfile << rhs;
    }
    return lf;
}
inline Logfile_only & operator<< (Logfile_only & lf, std::ostream & (*op) (std::ostream&)) {
    if (hlog.logfile.good()){
        hlog.logfile << op;
    }
    return lf;
}
extern Logfile_only hlogfile;

}
#endif /* LOGFILE_H */

