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
#include "error.h"
#include "logfile.h"

#include <mpi.h>
#include <cstdlib>
#include <iostream>

namespace hemo {

namespace {
int throwingGuards = 0;

bool mpiActive() {
  int initialized = 0, finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  return initialized && !finalized;
}

// Write "(HemoCell) (<kind>) [rank r] message" to stderr and, on the rank
// that writes the log file, to the log file.
void report(const char * kind, const std::string & message) {
  std::ostringstream line;
  line << "(HemoCell) (" << kind << ")";
  if (mpiActive()) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    line << " [rank " << rank << "]";
  }
  line << " " << message;
  std::cerr << line.str() << std::endl;
  if (hlog.logfile.is_open()) { hlog.logfile << line.str() << std::endl; }
}
} // namespace

void fatal(const std::string & message) {
  if (throwingGuards > 0) { throw FatalError(message); }
  report("Fatal error", message);
  std::cout << std::flush;
  if (mpiActive()) { MPI_Abort(MPI_COMM_WORLD, 1); }
  std::exit(1);
}

void warning(const std::string & message) {
  report("Warning", message);
}

FatalErrorsThrow::FatalErrorsThrow() {
  throwingGuards++;
}

FatalErrorsThrow::~FatalErrorsThrow() {
  throwingGuards--;
}

} // namespace hemo
