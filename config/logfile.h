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
#include <sstream>
#include <string>
namespace hemo {
/**
 * Log stream of HemoCell (`hlog`): writes to stdout and to the log file, on
 * the main MPI process only.
 *
 * The log file is opened by loadDirectories(), which for checkpointed runs
 * happens only in loadCheckPoint(). Messages logged before that (setup
 * parameters, warnings, errors before an exit) are printed immediately and
 * kept in `pending`, which open() writes to the top of the file.
 */
class Logfile {
public:
  Logfile() : logfile(*(new std::fstream())) {}
  Logfile(std::fstream & logfile_) : logfile(logfile_) {}
  std::fstream & logfile;
  std::string filename = "";
  /// True on the process that writes the log (MPI rank 0); set by HemoCell.
  bool mainProcess = false;
  /// Messages logged on the main process before the log file was opened.
  std::ostringstream pending;

  /// Open the log file (if not open yet) and write the pending messages.
  /// @return whether the log file is open afterwards.
  bool open(const std::string & filename_) {
    if (!logfile.is_open()) {
      filename = filename_;
      logfile.open(filename_, std::fstream::out);
      if (logfile.is_open()) {
        logfile << pending.str();
        pending.str("");
      }
    }
    return logfile.is_open();
  }
  /// Close the log file and stop logging (called when HemoCell is destroyed).
  void close() {
    if (logfile.is_open()) { logfile.close(); }
    pending.str("");
    mainProcess = false;
  }
};
template <typename Val>
Logfile & operator<<(Logfile & lf, Val const & rhs) {
  if (lf.logfile.is_open()) {
    std::cout << rhs;
    lf.logfile << rhs;
  } else if (lf.mainProcess) {
    std::cout << rhs;
    lf.pending << rhs;
  }
  return lf;
}
inline Logfile & operator<<(Logfile & lf, std::ostream & (*op)(std::ostream &)) {
  if (lf.logfile.is_open()) {
    std::cout << op << std::flush;
    lf.logfile << op;
  } else if (lf.mainProcess) {
    std::cout << op << std::flush;
    lf.pending << op;
  }
  return lf;
}
extern Logfile hlog;

/// Like hlog, but only to the log file (not stdout).
class Logfile_only {};
template <typename Val>
Logfile_only & operator<<(Logfile_only & lf, Val const & rhs) {
  if (hlog.logfile.is_open()) {
    hlog.logfile << rhs;
  } else if (hlog.mainProcess) {
    hlog.pending << rhs;
  }
  return lf;
}
inline Logfile_only & operator<<(Logfile_only & lf, std::ostream & (*op)(std::ostream &)) {
  if (hlog.logfile.is_open()) {
    hlog.logfile << op;
  } else if (hlog.mainProcess) {
    hlog.pending << op;
  }
  return lf;
}
extern Logfile_only hlogfile;

}
#endif /* LOGFILE_H */
