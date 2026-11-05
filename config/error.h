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
#ifndef HEMOCELL_ERROR_H
#define HEMOCELL_ERROR_H

#include <sstream>
#include <stdexcept>
#include <string>

namespace hemo {

/// Thrown by fatal() instead of ending the job while a FatalErrorsThrow
/// object exists.
class FatalError : public std::runtime_error {
public:
  explicit FatalError(const std::string & message) : std::runtime_error(message) {}
};

/**
 * @brief Report a fatal error and end the whole MPI job.
 *
 * Prints the message, prefixed with the MPI rank, to stderr on the calling
 * rank (`hlog` writes on rank 0 only, so an error found on another rank would
 * otherwise go unreported), adds it to the log file if this rank writes one,
 * and calls `MPI_Abort(MPI_COMM_WORLD, 1)` (`std::exit(1)` when MPI is not
 * initialised). Not collective: a single rank may call it.
 *
 * While a FatalErrorsThrow object exists it throws FatalError instead, so
 * tests can check error paths. Only use that where every rank reaches the
 * error, otherwise the other ranks wait forever in their next collective.
 *
 * @param[in] message Description of the error, without a trailing newline.
 */
[[noreturn]] void fatal(const std::string & message);

/// Report a warning from any rank: like fatal() it writes to stderr with the
/// rank (and to the log file on rank 0), but the run continues.
void warning(const std::string & message);

namespace detail {
inline void streamAll(std::ostringstream &) {}
template <typename First, typename... Rest>
void streamAll(std::ostringstream & out, const First & first, const Rest &... rest) {
  out << first;
  streamAll(out, rest...);
}
} // namespace detail

/// fatal() with its arguments streamed into the message,
/// e.g. `fatal("cell type ", index, " does not exist")`.
template <typename First, typename Second, typename... Rest>
[[noreturn]] void fatal(const First & first, const Second & second, const Rest &... rest) {
  std::ostringstream message;
  detail::streamAll(message, first, second, rest...);
  fatal(message.str());
}

/// warning() with its arguments streamed into the message.
template <typename First, typename Second, typename... Rest>
void warning(const First & first, const Second & second, const Rest &... rest) {
  std::ostringstream message;
  detail::streamAll(message, first, second, rest...);
  warning(message.str());
}

/// While an object of this class exists, fatal() throws FatalError instead of
/// ending the job (for tests of error paths). Not thread safe.
class FatalErrorsThrow {
public:
  FatalErrorsThrow();
  ~FatalErrorsThrow();
  FatalErrorsThrow(const FatalErrorsThrow &) = delete;
  FatalErrorsThrow & operator=(const FatalErrorsThrow &) = delete;
};

} // namespace hemo
#endif
