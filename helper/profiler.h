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

#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <string>
#include <map>
#include <logfile.h>

class Profiler {
public:
  Profiler(std::string name_);
  Profiler(std::string name_, Profiler & parent_);

  void start();
  void stop();
  void reset();
  void printStatistics();
  void outputStatistics();
  
  std::chrono::high_resolution_clock::duration elapsed();
  std::string elapsed_string();
  Profiler & operator[] (std::string);
  Profiler & getCurrent();


private:
  void stop_nowarn();
  template<typename T>
  void printStatistics_inner(int level, T & out);
  std::chrono::high_resolution_clock::duration total_time = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
  bool started = false;
  const std::string name;
  std::map<std::string,Profiler> timers;
  Profiler & parent;
  Profiler * current = this;
};

#endif /* PROFILER_H */

