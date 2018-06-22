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


#include "profiler.h"
#include "hemocell_internal.h"

Profiler::Profiler(std::string name_) :
name(name_), parent(*this), current(this)
{}

Profiler::Profiler(std::string name_, Profiler & parent_) :
name(name_), parent(parent_)
{}

void Profiler::start() {
  if (started) {
    hemo::hlog << "(Profiler) (Warning) Timer " << name << " has already been started" << std::endl;
  } else {
    //Not root node and parent not started
    if (&parent != this && !parent.started) {
      hemo::hlog << "(Profiler) (Warning) Starting timer " << name << " but parent has not been started, starting it for now but you should fix this in the code " << std::endl;
      parent.start();
    }
    start_time = std::chrono::high_resolution_clock::now();
    started = true;
  }
  
  //Check siblings for started
  for (std::pair<const std::string,Profiler> & timer_pair : parent.timers) {
    Profiler & timer = timer_pair.second;
    if (&timer == this) {continue;}
    if (timer.started) {
      hemo::hlog << "(Profiler) (Warning) Starting timer " << name << " but sibling " << timer.name << " has also been started, starting it for now but you should fix this in the code " << std::endl;
    }
  }
  
  //Set current
  current = this;
  Profiler * root = this;
  do  {
    root = &root->parent;
    root->current = this;
  } while (root != &root->parent);  
}

void Profiler::stop_nowarn() {
  if (started)
  {
    std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
    total_time = total_time + (stop_time - start_time);
    started = false;
  }
  
  //Stop all child timers
  for (std::pair<const std::string,Profiler> & timer_pair : timers) {
    Profiler & timer = timer_pair.second;
    timer.stop_nowarn();
  }
}

void Profiler::stop() {
  if (!started) {
    hemo::hlog << "(Profiler) (Warning) Timer " << name << " has not been started" << std::endl;
  } else {
    std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
    total_time = total_time + (stop_time - start_time);
    started = false;
  }
  
  //Stop all child timers
  for (std::pair<const std::string,Profiler> & timer_pair : timers) {
    Profiler & timer = timer_pair.second;
    timer.stop_nowarn();
  }
  
  //Adjust current timer
  current = &parent;
  Profiler * root = this;
  do  {
    root = &root->parent;
    root->current = &parent;
  } while (root != &root->parent);  
}

void Profiler::reset() {
  started = false;
  total_time = std::chrono::high_resolution_clock::duration::zero();
  
  //Reset all child timers
  for (std::pair<const std::string,Profiler> & timer_pair : timers) {
    Profiler & timer = timer_pair.second;
    timer.reset();
  }
}

std::chrono::high_resolution_clock::duration Profiler::elapsed() {
  if (!started) {
    return total_time;
  } else {
    return total_time + (std::chrono::high_resolution_clock::now() - start_time);
  }
}

std::string Profiler::elapsed_string() {
  if (!started) {
    return std::to_string(((double)std::chrono::duration_cast<std::chrono::milliseconds>(total_time).count())/1000.0);
  } else {
    return std::to_string(((double)std::chrono::duration_cast<std::chrono::milliseconds>(total_time).count())/1000.0);
  }
}

template<typename T>
void Profiler::printStatistics_inner(int level, T & out) {
  out << std::string(level,' ') << name << ": " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(total_time + (std::chrono::high_resolution_clock::now() - start_time)).count()/1000.0) << std::endl;
  //Print all child timers
  for (std::pair<const std::string,Profiler> & timer_pair : timers) {
    Profiler & timer = timer_pair.second;
    timer.printStatistics_inner(level+1, out);
  }
}

void Profiler::printStatistics() {
  hemo::hlog << "Hemocell Profiler Statistics (Only Process 0):" << std::endl;
  printStatistics_inner(0, hemo::hlog);
}

void Profiler::outputStatistics() {
  fstream sout;
  bool opened = false;
  int turn = 0;
  
  while (turn < global::mpi().getSize()) {
    if (turn == global::mpi().getRank()) {
      sout.open(hlog.filename + ".statistics",std::fstream::app);
      if (!sout.is_open()) {
        std::cout << "(Profiler) (Error) Opening " + hlog.filename << ".statistics, outputting everything to logfile instead" << std::endl;
        hemo::hlog << "Process " << global::mpi().getRank() << ":" << std::endl;
        printStatistics_inner(1,hlog.logfile);
      } else {
        opened = true;
        sout << "Process " << global::mpi().getRank() << ":" << std::endl;
        printStatistics_inner(1,sout);
      }
      if (opened) {
        sout.close();
      }
    }
    global::mpi().barrier();
    turn++;
  }  
}


Profiler & Profiler::operator[] (std::string name) {
  //try_emplace from c++17 would be sooo nice here, but gcc 5.4 does not yet have it
  if (timers.find(name) == timers.end()) {
    timers.insert(std::pair<std::string,Profiler>(name,Profiler(name, *this)));
  }
  return timers.at(name);
};

Profiler & Profiler::getCurrent() {
  if (this != &parent) {
    hemo::hlog << "(Profiler) (Warning) getCurrent called from non-root Profiler object, this will probably be incorrect" << std::endl;
  }
  return *current;
}