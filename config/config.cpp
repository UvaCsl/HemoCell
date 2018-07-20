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
#include "config.h"
#include "logfile.h"
#include "genericFunctions.h"
#include <stdexcept>
#include <sys/stat.h>

#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"

namespace hemo {
  Config::Config(string paramXmlFileName) 
  {  
    load(paramXmlFileName);
  }

  void Config::reload(std::string paramXmlFileName) {
    delete orig;
    load(paramXmlFileName);
  }
  
  void Config::load(std::string paramXmlFileName) {
    if (!file_exists(paramXmlFileName)) {
      plb::pcout << paramXmlFileName + " is not an existing config file, exiting ..." << std::endl;
      exit(1);
    }
    orig = new tinyxml2::XMLDocument();
    orig->LoadFile(paramXmlFileName.c_str());
   // Check if it is a fresh start or a checkpointed run 
    tinyxml2::XMLNode * first =orig->FirstChild()->NextSibling();
    if (first) {
      string firstField = first->Value(); 
      if (firstField == "Checkpoint") { checkpointed = 1; } //If the first field is not checkpoint but hemocell
      else { checkpointed = 0; }
    }
  }
  
  Config::~Config() {
    delete orig;
  }
  
  /*
   * Overload the overloaded operator to provide convenient access to
   * this["parameters"]["etc"]
   * Also hide the read() semantics, so it either returns a 
   */
  XMLElement  Config::operator[] (string name) const
  {
   // Set our direct access
   if (checkpointed) {
     return static_cast<XMLElement>(orig->FirstChildElement("Checkpoint"))["hemocell"][name];
   } else {
     return static_cast<XMLElement>(orig->FirstChildElement("hemocell"))[name];
   }
  }
  
  XMLElement XMLElement::operator[] (string name) const {
    tinyxml2::XMLElement * child = orig->FirstChildElement(name.c_str());
    if (child) {
      return static_cast<XMLElement>(child);
    } else {
      throw invalid_argument("XML child " + name + " does not exist.");
    }
  }

void loadDirectories(std::string configFileName, hemo::Config * cfg)  {
    //TODO This should be done through hemocell config, not some palabos global
  plb::global::directories().setOutputDir("./tmp/");
  plb::global::directories().setInputDir("./");
  plb::global::IOpolicy().activateParallelIO(true);
  plb::global::IOpolicy().setStlFilesHaveLowerBound(false);

  //setting outdir
  try {
    std::string outDir = (*cfg)["parameters"]["outputDirectory"].read<string>() + "/";
    if (outDir[0] != '/') {
          outDir = "./" + outDir;
    }
    plb::global::directories().setOutputDir(outDir);
  } catch (std::invalid_argument & exeption) {}
  mkpath((plb::global::directories().getOutputDir() + "/hdf5/").c_str(), 0777);

  //Setting logfile and logdir  
  plb::global::directories().setLogOutDir("./log/");
  try {
    std::string outDir = (*cfg)["parameters"]["logDirectory"].read<string>() + "/";
    if (outDir[0] != '/') {
          outDir = "./" + outDir;
    }
    plb::global::directories().setLogOutDir(outDir);
  } catch (std::invalid_argument & exeption) {}
  string logfilename = "logfile";
  try {
    logfilename = (*cfg)["parameters"]["logFile"].read<string>();
  } catch (std::invalid_argument & exeption) {}
  mkpath(plb::global::directories().getLogOutDir().c_str(), 0777);

  string base_filename =  plb::global::directories().getLogOutDir() + logfilename;
  string filename = base_filename;
  if (!file_exists(base_filename)) { 
    goto logfile_open_done; 
  }
  for (int i = 0; i < INT_MAX; i++) {
    filename = base_filename + "." + to_string(i);
    if (!file_exists(filename)) {
      goto logfile_open_done; 
    }
  }

  logfile_open_done:

  hemo::hlog.filename = filename;
  plb::global::mpi().barrier();
  if (plb::global::mpi().getRank() == 0) {
    hemo::hlog.logfile.open(filename , std::fstream::out);
    if (!hemo::hlog.logfile.is_open()) {
      plb::pcerr << "(HemoCell) (LogFile) Error opening logfile, exiting" << std::endl;
      exit(1);
    }
  }  
}

ConfigValues global;

void loadGlobalConfigValues(hemo::Config * cfg) {
  try {
   global.cellsDeletedInfo = (*cfg)["verbose"]["cellsDeletedInfo"].read<int>();
  } catch(std::invalid_argument & e) {}
  try {
   global.enableCEPACfield = (*cfg)["parameters"]["enableCEPACfield"].read<int>();
  } catch(std::invalid_argument & e) {}
  try {
   global.enableSolidifyMechanics = (*cfg)["parameters"]["enableSolidifyMechanics"].read<int>();
  } catch(std::invalid_argument & e) {}
}

}

