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
#include "constant_defaults.h"
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
    if (orig->LoadFile(paramXmlFileName.c_str()) != tinyxml2::XML_SUCCESS) {
      plb::pcout << "(HemoCell) (Config) Error parsing " << paramXmlFileName << ": " << orig->ErrorName();
      if (orig->GetErrorStr1()) { plb::pcout << " " << orig->GetErrorStr1(); }
      plb::pcout << " (line " << orig->GetErrorLineNum() << "), exiting ..." << std::endl;
      exit(1);
    }
    // Check if it is a fresh start (root <hemocell>) or a checkpoint (root
    // <Checkpoint> containing <hemocell>). The root is validated once here so
    // that operator[] does not need to check it on every access.
    tinyxml2::XMLElement * root = orig->RootElement();
    checkpointed = root && (string(root->Name()) == "Checkpoint");
    tinyxml2::XMLElement * hemocellRoot = checkpointed ? root->FirstChildElement("hemocell")
                                                       : orig->FirstChildElement("hemocell");
    if (!hemocellRoot) {
      plb::pcout << "(HemoCell) (Config) Error: " << paramXmlFileName << " has no <hemocell> root element"
                 << (checkpointed ? " inside <Checkpoint>" : "") << ", exiting ..." << std::endl;
      exit(1);
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
  
  bool XMLElement::parseValue(const std::string & text, bool & value) {
    std::string word;
    if (!parseValue(text, word)) {
      return false;
    }
    if (word == "1" || word == "true") { value = true; return true; }
    if (word == "0" || word == "false") { value = false; return true; }
    return false;
  }

  bool XMLElement::parseValue(const std::string & text, std::string & value) {
    const std::string whitespace = " \t\n\r\f\v";
    const std::string::size_type first = text.find_first_not_of(whitespace);
    if (first == std::string::npos) {
      return false;
    }
    const std::string::size_type last = text.find_last_not_of(whitespace);
    value = text.substr(first, last - first + 1);
    return true;
  }

  XMLElement XMLElement::operator[] (string name) const {
    tinyxml2::XMLElement * child = orig->FirstChildElement(name.c_str());
    if (child) {
      return static_cast<XMLElement>(child);
    } else {
      throw invalid_argument("XML child " + name + " does not exist.");
    }
  }

void loadDirectories(hemo::Config * cfg, bool edit_out_dir)  {
    //TODO This should be done through hemocell config, not some palabos global
  plb::global::directories().setInputDir("./");
  plb::global::IOpolicy().activateParallelIO(true);
  plb::global::IOpolicy().setStlFilesHaveLowerBound(false);

  //setting outdir
  if (edit_out_dir) {
    try {
      std::string outDir = (*cfg)["parameters"]["outputDirectory"].read<string>();
      while (outDir.size() > 1 && outDir[outDir.size()-1] == '/') {
        outDir.pop_back();
      }
      if (outDir[0] != '/') {
            outDir = "./" + outDir;
      }
      plb::global::directories().setOutputDir(outDir);
    } catch (std::invalid_argument & exeption) {
      plb::global::directories().setOutputDir("./tmp");
    }
    string filename = plb::global::directories().getOutputDir();
    if (!file_exists(filename)) {
      goto outdir_done;
    }
    for (int i = 0; i < INT_MAX; i++) {
      filename = plb::global::directories().getOutputDir() + "_" + to_string(i);
      if (!file_exists(filename)) {
        plb::global::directories().setOutputDir(filename);
        goto outdir_done; 
      }
    }
outdir_done:
    plb::global::directories().setOutputDir(plb::global::directories().getOutputDir() + "/");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (!plb::global::mpi().getRank()) {
    mkpath((plb::global::directories().getOutputDir() + "/hdf5/").c_str(), 0777);
  }
  //Setting logdir  
  try {
    std::string outDir = (*cfg)["parameters"]["logDirectory"].read<string>() + "/";
    outDir = plb::global::directories().getOutputDir() + "/" + outDir;
    plb::global::directories().setLogOutDir(outDir);
  } catch (std::invalid_argument & exeption) {
    plb::global::directories().setLogOutDir(plb::global::directories().getOutputDir() + "/log/");
  }

  //Setting logfile
  string logfilename = "logfile";
  try {
    logfilename = (*cfg)["parameters"]["logFile"].read<string>();
  } catch (std::invalid_argument & exeption) {}
  if (!plb::global::mpi().getRank()) {
    mkpath(plb::global::directories().getLogOutDir().c_str(), 0777);
  }
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
    if (!hemo::hlog.open(filename)) {
      plb::pcerr << "(HemoCell) (LogFile) Error opening logfile, exiting" << std::endl;
      exit(1);
    }
  }  
  
  //Setting CheckpointDirectory
  try { 
    hemo::global.checkpointDirectory = plb::global::directories().getOutputDir() + "/" + (*cfg)["parameters"]["checkpointDirectory"].read<string>() + "/";
  } catch (std::invalid_argument & exception) {
    hemo::global.checkpointDirectory = plb::global::directories().getOutputDir() + "/checkpoint/";
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
#ifndef SOLIDIFY_MECHANICS
   if (global.enableSolidifyMechanics) {
     hlog << "(Hemocell) (Config) Error EnableSolidifyMechanics is true but SOLIDIFY_MECHANICS (compile time) Is not defined" << std::endl;
     exit(1);
   }
#endif
  } catch(std::invalid_argument & e) {}
}

}

