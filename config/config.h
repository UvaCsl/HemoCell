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
#ifndef HEMOCELL_CONFIG_H
#define HEMOCELL_CONFIG_H

#include "profiler.h"

#include "external/tinyxml2/tinyxml2.h"
#include "constant_defaults.h"
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace hemo {

  class XMLElement {
    tinyxml2::XMLElement * orig;
  public:
    XMLElement(tinyxml2::XMLElement * orig_) : orig(orig_){}
    XMLElement operator[] (std::string name) const;

    /**
     * @brief Read the text of this element as a value of type T.
     *
     * Numbers are parsed with stream extraction and must make up the whole
     * text (surrounding whitespace is ignored), e.g. "1e6" is rejected for an
     * integer instead of silently becoming 1. Booleans accept 0/1/true/false.
     * Strings return the whole text with surrounding whitespace removed.
     *
     * @throws std::invalid_argument if the element is empty or its text cannot
     * be converted. Callers that treat a value as optional catch this exception
     * (the same one that is thrown for a missing element) and use a default.
     */
    template<typename T>
    T read() {
      const char * text = orig->GetText();
      if (!text) {
        throw std::invalid_argument(std::string("XML element <") + orig->Name() + "> is empty.");
      }
      T ret = T();
      if (!parseValue(text, ret)) {
        throw std::invalid_argument(std::string("Cannot convert value \"") + text + "\" of XML element <" + orig->Name() + ">.");
      }
      return ret;
    }
    
    tinyxml2::XMLElement * getOrig() { return orig;}

  private:
    template<typename V>
    static bool parseValue(const std::string & text, V & value) {
      std::istringstream stream(text);
      if (!(stream >> value)) {
        return false;
      }
      stream >> std::ws;
      return stream.eof(); // nothing but whitespace may follow the value
    }
    static bool parseValue(const std::string & text, bool & value);
    static bool parseValue(const std::string & text, std::string & value);
  };
  
  
  class Config {
    tinyxml2::XMLDocument * orig;
  public:
    /// True if the root element is <Checkpoint> (a checkpoint.xml written by HemoCell).
    bool checkpointed = false;

    Config(std::string paramXmlFilename); 
    ~Config();
    
    void reload(std::string paramXmlFilename);
  
    hemo::XMLElement operator[] (std::string name) const;
    
    tinyxml2::XMLNode* ShallowClone(tinyxml2::XMLDocument* document) const;
    bool ShallowEqual(const tinyxml2::XMLNode* compare ) const;
    bool Accept( tinyxml2::XMLVisitor* visitor ) const;
  private:
    void load(std::string paramXmlFilename);
  };

void loadDirectories(hemo::Config * cfg, bool edit_out_dir = true);

struct ConfigValues {
  bool hemoCellInitialized = false; // Keep track since two hemocells cannot run at the same time, because of static variables
  bool cellsDeletedInfo = false;

  bool enableCEPACfield = false;

  bool enableSolidifyMechanics = false;

  bool enableInteriorViscosity = false;
  
  std::string checkpointDirectory = "./checkpoint/";

  Profiler statistics = Profiler("HemoCell");
};

extern ConfigValues global;

void loadGlobalConfigValues(hemo::Config * cfg);

}
#endif
