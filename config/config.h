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

#include "external/tinyxml2/tinyxml2.h"
#include <string>
#include <iostream>
#include <sstream>

namespace hemo {

  class XMLElement {
    tinyxml2::XMLElement * orig;
  public:
    XMLElement(tinyxml2::XMLElement * orig_) : orig(orig_){}
    XMLElement operator[] (std::string name) const;
    
    template<typename T>
    T read() {
      std::stringstream value(orig->GetText());
      T ret = T();
      if (!(value>>ret)) {
        std::cout << "Cannot convert value from XML element" << std::endl;
      }
      return ret;
    }
    
    tinyxml2::XMLElement * getOrig() { return orig;}
    
  };
  
  
  class Config {
    tinyxml2::XMLDocument * orig;
  public:
    bool checkpointed;

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

void loadDirectories(std::string configFileName, hemo::Config * cfg);

struct ConfigValues {
    bool cellsDeletedInfo = false;
};

extern ConfigValues globalConfigValues;

void loadGlobalConfigValues(hemo::Config * cfg);

}
#endif
