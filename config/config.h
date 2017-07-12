#ifndef HEMOCELL_CONFIG_H
#define HEMOCELL_CONFIG_H

#include "external/tinyxml2/tinyxml2.h"
#include <string>
#include <iostream>

namespace hemo {

  class XMLElement {
    const tinyxml2::XMLElement * orig;
  public:
    XMLElement(const tinyxml2::XMLElement * orig_) : orig(orig_){}
    XMLElement operator[] (std::string name) const;
    
    template<typename T>
    T read() {
      std:: stringstream value(orig->GetText());
      T ret = T();
      if (!(value>>ret)) {
        std::cout << "Cannot convert value from XML element" << std::endl;
      }
      return ret;
    }
  };
  
  
  class Config {
    tinyxml2::XMLDocument * orig;
  public:
    bool checkpointed;

    Config(std::string paramXmlFilename); 
    
    hemo::XMLElement operator[] (std::string name) const;
    
    tinyxml2::XMLNode* ShallowClone(tinyxml2::XMLDocument* document) const;
    bool ShallowEqual(const tinyxml2::XMLNode* compare ) const;
    bool Accept( tinyxml2::XMLVisitor* visitor ) const;
  };

}
#endif
