#include "config.h"
#include <stdexcept>

namespace hemo {
  Config::Config(std::string paramXmlFileName) 
  {  
    orig = new tinyxml2::XMLDocument();
    orig->LoadFile(paramXmlFileName.c_str());
   // Check if it is a fresh start or a checkpointed run 
    tinyxml2::XMLNode * first =orig->FirstChild()->NextSibling();
    if (first) {
      std::string firstField = first->Value(); 
      if (firstField == "Checkpoint") { checkpointed = 1; } //If the first field is not checkpoint but hemocell
      else { checkpointed = 0; }
    }
  }

  /*
   * Overload the overloaded operator to provide convenient access to
   * this["parameters"]["etc"]
   * Also hide the read() semantics, so it either returns a 
   */
  XMLElement  Config::operator[] (std::string name) const
  {
   // Set our direct access
   if (checkpointed) {
     return static_cast<XMLElement>(orig->FirstChildElement("Checkpoint"))["hemocell"][name];
   } else {
     return static_cast<XMLElement>(orig->FirstChildElement("hemocell"))[name];
   }
  }
  
  XMLElement XMLElement::operator[] (std::string name) const {
    tinyxml2::XMLElement * child = orig->FirstChildElement(name.c_str());
    if (child) {
      return static_cast<XMLElement>(child);
    } else {
      throw std::invalid_argument("XML child " + name + " does not exist.");
    }
  }

  
}
