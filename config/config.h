#ifndef HEMOCELL_CONFIG_H
#define HEMOCELL_CONFIG_H

#include "hemocell_internal.h"
#include "TINYXML_xmlIO.h"
#include <string.h>

class Config : public plb::XMLreader {
  public:

		bool checkpointed;

    Config(std::string paramXmlFilename); 

    plb::XMLreaderProxy operator[] (std::string name) const;

};


#endif
