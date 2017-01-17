#include "config.h"

Config::Config(std::string paramXmlFileName) : plb::XMLreader(paramXmlFileName)
{ 

 // Check if it is a fresh start or a checkpointed run 
 std::string firstField = (*(this->getChildren(this->getFirstId())[0])).getName(); 

 if (firstField == "hemocell") { checkpointed = 0; } //If the first field is not checkpoint but hemocell
 else { checkpointed = 1; }


}

/*
 * Overload the overloaded operator to provide convenient access to
 * this["parameters"]["etc"]
 * Also hide the read() semantics, so it either returns a 
 */
plb::XMLreaderProxy Config::operator[] (std::string name) const
{
 // Set our direct access
 if (checkpointed) {
   return plb::XMLreader::operator[]("Checkpoint")["hemocell"][name];
 } else {
   return plb::XMLreader::operator[]("hemocell")[name];
 }
}
