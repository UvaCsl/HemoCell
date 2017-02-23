#if 0
#ifndef FCN_CHECKPOINT_H
#define FCN_CHECKPOINT_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellField3D.h"
#include "fcnGenericFunctions.h"
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <string>

using namespace std;
using namespace plb;

/* ******************* copyXMLreader2XMLwriter ***************************************** */
void copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer);
void copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer);



template<typename T, template<typename U> class Descriptor>
class FcnCheckpoint
{
public:
	FcnCheckpoint(XMLreader & documentXML)
	        { init(documentXML); } ;
	FcnCheckpoint(std::string paramXmlFileName)
	        { XMLreader document(paramXmlFileName); init(document); } ;
	~FcnCheckpoint() {} ;
	bool wasCheckpointed() { return isCheckpointed; } ;
	virtual void init(XMLreader & xmlr) ;
    void load(std::string paramXmlFileName, MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<HemoCellField* > & cellFields, plint & iter);
    void load(XMLreader & documentXML, MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<HemoCellField* > & cellFields, plint & iter);
	void save(MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<HemoCellField* > & cellFields, plint iter);

private:
    XMLwriter xmlw;
    bool isCheckpointed;
};



#include "checkPoint.hh"
#endif // FCN_CHECKPOINT_H
#endif