#ifndef FCN_CHECKPOINT_HH
#define FCN_CHECKPOINT_HH
#include "checkPoint.h"


template<typename T, template<typename U> class Descriptor>
void FcnCheckpoint<T, Descriptor>::init(XMLreader & xmlr) {
    std::string firstField = (*(xmlr.getChildren( xmlr.getFirstId() )[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    isCheckpointed = (firstField=="Checkpoint");

    XMLwriter& xmlMultiBlock = xmlw["Checkpoint"];
    xmlMultiBlock["General"]["Iteration"].set(0);

    if (!isCheckpointed) { copyXMLreader2XMLwriter(xmlr["ficsion"], xmlw["Checkpoint"]); }
    else { copyXMLreader2XMLwriter(xmlr["Checkpoint"]["ficsion"], xmlw["Checkpoint"]); }
}


template<typename T, template<typename U> class Descriptor>
void FcnCheckpoint<T, Descriptor>::load(std::string paramXmlFileName, MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<CellField3D<T, Descriptor>* > & cellFields, plint & iter) {
    XMLreader document(paramXmlFileName);
    load(document, lattice, cellFields, iter);
}


template<typename T, template<typename U> class Descriptor>
void FcnCheckpoint<T, Descriptor>::load(XMLreader & documentXML, MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<CellField3D<T, Descriptor>* > & cellFields, plint & iter) {
    std::string outDir = global::directories().getOutputDir();

    std::string firstField = (*(documentXML.getChildren( documentXML.getFirstId() )[0])).getName(); // VERY COMPLICATED! Hope I could find sth easier!
    isCheckpointed = (firstField=="Checkpoint");
    if (isCheckpointed) {
        documentXML["Checkpoint"]["General"]["Iteration"].read(iter);
        parallelIO::load(outDir + "lattice", lattice, true);
        for (int icf= 0; icf < cellFields.size(); ++icf) {
            parallelIO::load(outDir + cellFields[icf]->getIdentifier(), cellFields[icf]->getParticleField3D(), true);
            cellFields[icf]->createCellMap();
            cellFields[icf]->synchronizeCellQuantities();
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void FcnCheckpoint<T, Descriptor>::save(MultiBlockLattice3D<T, Descriptor> & lattice, std::vector<CellField3D<T, Descriptor>* > & cellFields, plint iter) {
    std::string outDir = global::directories().getOutputDir();
    /* Rename files, for safety reasons */
    renameFileToDotOld(outDir + "lattice.dat");
    renameFileToDotOld(outDir + "lattice.plb");
    renameFileToDotOld(outDir + "checkpoint.xml");
    for (int icf= 0; icf < cellFields.size(); ++icf) {
        renameFileToDotOld(outDir + cellFields[icf]->getIdentifier() + ".dat");
        renameFileToDotOld(outDir + cellFields[icf]->getIdentifier() + ".plb");
    }
    /* Save */
    parallelIO::save(lattice, "lattice", true);
    for (int icf= 0; icf < cellFields.size(); ++icf) {
        parallelIO::save(cellFields[icf]->getParticleField3D(), cellFields[icf]->getIdentifier(), true);
    }
    // Upon success, save xml and rename files!
    xmlw["Checkpoint"]["General"]["Iteration"].set(iter);
    xmlw.print(outDir + "checkpoint.xml");
}


/* ******************* copyXMLreader2XMLwriter ***************************************** */
void copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer) {
    std::string text = reader.getFirstText();
    if (!text.empty()) {
        writer[reader.getName()].setString(text);
    }
    std::vector<XMLreader*> const& children = reader.getChildren( reader.getFirstId() );
    for (pluint iNode=0; iNode<children.size(); ++iNode) {
        copyXMLreader2XMLwriter(*(children[iNode]), writer[reader.getName()]);
    }
}

void copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer) {
    std::string text;
    readerProxy.read(text);
    if (!text.empty()) {
        writer[readerProxy.getName()].setString(text);
    }
    std::vector<XMLreader*> const& children = readerProxy.getChildren();
    for (pluint iNode=0; iNode<children.size(); ++iNode) {
        copyXMLreader2XMLwriter(*(children[iNode]), writer[readerProxy.getName()]);
    }
}


inline bool file_exists (const std::string& name) {
    /* Checks if a file exists */

// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


int renameFileToDotOld(std::string fName) {
    plint renameStatus = 0;
    if (file_exists(fName)) {
        plint renameStatus = rename(fName.c_str(), (fName + ".old").c_str());
        if (renameStatus != 0) {
            pcout << fName << " error." << std::endl;
        }
    }
    return renameStatus;
}


#endif // FCN_CHECKPOINT_HH
