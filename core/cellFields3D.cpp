#ifndef CELLFIELDS3D_CPP
#define CELLFIELDS3D_CPP
#include "shellModel3D.h"
#include "cellFields3D.h"


CellFields3D::CellFields3D( MultiBlockLattice3D<double, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth) :
  lattice(lattice_)
{   
  pcout << "(CellField3D) particle envelope SHOULD BE larger then largest celldiameter [lu]: " << particleEnvelopeWidth << std::endl;

  MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
	MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );

  immersedParticles = new MultiParticleField3D<HEMOCELL_PARTICLE_FIELD>(
            particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
  immersedParticles->periodicity().toggleAll(true);
  immersedParticles->toggleInternalStatistics(false);
}

void CellFields3D::addCellType(TriangularSurfaceMesh<double> * meshElement, double hematocrit, ShellModel3D<double> *cellmodel, std::string name_)
{
  HemoCellField * cf = new HemoCellField(*this, Cell3D<double,DESCRIPTOR>(*meshElement));
  cf->hematocrit = hematocrit;
  cf->meshElement = meshElement;
  cf->name = name_;
  cellFields.push_back(*cf);
}

unsigned int CellFields3D::size() 
{
  return this->cellFields.size();
}

HemoCellField * CellFields3D::operator[](unsigned int index)
{
  return &(cellFields[index]);
}
//Initialization
//vector <CellField3D<double, DESCRIPTOR>> CellFields3D::getLegacyCellFieldsVector() {
//  vector <CellField3D<double, DESCRIPTOR>> cfv;
//  for (uint i = 0 ; i < cellFields.size() ; i++) {
//    cfv.push_back(CellField3D<double,DESCRIPTOR>(ce
     
//  }
//  return cfv;
//}

//SAVING FUNCTIONS
/* ******************* copyXMLreader2XMLwriter ***************************************** */
void CellFields3D::copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer) {
    std::string text = reader.getFirstText();
    if (!text.empty()) {
        writer[reader.getName()].setString(text);
    }
    std::vector<XMLreader*> const& children = reader.getChildren( reader.getFirstId() );
    for (pluint iNode=0; iNode<children.size(); ++iNode) {
        copyXMLreader2XMLwriter(*(children[iNode]), writer[reader.getName()]);
    }
}

void CellFields3D::copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer) {
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

void CellFields3D::setParticleUpdateScheme(double _cellTimeStep) 
{
  cellTimeStep = _cellTimeStep;
}

void CellFields3D::InitAfterLoadCheckpoint()
{
}

void CellFields3D::load(XMLreader * documentXML, plint & iter)
{
    std::string outDir = global::directories().getOutputDir();
    std::string firstField = (*(documentXML->getChildren( documentXML->getFirstId() )[0])).getName();
    bool isCheckpointed = (firstField=="Checkpoint");
    if (isCheckpointed) {
        (*documentXML)["Checkpoint"]["General"]["Iteration"].read(iter);
        parallelIO::load(outDir + "lattice", lattice, true);
        parallelIO::load(outDir + "particleField", *immersedParticles, true);

        InitAfterLoadCheckpoint();
    }
}

void CellFields3D::save(XMLreader * documentXML, plint iter)
{
    XMLreader xmlr = *documentXML;
    XMLwriter xmlw;
    std::string firstField = (*(xmlr.getChildren( xmlr.getFirstId() )[0])).getName(); 
    bool isCheckpointed = (firstField=="Checkpoint");
    if (!isCheckpointed) { copyXMLreader2XMLwriter(xmlr["hemocell"], xmlw["Checkpoint"]); }
    else { copyXMLreader2XMLwriter(xmlr["Checkpoint"]["hemocell"], xmlw["Checkpoint"]); }

    std::string outDir = global::directories().getOutputDir();
    /* Rename files, for safety reasons */
    if (global::mpi().isMainProcessor()) {
        renameFileToDotOld(outDir + "lattice.dat");
        renameFileToDotOld(outDir + "lattice.plb");
        renameFileToDotOld(outDir + "checkpoint.xml");
        renameFileToDotOld(outDir + "cellfields.dat");
        renameFileToDotOld(outDir + "cellfields.plb");
    } else {
        usleep(1000); // Sleep for 1000 milliseconds
    }
    /* Save XML & Data */
    xmlw["Checkpoint"]["General"]["Iteration"].set(iter);
    xmlw.print(outDir + "checkpoint.xml");
    parallelIO::save(lattice, "lattice", true);
    parallelIO::save(*immersedParticles, "particleField", true);
    // Upon success, save xml and rename files!

}

void readPositionsCellFields(std::string particlePosFile) {
}

void CellFields3D::advanceParticles() {
}

void CellFields3D::spreadForceIBM() {
}

void CellFields3D::interpolateVelocityIBM() {
}

void CellFields3D::applyConstitutiveModel() {
}

#endif


