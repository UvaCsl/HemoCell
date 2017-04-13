#ifndef CELLFIELDS3D_CPP
#define CELLFIELDS3D_CPP
#include "cellFields3D.h"


CellFields3D::CellFields3D( MultiBlockLattice3D<double, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth) :
  lattice(lattice_)
{   
  envelopeSize=particleEnvelopeWidth;
  pcout << "(CellField3D) particle envelope SHOULD BE larger then largest celldiameter [lu]: " << particleEnvelopeWidth << std::endl;

  MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
	MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );

  immersedParticles = new MultiParticleField3D<HEMOCELL_PARTICLE_FIELD>(
            particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
  immersedParticles->periodicity().toggleAll(false);
  immersedParticles->periodicity().toggle(0,true);
  immersedParticles->toggleInternalStatistics(false);
  InitAfterLoadCheckpoint();
}

HemoCellField * CellFields3D::addCellType(TriangularSurfaceMesh<double> & meshElement, double hematocrit, std::string name_)
{
  HemoCellField * cf = new HemoCellField(*this, Cell3D<double,DESCRIPTOR>(), meshElement);
  cf->hematocrit = hematocrit;
  cf->name = name_;
  cf->ctype = cellFields.size();
  cellFields.push_back(cf);
  return cf;
}

unsigned int CellFields3D::size() 
{
  return this->cellFields.size();
}

HemoCellField * CellFields3D::operator[](unsigned int index)
{
  return cellFields[index];
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

/*
 * Initialize variables that need to be loaded after a checkpoint AS WELL
 */
void CellFields3D::InitAfterLoadCheckpoint()
{
  std::vector<plint> const& blocks = immersedParticles->getLocalInfo().getBlocks();
  for (pluint iBlock=0; iBlock<blocks.size(); ++iBlock) {
    SmartBulk3D bulk(immersedParticles->getMultiBlockManagement(),blocks[iBlock]);
    Box3D blk = bulk.getBulk();
    immersedParticles->getComponent(blocks[iBlock]).setlocalDomain(blk);
    immersedParticles->getComponent(blocks[iBlock]).cellFields = this;
    immersedParticles->getComponent(blocks[iBlock]).atomicBlockId = blocks[iBlock];
    immersedParticles->getComponent(blocks[iBlock]).atomicLattice = &lattice.getComponent(blocks[iBlock]);
    immersedParticles->getComponent(blocks[iBlock]).envelopeSize = envelopeSize;
  }
}

void CellFields3D::load(XMLreader * documentXML, plint & iter)
{

    std::string outDir = global::directories().getOutputDir();
    std::string firstField = (*(documentXML->getChildren( documentXML->getFirstId() )[0])).getName();
    bool isCheckpointed = (firstField=="Checkpoint");
    if (isCheckpointed) {
    hemocellfunction = true;
        (*documentXML)["Checkpoint"]["General"]["Iteration"].read(iter);
        parallelIO::load(outDir + "lattice", lattice, true);
        parallelIO::load(outDir + "particleField", *immersedParticles, true);

        InitAfterLoadCheckpoint();
    hemocellfunction = false;
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

void CellFields3D::HemoInterpolateFluidVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  
  /*
    Box3D temp = blocks[0]->getBoundingBox();
    Dot3D tmp = blocks[0]->getLocation();
    cerr << "Box: " << temp.x0 << " " << temp.x1 << " " << temp.y0  << " "<< temp.y1  << " "<< temp.z0 <<" "<< temp.z1 <<std::endl;
    cerr << "Loc: " << tmp.x << " " << tmp.y << " " << tmp.z << std::endl;
    */
    
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->interpolateFluidVelocity(domain);
}
void CellFields3D::interpolateFluidVelocity() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoInterpolateFluidVelocity(),immersedParticles->getBoundingBox(),wrapper);
}

void CellFields3D::HemoSyncEnvelopes::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->syncEnvelopes();
}
void CellFields3D::syncEnvelopes() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  hemocellfunction=true;
  applyProcessingFunctional(new HemoSyncEnvelopes(),immersedParticles->getBoundingBox(),wrapper);
  hemocellfunction=false;
}

void CellFields3D::HemoAdvanceParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->advanceParticles();
}
void CellFields3D::advanceParticles() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoAdvanceParticles(),immersedParticles->getBoundingBox(),wrapper);
}

void CellFields3D::HemoSpreadParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->spreadParticleForce(domain);
}
void CellFields3D::spreadParticleForce() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSpreadParticleForce(),immersedParticles->getBoundingBox(),wrapper);
}

void CellFields3D::HemoApplyConstitutiveModel::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyConstitutiveModel();
}
void CellFields3D::applyConstitutiveModel() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoApplyConstitutiveModel(),immersedParticles->getBoundingBox(),wrapper);

}

void CellFields3D::HemoUnifyForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->unifyForceVectors();
}
void CellFields3D::unify_force_vectors() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoUnifyForceVectors(),immersedParticles->getBoundingBox(),wrapper);
}

void CellFields3D::HemoRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyRepulsionForce();
}
void CellFields3D::calculateRepulsionForce() {
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    applyProcessingFunctional(new HemoRepulsionForce(),immersedParticles->getBoundingBox(),wrapper);
}

void CellFields3D::HemoSeperateForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->separateForceVectors();
}
void CellFields3D::separate_force_vectors() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSeperateForceVectors(),immersedParticles->getBoundingBox(),wrapper);

}
void CellFields3D::HemoDeleteIncompleteCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->deleteIncompleteCells();
}
void CellFields3D::deleteIncompleteCells() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoDeleteIncompleteCells(),immersedParticles->getBoundingBox(),wrapper);

}

CellFields3D::HemoSeperateForceVectors * CellFields3D::HemoSeperateForceVectors::clone() const { return new CellFields3D::HemoSeperateForceVectors(*this);}
CellFields3D::HemoUnifyForceVectors *    CellFields3D::HemoUnifyForceVectors::clone() const    { return new CellFields3D::HemoUnifyForceVectors(*this);}
CellFields3D::HemoSpreadParticleForce *  CellFields3D::HemoSpreadParticleForce::clone() const { return new CellFields3D::HemoSpreadParticleForce(*this);}
CellFields3D::HemoInterpolateFluidVelocity * CellFields3D::HemoInterpolateFluidVelocity::clone() const { return new CellFields3D::HemoInterpolateFluidVelocity(*this);}
CellFields3D::HemoAdvanceParticles *     CellFields3D::HemoAdvanceParticles::clone() const { return new CellFields3D::HemoAdvanceParticles(*this);}
CellFields3D::HemoApplyConstitutiveModel * CellFields3D::HemoApplyConstitutiveModel::clone() const { return new CellFields3D::HemoApplyConstitutiveModel(*this);}
CellFields3D::HemoSyncEnvelopes *        CellFields3D::HemoSyncEnvelopes::clone() const { return new CellFields3D::HemoSyncEnvelopes(*this);}
CellFields3D::HemoRepulsionForce *        CellFields3D::HemoRepulsionForce::clone() const { return new CellFields3D::HemoRepulsionForce(*this);}
CellFields3D::HemoDeleteIncompleteCells *        CellFields3D::HemoDeleteIncompleteCells::clone() const { return new CellFields3D::HemoDeleteIncompleteCells(*this);}

void CellFields3D::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
   for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::dynamicVariables;
   }
}

MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & CellFields3D::getParticleField3D() { return *immersedParticles; };
CellFields3D::~CellFields3D() {
    delete immersedParticles;
}
#endif


