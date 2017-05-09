#ifndef CELLFIELDS3D_CPP
#define CELLFIELDS3D_CPP
#include "hemoCellFields.h"


HemoCellFields::HemoCellFields( MultiBlockLattice3D<double, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth) :
  lattice(lattice_)
{   
  envelopeSize=particleEnvelopeWidth;
  pcout << "(Hemocell) (HemoCellFields) (Init) particle envelope: " << particleEnvelopeWidth << " [lu]" << std::endl;

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

HemoCellField * HemoCellFields::addCellType(TriangularSurfaceMesh<double> & meshElement, std::string name_)
{
  HemoCellField * cf = new HemoCellField(*this, meshElement);
  cf->name = name_;
  cf->ctype = cellFields.size();
  cellFields.push_back(cf);
  return cf;
}

unsigned int HemoCellFields::size() 
{
  return this->cellFields.size();
}

HemoCellField * HemoCellFields::operator[](unsigned int index)
{
  return cellFields[index];
}

HemoCellField * HemoCellFields::operator[](string name)
{
  for (unsigned int i = 0; i < cellFields.size(); i++) {
      if (cellFields[i]->name == name) {
          return cellFields[i];
      }
  } 
  return NULL;
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
void HemoCellFields::copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer) {
    std::string text = reader.getFirstText();
    if (!text.empty()) {
        writer[reader.getName()].setString(text);
    }
    std::vector<XMLreader*> const& children = reader.getChildren( reader.getFirstId() );
    for (pluint iNode=0; iNode<children.size(); ++iNode) {
        copyXMLreader2XMLwriter(*(children[iNode]), writer[reader.getName()]);
    }
}

void HemoCellFields::copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer) {
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

void HemoCellFields::setParticleUpdateScheme(double _cellTimeStep) 
{
  cellTimeStep = _cellTimeStep;
}

/*
 * Initialize variables that need to be loaded after a checkpoint AS WELL
 */
void HemoCellFields::InitAfterLoadCheckpoint()
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
    
    //Calculate neighbours 
    immersedParticles->getSparseBlockStructure().findNeighbors(blocks[iBlock], envelopeSize,
                           immersedParticles->getComponent(blocks[iBlock]).neighbours);
  if (immersedParticles->getComponent(blocks[iBlock]).neighbours.size() > HEMOCELL_MAX_NEIGHBOURS) {
    cerr << "More neighbours of atomic block than allowed, exiting ...";
    exit(0);
  }
  }
}

void HemoCellFields::load(XMLreader * documentXML, unsigned int & iter)
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

void HemoCellFields::save(XMLreader * documentXML, unsigned int iter)
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

void HemoCellFields::HemoInterpolateFluidVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  
  /*
    Box3D temp = blocks[0]->getBoundingBox();
    Dot3D tmp = blocks[0]->getLocation();
    cerr << "Box: " << temp.x0 << " " << temp.x1 << " " << temp.y0  << " "<< temp.y1  << " "<< temp.z0 <<" "<< temp.z1 <<std::endl;
    cerr << "Loc: " << tmp.x << " " << tmp.y << " " << tmp.z << std::endl;
    */
    
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->interpolateFluidVelocity(domain);
}
void HemoCellFields::interpolateFluidVelocity() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyTimedProcessingFunctional(new HemoInterpolateFluidVelocity(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoSyncEnvelopes::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->syncEnvelopes();
}
void HemoCellFields::syncEnvelopes() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  hemocellfunction=true;
  applyTimedProcessingFunctional(new HemoSyncEnvelopes(),immersedParticles->getBoundingBox(),wrapper);
  hemocellfunction=false;
}

void HemoCellFields::HemoAdvanceParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->advanceParticles();
}
void HemoCellFields::advanceParticles() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyTimedProcessingFunctional(new HemoAdvanceParticles(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoSpreadParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->spreadParticleForce(domain);
}
void HemoCellFields::spreadParticleForce() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyTimedProcessingFunctional(new HemoSpreadParticleForce(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoApplyConstitutiveModel::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyConstitutiveModel();
}
void HemoCellFields::applyConstitutiveModel() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyTimedProcessingFunctional(new HemoApplyConstitutiveModel(),immersedParticles->getBoundingBox(),wrapper);

}

void HemoCellFields::HemoUnifyForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->unifyForceVectors();
}
void HemoCellFields::unify_force_vectors() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoUnifyForceVectors(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyRepulsionForce();
}
void HemoCellFields::calculateRepulsionForce() {
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    applyTimedProcessingFunctional(new HemoRepulsionForce(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoSeperateForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->separateForceVectors();
}
void HemoCellFields::separate_force_vectors() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSeperateForceVectors(),immersedParticles->getBoundingBox(),wrapper);

}
void HemoCellFields::HemoDeleteIncompleteCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->deleteIncompleteCells();
}
void HemoCellFields::deleteIncompleteCells() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoDeleteIncompleteCells(),immersedParticles->getBoundingBox(),wrapper);

}

HemoCellFields::HemoSeperateForceVectors * HemoCellFields::HemoSeperateForceVectors::clone() const { return new HemoCellFields::HemoSeperateForceVectors(*this);}
HemoCellFields::HemoUnifyForceVectors *    HemoCellFields::HemoUnifyForceVectors::clone() const    { return new HemoCellFields::HemoUnifyForceVectors(*this);}
HemoCellFields::HemoSpreadParticleForce *  HemoCellFields::HemoSpreadParticleForce::clone() const { return new HemoCellFields::HemoSpreadParticleForce(*this);}
HemoCellFields::HemoInterpolateFluidVelocity * HemoCellFields::HemoInterpolateFluidVelocity::clone() const { return new HemoCellFields::HemoInterpolateFluidVelocity(*this);}
HemoCellFields::HemoAdvanceParticles *     HemoCellFields::HemoAdvanceParticles::clone() const { return new HemoCellFields::HemoAdvanceParticles(*this);}
HemoCellFields::HemoApplyConstitutiveModel * HemoCellFields::HemoApplyConstitutiveModel::clone() const { return new HemoCellFields::HemoApplyConstitutiveModel(*this);}
HemoCellFields::HemoSyncEnvelopes *        HemoCellFields::HemoSyncEnvelopes::clone() const { return new HemoCellFields::HemoSyncEnvelopes(*this);}
HemoCellFields::HemoRepulsionForce *        HemoCellFields::HemoRepulsionForce::clone() const { return new HemoCellFields::HemoRepulsionForce(*this);}
HemoCellFields::HemoDeleteIncompleteCells *        HemoCellFields::HemoDeleteIncompleteCells::clone() const { return new HemoCellFields::HemoDeleteIncompleteCells(*this);}

void HemoCellFields::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
   for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::dynamicVariables;
   }
}

MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & HemoCellFields::getParticleField3D() { return *immersedParticles; };
HemoCellFields::~HemoCellFields() {
    delete immersedParticles;
}
#endif


