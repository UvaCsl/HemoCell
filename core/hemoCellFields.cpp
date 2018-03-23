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
#include "hemoCellFields.h"
#include "hemocell.h"


HemoCellFields::HemoCellFields( MultiBlockLattice3D<T, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth, HemoCell & hemocell_) :
  lattice(&lattice_), hemocell(hemocell_)
{   
  envelopeSize=particleEnvelopeWidth;
  pcout << "(Hemocell) (HemoCellFields) (Init) particle envelope: " << particleEnvelopeWidth << " [lu]" << std::endl;
  
  createParticleField();
}

void HemoCellFields::createParticleField(SparseBlockStructure3D* sbStructure, ThreadAttribution * tAttribution) {
  if (!sbStructure) {
    sbStructure = lattice->getSparseBlockStructure().clone();
  }
  if (!tAttribution) {
    tAttribution = lattice->getMultiBlockManagement().getThreadAttribution().clone();
  }
  plint refinement = lattice->getMultiBlockManagement().getRefinementLevel();

  immersedParticles = new MultiParticleField3D<HEMOCELL_PARTICLE_FIELD>(MultiBlockManagement3D(
      *sbStructure,
      tAttribution,
      envelopeSize,
      refinement ), defaultMultiBlockPolicy3D().getCombinedStatistics() );
  InitAfterLoadCheckpoint();

  immersedParticles->periodicity().toggle(0,lattice->periodicity().get(0));
  immersedParticles->periodicity().toggle(1,lattice->periodicity().get(1));
  immersedParticles->periodicity().toggle(2,lattice->periodicity().get(2));

  immersedParticles->toggleInternalStatistics(false);
  
  InitAfterLoadCheckpoint();
}

HemoCellField * HemoCellFields::addCellType(TriangularSurfaceMesh<T> & meshElement, std::string name_)
{
  HemoCellField * cf = new HemoCellField(*this, meshElement, name_, cellFields.size());
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

bool error_shown = false;

/*
 * Initialize variables that need to be loaded after a checkpoint AS WELL
 */
void HemoCellFields::InitAfterLoadCheckpoint()
{
  number_of_cells = getTotalNumberOfCells(*this);
  std::vector<plint> const& blocks = immersedParticles->getLocalInfo().getBlocks();
  for (pluint iBlock=0; iBlock<blocks.size(); ++iBlock) {
    SmartBulk3D bulk(immersedParticles->getMultiBlockManagement(),blocks[iBlock]);
    Box3D blk = bulk.getBulk();
    immersedParticles->getComponent(blocks[iBlock]).setlocalDomain(blk);
    immersedParticles->getComponent(blocks[iBlock]).cellFields = this;
    immersedParticles->getComponent(blocks[iBlock]).atomicBlockId = blocks[iBlock];
    immersedParticles->getComponent(blocks[iBlock]).atomicLattice = &lattice->getComponent(blocks[iBlock]);
    immersedParticles->getComponent(blocks[iBlock]).envelopeSize = envelopeSize;
    
    BlockLattice3D<T,DESCRIPTOR> * fluid = immersedParticles->getComponent(blocks[iBlock]).atomicLattice;
    for(unsigned int x = 0; x < fluid->getNx(); x++ ) {
      for(unsigned int y = 0; y < fluid->getNy(); y++ ) {
        for(unsigned int z = 0; z < fluid->getNz(); z++ ) {
          if (!fluid->get(x,y,z).getDynamics().isBoundary()) {
            immersedParticles->getComponent(blocks[iBlock]).nFluidCells++;
          }
        }
      }
    }
    
    //Calculate neighbours 
    immersedParticles->getSparseBlockStructure().findNeighbors(blocks[iBlock], envelopeSize,
                           immersedParticles->getComponent(blocks[iBlock]).neighbours);
  
    if (immersedParticles->getComponent(blocks[iBlock]).neighbours.size() > 27 && !error_shown) {
      error_shown = true;
      cerr << "(HemoCell) WARNING: The number of atomic neighbours is suspiciously high: " << immersedParticles->getComponent(blocks[iBlock]).neighbours.size() << " Usually it should be < 27 ! Check the atomic block structure!\n";
    }
  }
}

void HemoCellFields::load(XMLreader * documentXML, unsigned int & iter, Config * cfg)
{

    std::string outDir = global::directories().getOutputDir();
    if (cfg) {
      try {
        outDir = (*cfg)["parameters"]["checkpointDirectory"].read<string>() + "/";
        if (outDir[0] != '/') {
          outDir = "./" + outDir;
        }
      } catch (std::invalid_argument & exeption) {
      }
    }
    std::string firstField = (*(documentXML->getChildren( documentXML->getFirstId() )[0])).getName();
    bool isCheckpointed = (firstField=="Checkpoint");
    if (isCheckpointed) {
      (*documentXML)["Checkpoint"]["General"]["Iteration"].read(iter);
      parallelIO::load(outDir + "lattice", *lattice, true);
      parallelIO::load(outDir + "particleField", *immersedParticles, true);
      
    } else {
      pcout << "(HemoCell) (CellFields) loading checkpoint from non-checkpoint Config" << endl;
      parallelIO::load(outDir + "lattice", *lattice, true);
      parallelIO::load(outDir + "particleField", *immersedParticles, true);      
      
    }
    
    InitAfterLoadCheckpoint();
    syncEnvelopes();
    deleteIncompleteCells();
}

void HemoCellFields::save(XMLreader * documentXML, unsigned int iter, Config * cfg)
{
    XMLreader xmlr = *documentXML;
    XMLwriter xmlw;
    std::string firstField = (*(xmlr.getChildren( xmlr.getFirstId() )[0])).getName(); 
    bool isCheckpointed = (firstField=="Checkpoint");
    if (!isCheckpointed) { copyXMLreader2XMLwriter(xmlr["hemocell"], xmlw["Checkpoint"]); }
    else { copyXMLreader2XMLwriter(xmlr["Checkpoint"]["hemocell"], xmlw["Checkpoint"]); }

    std::string outDir = global::directories().getOutputDir();
    if (cfg) {
      try {
        outDir = (*cfg)["parameters"]["checkpointDirectory"].read<string>() + "/";
        if (outDir[0] != '/') {
          outDir = "./" + outDir;
        }
      } catch (std::invalid_argument & exeption) {
      }
    }
    mkpath(outDir.c_str(), 0777);

    
    /* Rename files, for safety reasons */
    if (global::mpi().isMainProcessor()) {
        renameFileToDotOld(outDir + "lattice.dat");
        renameFileToDotOld(outDir + "lattice.plb");
        renameFileToDotOld(outDir + "cellfields.dat");
        renameFileToDotOld(outDir + "cellfields.plb");
        renameFileToDotOld(outDir + "particleField.dat");
        renameFileToDotOld(outDir + "particleField.plb");
        renameFileToDotOld(outDir + "checkpoint.xml");
    } 
    
    global::mpi().barrier();
    
    /* Save XML & Data */
    xmlw["Checkpoint"]["General"]["Iteration"].set(iter);
    xmlw.print(outDir + "checkpoint.xml");
    parallelIO::save(*lattice, outDir + "lattice", true);
    parallelIO::save(*immersedParticles, outDir + "particleField", true);
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
  if (hemocell.iter % particleVelocityUpdateTimescale == 0) {

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyTimedProcessingFunctional(new HemoInterpolateFluidVelocity(),immersedParticles->getBoundingBox(),wrapper);
  
  }
}

void HemoCellFields::HemoSyncEnvelopes::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->syncEnvelopes();
}
void HemoCellFields::syncEnvelopes() {
  if (hemocell.iter % particleVelocityUpdateTimescale == 0) {
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(immersedParticles);
    applyTimedProcessingFunctional(new HemoSyncEnvelopes(),immersedParticles->getBoundingBox(),wrapper);
  }
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
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyConstitutiveModel(forced);
}
void HemoCellFields::applyConstitutiveModel(bool forced) {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  HemoApplyConstitutiveModel * fnct = new HemoApplyConstitutiveModel();
  fnct->forced = forced;
  applyTimedProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
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
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyRepulsionForce(forced);
}
void HemoCellFields::applyRepulsionForce(bool forced) {
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    HemoRepulsionForce * fnct = new HemoRepulsionForce();
    fnct->forced = forced;
    applyTimedProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
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
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->deleteIncompleteCells(verbose);
}
void HemoCellFields::deleteIncompleteCells(bool verbose) {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  HemoDeleteIncompleteCells * fnct = new HemoDeleteIncompleteCells();
  fnct->verbose = verbose;
  applyTimedProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);

}
void HemoCellFields::HemoGetParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD * pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  Box3D localDomain;
  intersect(domain,pf->localDomain,localDomain);
  pf->findParticles(localDomain,particles);
}
void HemoCellFields::getParticles(vector<HemoCellParticle *> & particles, Box3D& domain) {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoGetParticles(particles),domain,wrapper);
}
void HemoCellFields::HemoSetParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD * pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  Box3D localDomain;
  intersect(domain,pf->localDomain,localDomain);
  for (HemoCellParticle & particle : particles ) {
    pf->addParticle(localDomain,&particle);
  }
}
void HemoCellFields::addParticles(vector<HemoCellParticle> & particles) {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSetParticles(particles),immersedParticles->getBoundingBox(),wrapper);
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
HemoCellFields::HemoGetParticles *        HemoCellFields::HemoGetParticles::clone() const { return new HemoCellFields::HemoGetParticles(*this);}
HemoCellFields::HemoSetParticles *        HemoCellFields::HemoSetParticles::clone() const { return new HemoCellFields::HemoSetParticles(*this);}


void HemoCellFields::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
   for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::dynamicVariables;
   }
}

MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & HemoCellFields::getParticleField3D() { return *immersedParticles; };
HemoCellFields::~HemoCellFields() {
    delete immersedParticles;
}



