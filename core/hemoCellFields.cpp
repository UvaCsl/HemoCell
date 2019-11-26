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
#include <mpi.h>

#include "hemoCellFields.h"
#include "hemocell.h"
#include "readPositionsBloodCells.h"
#include "constantConversion.h"
#include "bindingField.h"

#include "palabos3D.h"
#include "palabos3D.hh"


namespace hemo {

 
HemoCellFields::HemoCellFields( MultiBlockLattice3D<T, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth, HemoCell & hemocell_) :
  lattice(&lattice_), hemocell(hemocell_)
{   
  envelopeSize=particleEnvelopeWidth;
  hlog << "(Hemocell) (HemoCellFields) (Init) particle envelope: " << particleEnvelopeWidth << " [lu]" << std::endl;
  if (hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth() < 2) {
    hlog << "(Hemocell) (ERROR) fluid envelope is less than 2, this will cause incorrect forces over the block boundaries" <<endl;
    exit(1);
  }
  createParticleField();
  if (global.enableCEPACfield) {
    createCEPACfield();
  } 
  InitAfterLoadCheckpoint();
}

HemoCellFields::~HemoCellFields() {
  if (CEPACfield) {
    delete CEPACfield;
  }
  if (immersedParticles) {
    delete immersedParticles;
  }
  for (HemoCellField * field : cellFields) {
    delete field;
  }
  if (large_communicator) {
    delete large_communicator;
  }  
}

void HemoCellFields::createParticleField(SparseBlockStructure3D* sbStructure, ThreadAttribution * tAttribution) {
  bool del_sbStruct = false;
  if (!sbStructure) {
    sbStructure = hemocell.domain_lattice->getSparseBlockStructure().clone();
    del_sbStruct = true;
  }
  if (!tAttribution) {
    tAttribution = hemocell.domain_lattice->getMultiBlockManagement().getThreadAttribution().clone();
  }
  plint refinement = lattice->getMultiBlockManagement().getRefinementLevel();

  if (hemocell.preInlet) {
    preinlet_immersedParticles = new MultiParticleField3D<HEMOCELL_PARTICLE_FIELD>(MultiBlockManagement3D(
      *hemocell.preinlet_lattice->getSparseBlockStructure().clone(),
      hemocell.preinlet_lattice->getMultiBlockManagement().getThreadAttribution().clone(),
      envelopeSize,
      refinement ), plb::defaultMultiBlockPolicy3D().getCombinedStatistics() );
  }
  domain_immersedParticles = new MultiParticleField3D<HEMOCELL_PARTICLE_FIELD>(MultiBlockManagement3D(
      *sbStructure,
      tAttribution,
      envelopeSize,
      refinement ), plb::defaultMultiBlockPolicy3D().getCombinedStatistics() );
  if (hemocell.partOfpreInlet) {
    immersedParticles = preinlet_immersedParticles;
  } else {
    immersedParticles = domain_immersedParticles;
  }
  InitAfterLoadCheckpoint();

  immersedParticles->periodicity().toggle(0,lattice->periodicity().get(0));
  immersedParticles->periodicity().toggle(1,lattice->periodicity().get(1));
  immersedParticles->periodicity().toggle(2,lattice->periodicity().get(2));

  immersedParticles->toggleInternalStatistics(false);
  
  if (del_sbStruct) {
    delete sbStructure;
  }
  
  InitAfterLoadCheckpoint();
}

void HemoCellFields::createCEPACfield() {
  SparseBlockStructure3D* sbStructure = lattice->getSparseBlockStructure().clone();
  ThreadAttribution * tAttribution = lattice->getMultiBlockManagement().getThreadAttribution().clone();
  plint refinement = lattice->getMultiBlockManagement().getRefinementLevel();
  lattice->getBlockCommunicator();
  CEPACfield = new MultiBlockLattice3D<T,CEPAC_DESCRIPTOR>(
          MultiBlockManagement3D( *sbStructure,
                                  tAttribution,
                                  envelopeSize,
                                  refinement ),
          plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
          plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
          plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T,CEPAC_DESCRIPTOR>(),
          new plb::AdvectionDiffusionBGKdynamics<T,CEPAC_DESCRIPTOR>(param::tau_CEPAC)
          );
  
  CEPACfield->periodicity().toggle(0,lattice->periodicity().get(0));
  CEPACfield->periodicity().toggle(1,lattice->periodicity().get(1));
  CEPACfield->periodicity().toggle(2,lattice->periodicity().get(2));

  CEPACfield->toggleInternalStatistics(false);
  
  integrateProcessingFunctional ( // instead of integrateProcessingFunctional
    new LatticeToPassiveAdvDiff3D<T,DESCRIPTOR,CEPAC_DESCRIPTOR>(1),
    lattice->getBoundingBox(), *lattice, *CEPACfield, 1);

}

HemoCellField * HemoCellFields::addCellType(std::string name_, int constructType)
{
  HemoCellField * cf = new HemoCellField(*this, name_, cellFields.size(), constructType);
  cellFields.push_back(cf);
  return cf;
}

unsigned int HemoCellFields::size() 
{
  return this->cellFields.size();
}

HemoCellField * HemoCellFields::operator[](unsigned int index)
{
  if (index >= cellFields.size()) {
    hlog << "(HemoCellFields) Error, cellindex " << index << " requested, but there are only " << cellFields.size() << " celltypes." << endl;
    exit(1);
  }
  return cellFields[index];
}

HemoCellField * HemoCellFields::operator[](string name)
{
  for (unsigned int i = 0; i < cellFields.size(); i++) {
      if (cellFields[i]->name == name) {
          return cellFields[i];
      }
  } 
  hlog << "(Error) (CellFields) " << name << "Celltype requested but it does not exist" << endl;
	exit(1);	
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
    if (global.enableCEPACfield && CEPACfield) {
      immersedParticles->getComponent(blocks[iBlock]).CEPAClattice = &CEPACfield->getComponent(blocks[iBlock]);
    }
    immersedParticles->getComponent(blocks[iBlock]).envelopeSize = envelopeSize;
    
    BlockLattice3D<T,DESCRIPTOR> * fluid = immersedParticles->getComponent(blocks[iBlock]).atomicLattice;
    immersedParticles->getComponent(blocks[iBlock]).nFluidCells = 0;
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
  
    if (immersedParticles->getComponent(blocks[iBlock]).neighbours.size() > max_neighbours) {
      max_neighbours = immersedParticles->getComponent(blocks[iBlock]).neighbours.size();
    }
  }
}

void HemoCellFields::load(XMLreader * documentXML, unsigned int & iter, Config * cfg)
{

    std::string firstField = (*(documentXML->getChildren( documentXML->getFirstId() )[0])).getName();
    bool isCheckpointed = (firstField=="Checkpoint");
    if (isCheckpointed) {
      (*documentXML)["Checkpoint"]["General"]["Iteration"].read(iter);
      std::string outDir;
      (*documentXML)["Checkpoint"]["General"]["OutDirectory"].read(outDir);
      plb::global::directories().setOutputDir(outDir);
      loadDirectories(cfg,false);

      std::string & chkDir = hemo::global.checkpointDirectory;

      if (hemocell.preInlet) {
        plb::parallelIO::load(chkDir + "PRE_lattice", *hemocell.preinlet_lattice, true);
        plb::parallelIO::load(chkDir + "PRE_particleField", *preinlet_immersedParticles, true);
      }
      plb::parallelIO::load(chkDir + "lattice", *hemocell.domain_lattice, true);
      plb::parallelIO::load(chkDir + "particleField", *domain_immersedParticles, true);
    } else {
      pcout << "(HemoCell) (CellFields) loading checkpoint from non-checkpoint Config" << endl;
      std::string & chkDir = hemo::global.checkpointDirectory;
      if (hemocell.preInlet) {
        plb::parallelIO::load(chkDir + "PRE_lattice", *hemocell.preinlet_lattice, true);
        plb::parallelIO::load(chkDir + "PRE_particleField", *preinlet_immersedParticles, true);
      }
      plb::parallelIO::load(chkDir + "lattice", *hemocell.domain_lattice, true);
      plb::parallelIO::load(chkDir + "particleField", *domain_immersedParticles, true);
    
    }
    
    InitAfterLoadCheckpoint();
    syncEnvelopes();
    deleteIncompleteCells();
}

void HemoCellFields::save(XMLreader *xmlr, unsigned int iter, Config * cfg)
{
    XMLwriter xmlw;
    std::string firstField = (*(xmlr->getChildren( xmlr->getFirstId() )[0])).getName(); 
    bool isCheckpointed = (firstField=="Checkpoint");
    if (!isCheckpointed) { copyXMLreader2XMLwriter((*xmlr)["hemocell"], xmlw["Checkpoint"]); }
    else { copyXMLreader2XMLwriter((*xmlr)["Checkpoint"]["hemocell"], xmlw["Checkpoint"]); }

    std::string & outDir = hemo::global.checkpointDirectory;

    mkpath(outDir.c_str(), 0777);

    
    /* Rename files, for safety reasons */
    if (global::mpi().isMainProcessor()) {
        renameFileToDotOld(outDir + "lattice.dat");
        renameFileToDotOld(outDir + "lattice.plb");
        renameFileToDotOld(outDir + "particleField.dat");
        renameFileToDotOld(outDir + "particleField.plb");
        renameFileToDotOld(outDir + "checkpoint.xml");
        if (hemocell.preInlet) {
          renameFileToDotOld(outDir + "PRE_lattice.dat");
          renameFileToDotOld(outDir + "PRE_lattice.plb");
          renameFileToDotOld(outDir + "PRE_particleField.dat");
          renameFileToDotOld(outDir + "PRE_particleField.plb");
        }
    } 
    
    global::mpi().barrier();
    
    /* Save XML & Data */
    xmlw["Checkpoint"]["General"]["Iteration"].set(iter);
    xmlw["Checkpoint"]["General"]["OutDirectory"].set(plb::global::directories().getOutputDir());
    xmlw.print(outDir + "checkpoint.xml");

    if (hemocell.preInlet) {
      plb::parallelIO::save(*hemocell.preinlet_lattice, outDir + "PRE_lattice", true);
      plb::parallelIO::save(*preinlet_immersedParticles, outDir + "PRE_particleField", true);
    }

    plb::parallelIO::save(*hemocell.domain_lattice, outDir + "lattice", true);
    plb::parallelIO::save(*domain_immersedParticles, outDir + "particleField", true);
}

void readPositionsCellFields(std::string particlePosFile) {
}


//void HemoCellFields::

void HemoCellFields::HemoFindInternalParticleGridPoints::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->findInternalParticleGridPoints(domain);
}

void HemoCellFields::findInternalParticleGridPoints() {
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(immersedParticles);
    applyProcessingFunctional(new HemoFindInternalParticleGridPoints(),immersedParticles->getBoundingBox(),wrapper);
}


void HemoCellFields::HemoInternalGridPointsMembrane::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->internalGridPointsMembrane(domain);
}

void HemoCellFields::internalGridPointsMembrane() {
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(immersedParticles);
    applyProcessingFunctional(new HemoInternalGridPointsMembrane(),immersedParticles->getBoundingBox(),wrapper);

}


void HemoCellFields::HemoInterpolateFluidVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->interpolateFluidVelocity(domain);
}
void HemoCellFields::interpolateFluidVelocity() {
  global.statistics.getCurrent()["interpolateFluidVelocity"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoInterpolateFluidVelocity(),immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::calculateCommunicationStructure() {
  MultiBlockManagement3D management_temp(immersedParticles->getMultiBlockManagement());
  ParallelBlockCommunicator3D * communicator = dynamic_cast<ParallelBlockCommunicator3D const *>(&immersedParticles->getBlockCommunicator())->clone();
  communicator->duplicateOverlaps(management_temp,immersedParticles->periodicity());
  large_communicator = new CommunicationStructure3D(*communicator->communication);
  immersedParticles->getMultiBlockManagement().changeEnvelopeWidth(3);
  immersedParticles->signalPeriodicity();
  immersedParticles->getBlockCommunicator().duplicateOverlaps(*immersedParticles,modif::hemocell_no_comm);
  delete communicator;
}

void HemoCellFields::HemoSyncEnvelopes::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->syncEnvelopes();
}
void HemoCellFields::syncEnvelopes() {
  global.statistics.getCurrent()["syncEnvelopes"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks() ) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    pf.removeParticles_inverse(pf.localDomain);
  }
  immersedParticles->getBlockCommunicator().duplicateOverlaps(*immersedParticles,modif::hemocell);
  
  if (large_communicator) {
  
    CommunicationStructure3D * comms = large_communicator;
    std::set<int> recv_procs, send_procs;
    std::map<int,vector<CommunicationInfo3D const *>> recv_infos, send_infos;
    for (CommunicationInfo3D const& info : comms->sendPackage) {
      send_procs.insert(info.toProcessId);
      send_infos[info.toProcessId].push_back(&info);
    }
    for (CommunicationInfo3D const& info : comms->recvPackage) {
      recv_procs.insert(info.fromProcessId);
      recv_infos[info.fromProcessId].push_back(&info);
    }

    set<int> locals;
    for (plint lbid : immersedParticles->getLocalInfo().getBlocks() ) {
      HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
      for(HemoCellParticle & particle: pf.particles) {
        locals.insert(particle.sv.cellId);
      }
    }
    vector<int> locals_v;
    locals_v.insert(locals_v.end(),locals.begin(),locals.end());
         
    vector<int> recv_procs_v;
    recv_procs_v.insert(recv_procs_v.end(),recv_procs.begin(),recv_procs.end());
    vector<MPI_Request> reqs(recv_procs.size());

    for (unsigned int i = 0 ; i < recv_procs_v.size() ; i ++) {
      MPI_Isend(&locals_v[0],locals_v.size(),MPI_INT,recv_procs_v[i],24,MPI_COMM_WORLD,&reqs[i]);
    }
    sendBuffers.resize(send_procs.size());
    for (unsigned int i = 0 ; i < send_procs.size() ; i ++) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE,24,MPI_COMM_WORLD,&status);
      int count;
      MPI_Get_count(&status,MPI_INT,&count);
      vector<int> requested_ids(count);
      vector<NoInitChar> & sendBuffer = sendBuffers[i];
      sendBuffer.clear();
      int offset = 0 ;
      MPI_Recv(&requested_ids[0],count,MPI_INT,status.MPI_SOURCE,24,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      for (CommunicationInfo3D const * info : send_infos[status.MPI_SOURCE] ) {
        HemoCellParticleField & pf = immersedParticles->getComponent(info->fromBlockId);
        int offset_p = pf.getDataTransfer().getOffset(info->absoluteOffset);
        const map<int,vector<int>> & ppc = pf.get_particles_per_cell();
        
        for (int id : requested_ids) {
          if (((offset_p < 0) && (id > INT_MAX+offset_p)) ||
              ((offset_p > 0) && (id < INT_MIN+offset_p))) {
            cout << "(HemoCellFields syncEnvelopes) Almost invoking overflow in periodic particle communication, resetting ID to base ID instead, this will most likely delete the particle" << endl;
            id = base_cell_id(id);
          } else {
            id = id - offset_p;
          }
          if (ppc.find(id) == ppc.end()) { continue; }
          for (int pid : ppc.at(id)) {
            if (pid <= -1) { continue; }
            if (pid >= (int) pf.particles.size()) { continue; }
            sendBuffer.resize(sendBuffer.size()+sizeof(HemoCellParticle::serializeValues_t));
            *((HemoCellParticle::serializeValues_t*)&sendBuffer[offset]) = pf.particles[pid].sv;
            offset += sizeof(HemoCellParticle::serializeValues_t);
          }         
        }
      }
      reqs.emplace_back();
      MPI_Isend(sendBuffer.data(),sendBuffer.size(),MPI_CHAR,status.MPI_SOURCE,42,MPI_COMM_WORLD,&reqs.back());
    }

    vector<MPI_Request> recv_reqs(recv_procs.size());
    recvBuffers.resize(recv_procs.size());
    for (unsigned int i = 0 ; i < recv_procs.size() ; i ++) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE,42,MPI_COMM_WORLD,&status);
      int count;
      MPI_Get_count(&status,MPI_CHAR,&count);
      vector<NoInitChar> & recv_buffer = recvBuffers[i];
      recv_buffer.resize(count);
      MPI_Irecv(recv_buffer.data(),count,MPI_CHAR,status.MPI_SOURCE,42,MPI_COMM_WORLD,&recv_reqs[i]);
    }

    for (unsigned int i = 0 ; i < recv_procs.size() ; i ++) {
      int index;
      MPI_Status status;
      if (MPI_SUCCESS != MPI_Waitany(recv_reqs.size(),&recv_reqs[0],&index,&status)) {
        hlog << "(HemoCellFields) (syncenvelopes) error returned in WaitAny " << status.MPI_ERROR << endl;
        exit(1);
      }
      recv_reqs[index] = MPI_REQUEST_NULL;
      //Get Offsets and Destinations
      for (CommunicationInfo3D const * info : recv_infos[status.MPI_SOURCE]) {
        HEMOCELL_PARTICLE_FIELD& toBlock = immersedParticles->getComponent(info->toBlockId);
        toBlock.getDataTransfer().receive (info->toDomain, recvBuffers[index], info->absoluteOffset );
      }
    }

    MPI_Waitall(reqs.size(),reqs.data(),MPI_STATUSES_IGNORE);
    
    // 3. Local copies which require no communication.
    for (unsigned iSendRecv=0; iSendRecv<comms->sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo3D const& info = comms->sendRecvPackage[iSendRecv];
        AtomicBlock3D const& fromBlock = immersedParticles->getComponent(info.fromBlockId);
        AtomicBlock3D& toBlock = immersedParticles->getComponent(info.toBlockId);

        toBlock.getDataTransfer().attribute (
                info.toDomain, 0, 0, 0 , fromBlock,
                modif::hemocell, info.absoluteOffset );
    }
    
  }
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoAdvanceParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->advanceParticles();
}
void HemoCellFields::advanceParticles() {
  global.statistics.getCurrent()["advanceParticles"].start();
    
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoAdvanceParticles(),immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSpreadParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->spreadParticleForce(domain);
}
void HemoCellFields::spreadParticleForce() {
  global.statistics.getCurrent()["spreadParticleForce"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSpreadParticleForce(),immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoApplyConstitutiveModel::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyConstitutiveModel(forced);
}
void HemoCellFields::applyConstitutiveModel(bool forced) {
  global.statistics.getCurrent()["applyConstitutiveModel"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  HemoApplyConstitutiveModel * fnct = new HemoApplyConstitutiveModel();
  fnct->forced = forced;
  applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoUnifyForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->unifyForceVectors();
}
void HemoCellFields::unify_force_vectors() {
  global.statistics.getCurrent()["unifyForceVectors"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoUnifyForceVectors(),immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyRepulsionForce();
}
void HemoCellFields::applyRepulsionForce() {
  global.statistics.getCurrent()["repulsionForce"].start();

  vector<MultiBlock3D*>wrapper;
  wrapper.push_back(immersedParticles);
  HemoRepulsionForce * fnct = new HemoRepulsionForce();
  applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoBoundaryRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->applyBoundaryRepulsionForce();
}
void HemoCellFields::applyBoundaryRepulsionForce() {
  global.statistics.getCurrent()["boundaryRepulsionForce"].start();

  vector<MultiBlock3D*>wrapper;
  wrapper.push_back(immersedParticles);
  HemoBoundaryRepulsionForce * fnct = new HemoBoundaryRepulsionForce();
  applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoPopulateBoundaryParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->populateBoundaryParticles();
}
void HemoCellFields::populateBoundaryParticles() {
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    HemoPopulateBoundaryParticles * fnct = new HemoPopulateBoundaryParticles();
    applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoPopulateBindingSites::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->populateBindingSites(domain);
}
void HemoCellFields::populateBindingSites(plb::Box3D * box) {
  //Initialize bindingField before entering functional
  bindingFieldHelper::get(*this);
  
  vector<MultiBlock3D*>wrapper;
  wrapper.push_back(immersedParticles);
  Box3D domain;
  if(box) {
    domain = *box;
  } else {
    domain = immersedParticles->getBoundingBox();
  }
  HemoPopulateBindingSites * fnct = new HemoPopulateBindingSites();
  applyProcessingFunctional(fnct,domain,wrapper);
}
void HemoCellFields::HemoupdateResidenceTime::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->updateResidenceTime(rtime);
}
void HemoCellFields::updateResidenceTime(unsigned int rtime) {
    global.statistics.getCurrent()["updateResidenceTime"].start();
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    HemoupdateResidenceTime * fnct = new HemoupdateResidenceTime();
    fnct->rtime = rtime;
    applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
    global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSeperateForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->separateForceVectors();
}
void HemoCellFields::separate_force_vectors() {
  global.statistics.getCurrent()["separateForceVectors"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSeperateForceVectors(),immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}
void HemoCellFields::HemoDeleteIncompleteCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->deleteIncompleteCells(verbose);
}
void HemoCellFields::deleteIncompleteCells(bool verbose) {
  global.statistics.getCurrent()["deleteIncompleteCells"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  HemoDeleteIncompleteCells * fnct = new HemoDeleteIncompleteCells();
  fnct->verbose = verbose;
  applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
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
    pf->addParticle(&particle);
  }
}
void HemoCellFields::addParticles(vector<HemoCellParticle> & particles) {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSetParticles(particles),immersedParticles->getBoundingBox(),wrapper);
}


void HemoCellFields::HemoDeleteNonLocalParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD * pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  pf->removeParticles_inverse(pf->localDomain.enlarge(envelopeSize));
}
void HemoCellFields::deleteNonLocalParticles(int envelope) {
  global.statistics.getCurrent()["deleteNonLocalParticles"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoDeleteNonLocalParticles(envelope),immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSolidifyCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD * pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  pf->solidifyCells();
}
void HemoCellFields::solidifyCells() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSolidifyCells(),immersedParticles->getBoundingBox(),wrapper);
}


HemoCellFields::HemoInternalGridPointsMembrane *  HemoCellFields::HemoInternalGridPointsMembrane::clone() const { return new HemoCellFields::HemoInternalGridPointsMembrane(*this);}
HemoCellFields::HemoFindInternalParticleGridPoints *  HemoCellFields::HemoFindInternalParticleGridPoints::clone() const { return new HemoCellFields::HemoFindInternalParticleGridPoints(*this);}
HemoCellFields::HemoSeperateForceVectors * HemoCellFields::HemoSeperateForceVectors::clone() const { return new HemoCellFields::HemoSeperateForceVectors(*this);}
HemoCellFields::HemoUnifyForceVectors *    HemoCellFields::HemoUnifyForceVectors::clone() const    { return new HemoCellFields::HemoUnifyForceVectors(*this);}
HemoCellFields::HemoSpreadParticleForce *  HemoCellFields::HemoSpreadParticleForce::clone() const { return new HemoCellFields::HemoSpreadParticleForce(*this);}
HemoCellFields::HemoInterpolateFluidVelocity * HemoCellFields::HemoInterpolateFluidVelocity::clone() const { return new HemoCellFields::HemoInterpolateFluidVelocity(*this);}
HemoCellFields::HemoAdvanceParticles *     HemoCellFields::HemoAdvanceParticles::clone() const { return new HemoCellFields::HemoAdvanceParticles(*this);}
HemoCellFields::HemoApplyConstitutiveModel * HemoCellFields::HemoApplyConstitutiveModel::clone() const { return new HemoCellFields::HemoApplyConstitutiveModel(*this);}
HemoCellFields::HemoSyncEnvelopes *        HemoCellFields::HemoSyncEnvelopes::clone() const { return new HemoCellFields::HemoSyncEnvelopes(*this);}
HemoCellFields::HemoRepulsionForce *        HemoCellFields::HemoRepulsionForce::clone() const { return new HemoCellFields::HemoRepulsionForce(*this);}
HemoCellFields::HemoBoundaryRepulsionForce *        HemoCellFields::HemoBoundaryRepulsionForce::clone() const { return new HemoCellFields::HemoBoundaryRepulsionForce(*this);}
HemoCellFields::HemoDeleteIncompleteCells *        HemoCellFields::HemoDeleteIncompleteCells::clone() const { return new HemoCellFields::HemoDeleteIncompleteCells(*this);}
HemoCellFields::HemoGetParticles *        HemoCellFields::HemoGetParticles::clone() const { return new HemoCellFields::HemoGetParticles(*this);}
HemoCellFields::HemoSetParticles *        HemoCellFields::HemoSetParticles::clone() const { return new HemoCellFields::HemoSetParticles(*this);}
HemoCellFields::HemoPopulateBoundaryParticles *        HemoCellFields::HemoPopulateBoundaryParticles::clone() const { return new HemoCellFields::HemoPopulateBoundaryParticles(*this);}
HemoCellFields::HemoDeleteNonLocalParticles *        HemoCellFields::HemoDeleteNonLocalParticles::clone() const { return new HemoCellFields::HemoDeleteNonLocalParticles(*this);}
HemoCellFields::HemoSolidifyCells *        HemoCellFields::HemoSolidifyCells::clone() const { return new HemoCellFields::HemoSolidifyCells(*this);}
HemoCellFields::HemoPopulateBindingSites * HemoCellFields::HemoPopulateBindingSites::clone() const { return new HemoCellFields::HemoPopulateBindingSites(*this);}
HemoCellFields::HemoupdateResidenceTime * HemoCellFields::HemoupdateResidenceTime::clone() const { return new HemoCellFields::HemoupdateResidenceTime(*this);}


void HemoCellFields::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
   for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::hemocell;
   }
}

MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & HemoCellFields::getParticleField3D() { return *immersedParticles; };

}
