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
#include <algorithm>

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
    preinlet_immersedParticles = new MultiParticleField3D<HemoCellParticleField>(MultiBlockManagement3D(
      *hemocell.preinlet_lattice->getSparseBlockStructure().clone(),
      hemocell.preinlet_lattice->getMultiBlockManagement().getThreadAttribution().clone(),
      envelopeSize,
      refinement ), plb::defaultMultiBlockPolicy3D().getCombinedStatistics() );
  }
  domain_immersedParticles = new MultiParticleField3D<HemoCellParticleField>(MultiBlockManagement3D(
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->findInternalParticleGridPoints(domain);
}

void HemoCellFields::findInternalParticleGridPoints() {
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(immersedParticles);
    applyProcessingFunctional(new HemoFindInternalParticleGridPoints(),immersedParticles->getBoundingBox(),wrapper);
}


void HemoCellFields::HemoInternalGridPointsMembrane::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->internalGridPointsMembrane(domain);
}

void HemoCellFields::internalGridPointsMembrane() {
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(immersedParticles);
    applyProcessingFunctional(new HemoInternalGridPointsMembrane(),immersedParticles->getBoundingBox(),wrapper);

}


void HemoCellFields::HemoInterpolateFluidVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->interpolateFluidVelocity(domain);
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->syncEnvelopes();
}
void HemoCellFields::syncEnvelopes() {
  global.statistics.getCurrent()["syncEnvelopes"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks() ) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    // enlarge(3) keeps Palabos 3-LU envelope copies alive through removeParticles_inverse.
    // Without this, the dead-zone: vertex at x∈(202.0,202.5] maps to B0 copy at x∈(-1,-0.5],
    // which localDomain threshold (> -0.5) removes. Fix 2 then finds the copy via cid±N lookup.
    pf.removeParticles_inverse(pf.localDomain.enlarge(3));
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

    // Fixed pull protocol: targeted, race-free particle sync.
    // Phase 1: exchange needs (incomplete cellIds) with neighbors using specific-source
    //          MPI_Probe — eliminates the ANY_SOURCE race of the original pull protocol.
    // Phase 2: each rank sends only the particles that were explicitly requested.

    vector<int> send_procs_v(send_procs.begin(), send_procs.end());
    vector<int> recv_procs_v(recv_procs.begin(), recv_procs.end());

    // Phase 1a: scan local ppc for incomplete cells (after the 3-LU Palabos sync above).
    set<int> local_incomplete;
    map<plint, set<int>> block_incomplete;  // lbid -> set of incomplete cids on that block
    for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
      const map<int,vector<int>> & ppc = pf.get_particles_per_cell();
      for (const auto & kv : ppc) {
        for (int idx : kv.second) {
          if (idx == -1) {
            local_incomplete.insert(kv.first);
            block_incomplete[lbid].insert(kv.first);
            break;
          }
        }
      }
    }
    // base_to_needed_cids: maps base_cell_id(cid) -> set of cids that need supply.
    // Enables matching periodic cellIds: sender packs sv.cellId=base_cid; receiver
    // remaps to the requested periodic cid before insertion (no absoluteOffset applied).
    map<int, set<int>> base_to_needed_cids;
    for (int cid : local_incomplete)
      base_to_needed_cids[base_cell_id(cid)].insert(cid);
    // Domain-wide BB expansion for targeted delivery — accepts any stretch or periodic offset.
    Box3D domain_bb = immersedParticles->getBoundingBox();
    plint expand = std::max({domain_bb.x1 - domain_bb.x0 + 1,
                             domain_bb.y1 - domain_bb.y0 + 1,
                             domain_bb.z1 - domain_bb.z0 + 1});
    // Domain lengths and period_crosses table used in Fix 2 (cid±shift lookup).
    // getOffset() encodes crossing direction as: ±1*N for x, ±limit_y*N for y, ±limit_z*N for z.
    // For each k_shift, the position correction undoes the absoluteOffset Palabos applied.
    // k=+1 (x>0 crossing): absoluteOffset.x=+domain_x was added → correction=-domain_x, etc.
    plint domain_x = domain_bb.x1 - domain_bb.x0 + 1;
    plint domain_y = domain_bb.y1 - domain_bb.y0 + 1;
    plint domain_z = domain_bb.z1 - domain_bb.z0 + 1;
    struct PeriodCross { int k; T dx, dy, dz; };
    vector<PeriodCross> period_crosses;
    period_crosses.push_back({ +1, -(T)domain_x,       0,           0       });
    period_crosses.push_back({ -1, +(T)domain_x,       0,           0       });
    if (periodicity_limit_offset_y != 0) {
      period_crosses.push_back({ +periodicity_limit_offset_y,  0, -(T)domain_y,  0 });
      period_crosses.push_back({ -periodicity_limit_offset_y,  0, +(T)domain_y,  0 });
    }
    if (periodicity_limit_offset_z != 0) {
      period_crosses.push_back({ +periodicity_limit_offset_z,  0,  0, -(T)domain_z });
      period_crosses.push_back({ -periodicity_limit_offset_z,  0,  0, +(T)domain_z });
    }
    // Combined-axis crossings: a vertex that simultaneously crosses two or three
    // periodic boundaries in one step gets a k_net that is the SUM of single-axis
    // k values.  Palabos applies the combined absoluteOffset to the position and
    // getOffset() encodes both axes into one integer, so Fix 2 must look for
    // cid ± k_combined×N.  Position correction undoes BOTH shifts at once.
    // Example (xy-diagonal, k=-101 = k_x=-1+k_y=-100): vertex at global (103,103)
    // in block 31 is pruned by removeParticles_inverse, sent by Palabos to block 0
    // as cid=-122463 at position (3,3). Without k=-101 in the table Fix 2 cannot
    // recover it, and deleteIncompleteCells fires.  With it, block 0's ppc[-122463]
    // is found, position corrected to (103,103), and delivered to block 31.
    if (periodicity_limit_offset_y != 0) {
      // xy combined (4 sign combinations)
      const int ky = periodicity_limit_offset_y;
      period_crosses.push_back({ +1+ky, -(T)domain_x, -(T)domain_y,  0 });
      period_crosses.push_back({ -1+ky, +(T)domain_x, -(T)domain_y,  0 });
      period_crosses.push_back({ +1-ky, -(T)domain_x, +(T)domain_y,  0 });
      period_crosses.push_back({ -1-ky, +(T)domain_x, +(T)domain_y,  0 });
    }
    if (periodicity_limit_offset_z != 0) {
      // xz combined (4 sign combinations)
      const int kz = periodicity_limit_offset_z;
      period_crosses.push_back({ +1+kz, -(T)domain_x,  0, -(T)domain_z });
      period_crosses.push_back({ -1+kz, +(T)domain_x,  0, -(T)domain_z });
      period_crosses.push_back({ +1-kz, -(T)domain_x,  0, +(T)domain_z });
      period_crosses.push_back({ -1-kz, +(T)domain_x,  0, +(T)domain_z });
    }
    if (periodicity_limit_offset_y != 0 && periodicity_limit_offset_z != 0) {
      const int ky = periodicity_limit_offset_y;
      const int kz = periodicity_limit_offset_z;
      // yz combined (4 sign combinations)
      period_crosses.push_back({ +ky+kz,  0, -(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ -ky+kz,  0, +(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ +ky-kz,  0, -(T)domain_y, +(T)domain_z });
      period_crosses.push_back({ -ky-kz,  0, +(T)domain_y, +(T)domain_z });
      // xyz combined (8 sign combinations)
      period_crosses.push_back({ +1+ky+kz, -(T)domain_x, -(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ -1+ky+kz, +(T)domain_x, -(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ +1-ky+kz, -(T)domain_x, +(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ -1-ky+kz, +(T)domain_x, +(T)domain_y, -(T)domain_z });
      period_crosses.push_back({ +1+ky-kz, -(T)domain_x, -(T)domain_y, +(T)domain_z });
      period_crosses.push_back({ -1+ky-kz, +(T)domain_x, -(T)domain_y, +(T)domain_z });
      period_crosses.push_back({ +1-ky-kz, -(T)domain_x, +(T)domain_y, +(T)domain_z });
      period_crosses.push_back({ -1-ky-kz, +(T)domain_x, +(T)domain_y, +(T)domain_z });
    }
    vector<int> needs_v(local_incomplete.begin(), local_incomplete.end());

    // Phase 1b: send my needs to suppliers (recv_procs); receive needs from requesters
    // (send_procs). MPI_Probe uses specific source on both sides — no ANY_SOURCE, no race.
    vector<MPI_Request> needs_reqs(recv_procs_v.size());
    for (unsigned int i = 0; i < recv_procs_v.size(); i++) {
      MPI_Isend(needs_v.data(), (int)needs_v.size(), MPI_INT,
                recv_procs_v[i], 24, MPI_COMM_WORLD, &needs_reqs[i]);
    }
    map<int, vector<int>> requested_by;
    for (int r : send_procs_v) {
      MPI_Status status;
      MPI_Probe(r, 24, MPI_COMM_WORLD, &status);  // specific source: safe, no race
      int cnt;
      MPI_Get_count(&status, MPI_INT, &cnt);
      requested_by[r].resize(cnt);
      MPI_Recv(requested_by[r].data(), cnt, MPI_INT,
               r, 24, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Waitall((int)needs_reqs.size(), needs_reqs.data(), MPI_STATUSES_IGNORE);

    // Phase 2a: build targeted supply buffers and send to requesters (send_procs).
    // Periodic fallback: cid outside [0, number_of_cells) means a periodic re-entry;
    // look up base_cell_id and remap sv.cellId so the receiver's ppc inserts under the
    // correct (offset) key.
    sendBuffers.resize(send_procs_v.size());
    vector<MPI_Request> supply_reqs(send_procs_v.size());
    for (unsigned int i = 0; i < send_procs_v.size(); i++) {
      int r = send_procs_v[i];
      vector<NoInitChar> & buf = sendBuffers[i];
      buf.clear();
      for (int cid : requested_by[r]) {
        for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
          HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
          const map<int,vector<int>> & ppc = pf.get_particles_per_cell();
          auto it = ppc.find(cid);
          if (it != ppc.end()) {
            for (int pid : it->second) {
              if (pid < 0 || pid >= (int)pf.particles.size()) continue;
              size_t off = buf.size();
              buf.resize(off + sizeof(HemoCellParticle::serializeValues_t));
              *((HemoCellParticle::serializeValues_t*)&buf[off]) = pf.particles[pid].sv;
            }
          }
          if (cid < 0 || cid >= number_of_cells) {
            int base_cid = base_cell_id(cid);
            auto it2 = ppc.find(base_cid);
            if (it2 != ppc.end()) {
              for (int pid : it2->second) {
                if (pid < 0 || pid >= (int)pf.particles.size()) continue;
                // Send with original sv.cellId (= base_cid) and unmodified position.
                // Phase 2b targeted delivery remaps sv.cellId to the requested periodic
                // cid on the receiver side using base_to_needed_cids, without applying
                // any absoluteOffset to position. Full-domain BB expansion accepts the
                // absolute position regardless of stretch or periodic configuration.
                size_t off = buf.size();
                buf.resize(off + sizeof(HemoCellParticle::serializeValues_t));
                *((HemoCellParticle::serializeValues_t*)&buf[off]) = pf.particles[pid].sv;
              }
            }
          }
          // Fix 2: periodic-offset copy lookup (cid ± shift) kept alive by enlarge(3).
          // A neighbor block may hold the vertex under cid±shift*N after Palabos applied
          // absoluteOffset to position and ±shift*N to cellId. Remap cellId to the
          // requested cid and undo the position offset (per period_crosses table).
          for (const PeriodCross & pc : period_crosses) {
            int cid_shifted = cid + pc.k * number_of_cells;
            auto it_sh = ppc.find(cid_shifted);
            if (it_sh != ppc.end()) {
              for (int pid : it_sh->second) {
                if (pid < 0 || pid >= (int)pf.particles.size()) continue;
                HemoCellParticle::serializeValues_t sv2 = pf.particles[pid].sv;
                sv2.cellId = cid;
                sv2.position[0] += pc.dx;
                sv2.position[1] += pc.dy;
                sv2.position[2] += pc.dz;
                size_t off = buf.size();
                buf.resize(off + sizeof(HemoCellParticle::serializeValues_t));
                *((HemoCellParticle::serializeValues_t*)&buf[off]) = sv2;
              }
            }
          }
        }
      }
      MPI_Isend(buf.data(), (int)buf.size(), MPI_CHAR,
                r, 42, MPI_COMM_WORLD, &supply_reqs[i]);
    }

    // Phase 2b: receive supply from suppliers (recv_procs) and deliver to local blocks.
    // Targeted delivery: parse buffer directly, deliver only to blocks that need each
    // cellId. sv.cellId from sender is base_cid (never pre-offset by Phase 2a). For
    // periodic cids, remap using base_to_needed_cids. Position stays at original absolute
    // coords — full-domain BB expansion accepts any stretch without position shifting.
    // MPI_Probe uses specific source — no ANY_SOURCE, no race.
    recvBuffers.resize(recv_procs_v.size());
    for (unsigned int i = 0; i < recv_procs_v.size(); i++) {
      int r = recv_procs_v[i];
      MPI_Status status;
      MPI_Probe(r, 42, MPI_COMM_WORLD, &status);  // specific source: safe, no race
      int cnt;
      MPI_Get_count(&status, MPI_CHAR, &cnt);
      recvBuffers[i].resize(cnt);
      MPI_Recv(recvBuffers[i].data(), cnt, MPI_CHAR,
               r, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int n_sv = cnt / (int)sizeof(HemoCellParticle::serializeValues_t);
      for (int pi = 0; pi < n_sv; pi++) {
        HemoCellParticle::serializeValues_t sv;
        memcpy(&sv, &recvBuffers[i][pi * sizeof(HemoCellParticle::serializeValues_t)], sizeof(sv));
        // Non-periodic: sv.cellId is directly the needed cid.
        if (local_incomplete.count(sv.cellId)) {
          for (auto & blk_kv : block_incomplete) {
            if (!blk_kv.second.count(sv.cellId)) continue;
            HemoCellParticleField & pf = immersedParticles->getComponent(blk_kv.first);
            Box3D saved = pf.getBoundingBox();
            pf.getBoundingBox() = saved.enlarge(expand);
            pf.addParticle(sv);
            pf.getBoundingBox() = saved;
          }
        }
        // Periodic fallback: sv.cellId is base_cid; remap to each requesting periodic cid.
        auto bt = base_to_needed_cids.find(sv.cellId);
        if (bt != base_to_needed_cids.end()) {
          for (int target_cid : bt->second) {
            if (target_cid == sv.cellId) continue;  // non-periodic already handled above
            HemoCellParticle::serializeValues_t sv2 = sv;
            sv2.cellId = target_cid;
            for (auto & blk_kv : block_incomplete) {
              if (!blk_kv.second.count(target_cid)) continue;
              HemoCellParticleField & pf = immersedParticles->getComponent(blk_kv.first);
              Box3D saved = pf.getBoundingBox();
              pf.getBoundingBox() = saved.enlarge(expand);
              pf.addParticle(sv2);
              pf.getBoundingBox() = saved;
            }
          }
        }
      }
    }
    MPI_Waitall((int)supply_reqs.size(), supply_reqs.data(), MPI_STATUSES_IGNORE);
    
    // 3. Local copies which require no communication.
    for (unsigned iSendRecv=0; iSendRecv<comms->sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo3D const& info = comms->sendRecvPackage[iSendRecv];
        AtomicBlock3D const& fromBlock = immersedParticles->getComponent(info.fromBlockId);
        AtomicBlock3D& toBlock = immersedParticles->getComponent(info.toBlockId);

        toBlock.getDataTransfer().attribute (
                info.toDomain, 0, 0, 0 , fromBlock,
                modif::hemocell, info.absoluteOffset );
    }

    // Same-rank targeted supply: fill remaining incomplete cells from other local blocks.
    // Mirrors Phase 2b targeted delivery but for same-rank block pairs — handles stretched
    // cells and periodic crossings that attribute() misses due to its 3-LU BB limit.
    if (!block_incomplete.empty()) {
      for (auto & blk_kv : block_incomplete) {
        HemoCellParticleField & dst = immersedParticles->getComponent(blk_kv.first);
        for (plint src_lbid : immersedParticles->getLocalInfo().getBlocks()) {
          if (src_lbid == blk_kv.first) continue;
          HemoCellParticleField & src = immersedParticles->getComponent(src_lbid);
          const map<int,vector<int>> & src_ppc = src.get_particles_per_cell();
          for (int cid : blk_kv.second) {
            // Direct match: src has this cid (non-periodic or already-offset periodic copy)
            auto it = src_ppc.find(cid);
            if (it != src_ppc.end()) {
              for (int pid : it->second) {
                if (pid < 0 || pid >= (int)src.particles.size()) continue;
                Box3D saved = dst.getBoundingBox();
                dst.getBoundingBox() = saved.enlarge(expand);
                dst.addParticle(src.particles[pid].sv);
                dst.getBoundingBox() = saved;
              }
            }
            // Periodic fallback: cid is a periodic re-entry; src may have base_cid instead
            if (cid < 0 || cid >= number_of_cells) {
              int base = base_cell_id(cid);
              auto it2 = src_ppc.find(base);
              if (it2 != src_ppc.end()) {
                for (int pid : it2->second) {
                  if (pid < 0 || pid >= (int)src.particles.size()) continue;
                  HemoCellParticle::serializeValues_t sv2 = src.particles[pid].sv;
                  sv2.cellId = cid;  // remap base_cid to requested periodic cid
                  Box3D saved = dst.getBoundingBox();
                  dst.getBoundingBox() = saved.enlarge(expand);
                  dst.addParticle(sv2);
                  dst.getBoundingBox() = saved;
                }
              }
            }
            // Fix 2 (same-rank): periodic-offset copy lookup using period_crosses table.
            // Mirrors Phase 2a Fix 2: src block holds vertex under cid±shift*N (Palabos
            // envelope copy). Remap cellId and undo position offset before delivering.
            for (const PeriodCross & pc : period_crosses) {
              int cid_shifted = cid + pc.k * number_of_cells;
              auto it_sh = src_ppc.find(cid_shifted);
              if (it_sh != src_ppc.end()) {
                for (int pid : it_sh->second) {
                  if (pid < 0 || pid >= (int)src.particles.size()) continue;
                  HemoCellParticle::serializeValues_t sv2 = src.particles[pid].sv;
                  sv2.cellId = cid;
                  sv2.position[0] += pc.dx;
                  sv2.position[1] += pc.dy;
                  sv2.position[2] += pc.dz;
                  Box3D saved = dst.getBoundingBox();
                  dst.getBoundingBox() = saved.enlarge(expand);
                  dst.addParticle(sv2);
                  dst.getBoundingBox() = saved;
                }
              }
            }
          }
        }
      }
    }

  // Diagnostic: re-scan for cells still incomplete after ALL Phase 1/2 delivery.
  // If any remain, Phase 1/2 failed to supply them — either no neighbor rank has the
  // vertex in its ppc, or addParticle rejected it. TCR will attempt next; if TCR also
  // fails (no rank has the vertex anywhere), cascade deletion follows.
  {
    int still_count = 0;
    int me_rank = global::mpi().getRank();
    for (const auto& blk_kv : block_incomplete) {
      HemoCellParticleField & pf2 = immersedParticles->getComponent(blk_kv.first);
      const auto & ppc2 = pf2.get_particles_per_cell();
      for (int cid : blk_kv.second) {
        auto it2 = ppc2.find(cid);
        if (it2 == ppc2.end()) { still_count++; continue; }
        for (int idx : it2->second) {
          if (idx == -1) { still_count++; break; }
        }
      }
    }
    if (still_count > 0) {
      hlog << "(Phase1/2-diag) iter=" << hemocell.iter << " rank=" << me_rank
           << ": " << still_count << " cid(s) still incomplete after Phase 1/2 supply\n";
    }
  }

  }
  global.statistics.getCurrent().stop();
}

void HemoCellFields::syncEnvelopesTargetedRepair() {
  // Step 1: scan local blocks for incomplete cells (-1 in ppc)
  map<plint, set<int>> block_incomplete;
  set<int> local_incomplete;
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    const map<int,vector<int>> & ppc = pf.get_particles_per_cell();
    for (const auto & kv : ppc) {
      for (int idx : kv.second) {
        if (idx == -1) {
          block_incomplete[lbid].insert(kv.first);
          local_incomplete.insert(kv.first);
          break;
        }
      }
    }
  }

  // Fast-exit: one Allreduce to check if any rank has incomplete cells.
  // Use SUM so global_count = total incomplete cells across all ranks (useful for logging).
  int local_count = (int)local_incomplete.size();
  int global_count = 0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (global_count == 0) return;

  // Log: rank 0 only, first occurrence and every 10000 iters. Uses global count so
  // it fires even when rank 0 itself has no incomplete cells.
  if (global::mpi().isMainProcessor()) {
    static bool tcr_first = true;
    static unsigned int tcr_last_log = 0;
    if (tcr_first || hemocell.iter - tcr_last_log >= 10000) {
      tcr_first = false;
      tcr_last_log = hemocell.iter;
      hlog << "(HemoCell) (TCR) iter=" << hemocell.iter
           << ": " << global_count << " incomplete cell(s) across all ranks, firing repair\n";
    }
  }

  int nranks = global::mpi().getSize();
  int me = global::mpi().getRank();

  // Allgather 1: share each rank's incomplete cellId needs globally
  vector<int> local_inc_v(local_incomplete.begin(), local_incomplete.end());
  int my_need_count = (int)local_inc_v.size();
  vector<int> all_need_counts(nranks);
  MPI_Allgather(&my_need_count, 1, MPI_INT,
                all_need_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  vector<int> need_displs(nranks + 1, 0);
  for (int r = 0; r < nranks; r++) need_displs[r+1] = need_displs[r] + all_need_counts[r];
  vector<int> all_needs(need_displs[nranks]);
  MPI_Allgatherv(local_inc_v.data(), my_need_count, MPI_INT,
                 all_needs.data(), all_need_counts.data(), need_displs.data(),
                 MPI_INT, MPI_COMM_WORLD);

  // Build per-rank need sets and global incomplete set
  set<int> global_incomplete_set;
  map<int, set<int>> needed_by_rank;
  for (int r = 0; r < nranks; r++) {
    for (int i = need_displs[r]; i < need_displs[r+1]; i++) {
      needed_by_rank[r].insert(all_needs[i]);
      global_incomplete_set.insert(all_needs[i]);
    }
  }

  // Supply scan: collect serialized particles for globally-incomplete cellIds.
  // For periodic re-entries (cid outside [0, number_of_cells)), also check the
  // base (non-offset) cellId and remap sv.cellId before sending so the receiver
  // inserts the particle under the correct ppc key.
  map<int, vector<HemoCellParticle::serializeValues_t>> my_supply;
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    const map<int,vector<int>> & ppc = pf.get_particles_per_cell();
    for (int cid : global_incomplete_set) {
      // Direct lookup: handles non-periodic cids and already-offset periodic copies
      auto it = ppc.find(cid);
      if (it != ppc.end()) {
        for (int pid : it->second) {
          if (pid == -1 || pid >= (int)pf.particles.size()) continue;
          my_supply[cid].push_back(pf.particles[pid].sv);
        }
      }
      // Periodic fallback: cid outside [0, number_of_cells) means this is a
      // periodic re-entry. The original owning block stores the cell under
      // base_cell_id(cid). Find those particles and remap cellId so the
      // receiver's update_ppc inserts under ppc[cid].
      if (cid < 0 || cid >= number_of_cells) {
        int base_cid = base_cell_id(cid);
        auto it2 = ppc.find(base_cid);
        if (it2 != ppc.end()) {
          for (int pid : it2->second) {
            if (pid == -1 || pid >= (int)pf.particles.size()) continue;
            HemoCellParticle::serializeValues_t sv = pf.particles[pid].sv;
            sv.cellId = cid;
            my_supply[cid].push_back(sv);
          }
        }
      }
    }
  }

  vector<int> my_supply_cids;
  for (const auto & kv : my_supply) my_supply_cids.push_back(kv.first);

  // Allgather 2: share supply (which cellIds each rank can provide) globally
  int my_supply_count = (int)my_supply_cids.size();
  vector<int> all_supply_counts(nranks);
  MPI_Allgather(&my_supply_count, 1, MPI_INT,
                all_supply_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  vector<int> supply_displs(nranks + 1, 0);
  for (int r = 0; r < nranks; r++) supply_displs[r+1] = supply_displs[r] + all_supply_counts[r];
  vector<int> all_supply(supply_displs[nranks]);
  MPI_Allgatherv(my_supply_cids.data(), my_supply_count, MPI_INT,
                 all_supply.data(), all_supply_counts.data(), supply_displs.data(),
                 MPI_INT, MPI_COMM_WORLD);

  // Build supply_by_rank[r] = set of cellIds rank r can provide
  map<int, set<int>> supply_by_rank;
  for (int r = 0; r < nranks; r++) {
    for (int i = supply_displs[r]; i < supply_displs[r+1]; i++) {
      supply_by_rank[r].insert(all_supply[i]);
    }
  }

  // Diagnostic: identify globally-incomplete cids with NO supplier on any rank.
  // If any exist, those vertices were deleted from all blocks (cascade deletion scenario).
  // Phase 1/2 and TCR cannot repair them — deleteIncompleteCells will follow.
  if (global::mpi().isMainProcessor()) {
    for (int cid : global_incomplete_set) {
      bool any_supplier = false;
      for (int r = 0; r < nranks; r++) {
        if (supply_by_rank[r].count(cid)) { any_supplier = true; break; }
      }
      if (!any_supplier) {
        hlog << "(TCR-diag) iter=" << hemocell.iter << " cid=" << cid
             << " (base=" << base_cell_id(cid) << ")"
             << " has NO supply on ANY rank — vertex permanently lost (cascade deletion)\n";
      }
    }
  }

  // Expand bounding box on all local blocks to accept particles from any domain position.
  // addParticle gates on getBoundingBox(); particle_grid update inside has its own
  // bounds guard so no crash risk from the expanded box.
  Box3D domain_bb = immersedParticles->getBoundingBox();
  plint expand = std::max({domain_bb.x1 - domain_bb.x0 + 1,
                           domain_bb.y1 - domain_bb.y0 + 1,
                           domain_bb.z1 - domain_bb.z0 + 1});
  map<plint, Box3D> saved_bb;
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    saved_bb[lbid] = pf.getBoundingBox();
    Box3D & bb = pf.getBoundingBox();
    bb = Box3D(bb.x0 - expand, bb.x1 + expand,
               bb.y0 - expand, bb.y1 + expand,
               bb.z0 - expand, bb.z1 + expand);
  }

  // Send: for each rank that needs something I supply, send those particles.
  // Send/receive sets are symmetric by construction (Allgather 1 and 2 agree),
  // so there is no MPI_Probe race and no deadlock.
  vector<vector<char>> send_bufs;
  vector<MPI_Request> send_reqs;
  for (int r = 0; r < nranks; r++) {
    if (r == me) continue;
    vector<int> to_send;
    for (const auto & kv : my_supply) {
      if (needed_by_rank[r].count(kv.first)) to_send.push_back(kv.first);
    }
    if (to_send.empty()) continue;
    vector<char> buf;
    for (int cid : to_send) {
      for (const auto & sv : my_supply[cid]) {
        size_t off = buf.size();
        buf.resize(off + sizeof(HemoCellParticle::serializeValues_t));
        memcpy(&buf[off], &sv, sizeof(HemoCellParticle::serializeValues_t));
      }
    }
    send_bufs.push_back(std::move(buf));
    send_reqs.push_back(MPI_REQUEST_NULL);
    MPI_Isend(send_bufs.back().data(), (int)send_bufs.back().size(), MPI_CHAR,
              r, 46, MPI_COMM_WORLD, &send_reqs.back());
  }

  // Determine which ranks will send to me (their supply intersects my needs)
  vector<int> supplier_ranks;
  if (!local_incomplete.empty()) {
    for (int r = 0; r < nranks; r++) {
      if (r == me) continue;
      for (int cid : local_incomplete) {
        if (supply_by_rank[r].count(cid)) { supplier_ranks.push_back(r); break; }
      }
    }
  }

  // Receive from each supplier and deliver to incomplete blocks
  for (int r : supplier_ranks) {
    MPI_Status status;
    MPI_Probe(r, 46, MPI_COMM_WORLD, &status);
    int count;
    MPI_Get_count(&status, MPI_CHAR, &count);
    vector<char> rbuf(count);
    MPI_Recv(rbuf.data(), count, MPI_CHAR, r, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int n_particles = count / (int)sizeof(HemoCellParticle::serializeValues_t);
    for (int pi = 0; pi < n_particles; pi++) {
      HemoCellParticle::serializeValues_t sv;
      memcpy(&sv, &rbuf[pi * sizeof(HemoCellParticle::serializeValues_t)], sizeof(sv));
      HemoCellParticle p(sv);
      for (const auto & kv : block_incomplete) {
        if (!kv.second.count(sv.cellId)) continue;
        HemoCellParticleField & pf = immersedParticles->getComponent(kv.first);
        pf.addParticle(&p);
      }
    }
  }

  // Same-rank self-supply: fill local incomplete blocks from local particle supply
  for (const auto & supply_kv : my_supply) {
    int cid = supply_kv.first;
    for (const auto & blk_kv : block_incomplete) {
      if (!blk_kv.second.count(cid)) continue;
      HemoCellParticleField & pf = immersedParticles->getComponent(blk_kv.first);
      for (const auto & sv : supply_kv.second) {
        HemoCellParticle p(sv);
        pf.addParticle(&p);
      }
    }
  }

  // Wait for all sends to complete
  if (!send_reqs.empty())
    MPI_Waitall((int)send_reqs.size(), send_reqs.data(), MPI_STATUSES_IGNORE);

  // Restore bounding boxes
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    immersedParticles->getComponent(lbid).getBoundingBox() = saved_bb[lbid];
  }
}

void HemoCellFields::HemoAdvanceParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->advanceParticles();
}
void HemoCellFields::advanceParticles() {
  global.statistics.getCurrent()["advanceParticles"].start();
    
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoAdvanceParticles(),immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSpreadParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->spreadParticleForce(domain);
}
void HemoCellFields::spreadParticleForce() {
  global.statistics.getCurrent()["spreadParticleForce"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSpreadParticleForce(),immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoApplyConstitutiveModel::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->applyConstitutiveModel(forced);
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->unifyForceVectors();
}
void HemoCellFields::unify_force_vectors() {
  global.statistics.getCurrent()["unifyForceVectors"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoUnifyForceVectors(),immersedParticles->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->applyRepulsionForce();
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->applyBoundaryRepulsionForce();
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->populateBoundaryParticles();
}
void HemoCellFields::populateBoundaryParticles() {
    vector<MultiBlock3D*>wrapper;
    wrapper.push_back(immersedParticles);
    HemoPopulateBoundaryParticles * fnct = new HemoPopulateBoundaryParticles();
    applyProcessingFunctional(fnct,immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoPopulateBindingSites::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->populateBindingSites(domain);
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->updateResidenceTime(rtime);
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
    dynamic_cast<HemoCellParticleField*>(blocks[0])->separateForceVectors();
}
void HemoCellFields::separate_force_vectors() {
  global.statistics.getCurrent()["separateForceVectors"].start();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSeperateForceVectors(),immersedParticles->getBoundingBox(),wrapper);

  global.statistics.getCurrent().stop();
}
void HemoCellFields::HemoDeleteIncompleteCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HemoCellParticleField*>(blocks[0])->deleteIncompleteCells(verbose);
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
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField*>(blocks[0]);
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
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField*>(blocks[0]);
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
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField*>(blocks[0]);
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
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField*>(blocks[0]);
  pf->solidifyCells();
}
void HemoCellFields::solidifyCells() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSolidifyCells(),immersedParticles->getBoundingBox(),wrapper);
}

void HemoCellFields::HemoPrepareSolidification::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField*>(blocks[0]);
  pf->prepareSolidification();
}
void HemoCellFields::prepareSolidification() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoPrepareSolidification(),immersedParticles->getBoundingBox(),wrapper);
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
HemoCellFields::HemoPrepareSolidification *        HemoCellFields::HemoPrepareSolidification::clone() const { return new HemoCellFields::HemoPrepareSolidification(*this);}
HemoCellFields::HemoPopulateBindingSites * HemoCellFields::HemoPopulateBindingSites::clone() const { return new HemoCellFields::HemoPopulateBindingSites(*this);}
HemoCellFields::HemoupdateResidenceTime * HemoCellFields::HemoupdateResidenceTime::clone() const { return new HemoCellFields::HemoupdateResidenceTime(*this);}


void HemoCellFields::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
   for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::hemocell;
   }
}

MultiParticleField3D<HemoCellParticleField> & HemoCellFields::getParticleField3D() { return *immersedParticles; };

}
