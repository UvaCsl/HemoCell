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
#include <numeric>

#include "hemoCellFields.h"
#include "hemocell.h"
#include "readPositionsBloodCells.h"
#include "constantConversion.h"
#include "bindingField.h"
#include "interiorViscosity.h"

#include "palabos3D.h"
#include "palabos3D.hh"


namespace hemo {


HemoCellFields::HemoCellFields(MultiBlockLattice3D<T, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth,
                               HemoCell & hemocell_)
    : lattice(&lattice_), hemocell(hemocell_) {
  envelopeSize = particleEnvelopeWidth;
  hlog << "(Hemocell) (HemoCellFields) (Init) particle envelope: " << particleEnvelopeWidth << " [lu]" << std::endl;
  if (hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth() < 2) {
    fatal("The fluid envelope is narrower than 2 nodes, which gives wrong forces at block boundaries.");
  }
  MPI_Comm_dup(MPI_COMM_WORLD, &particleComm);
  createParticleField();
  if (global.enableCEPACfield) { createCEPACfield(); }
  InitAfterLoadCheckpoint();
}

HemoCellFields::~HemoCellFields() {
  if (CEPACfield) { delete CEPACfield; }
  if (immersedParticles) { delete immersedParticles; }
  for (HemoCellField * field : cellFields) { delete field; }
  if (large_communicator) { delete large_communicator; }
  // After the particle fields, whose destructors may still use them.
  delete bindingHelper;
  delete interiorViscosityHelper;
  int mpiFinalized = 0;
  MPI_Finalized(&mpiFinalized);
  if (particleComm != MPI_COMM_NULL && !mpiFinalized) { MPI_Comm_free(&particleComm); }
}

void HemoCellFields::createParticleField(SparseBlockStructure3D * sbStructure, ThreadAttribution * tAttribution) {
  // A new particle field invalidates the cached large-envelope communication
  // structure (it refers to block IDs and ranks of the old field).
  // HemoCell::iterate() rebuilds it before the next step.
  if (large_communicator) {
    delete large_communicator;
    large_communicator = 0;
  }
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
    preinlet_immersedParticles = new MultiParticleField3D<HemoCellParticleField>(
        MultiBlockManagement3D(*hemocell.preinlet_lattice->getSparseBlockStructure().clone(),
                               hemocell.preinlet_lattice->getMultiBlockManagement().getThreadAttribution().clone(),
                               envelopeSize, refinement),
        plb::defaultMultiBlockPolicy3D().getCombinedStatistics());
  }
  domain_immersedParticles = new MultiParticleField3D<HemoCellParticleField>(
      MultiBlockManagement3D(*sbStructure, tAttribution, envelopeSize, refinement),
      plb::defaultMultiBlockPolicy3D().getCombinedStatistics());
  if (hemocell.partOfpreInlet) {
    immersedParticles = preinlet_immersedParticles;
  } else {
    immersedParticles = domain_immersedParticles;
  }
  InitAfterLoadCheckpoint();

  immersedParticles->periodicity().toggle(0, lattice->periodicity().get(0));
  immersedParticles->periodicity().toggle(1, lattice->periodicity().get(1));
  immersedParticles->periodicity().toggle(2, lattice->periodicity().get(2));

  immersedParticles->toggleInternalStatistics(false);

  if (del_sbStruct) { delete sbStructure; }

  InitAfterLoadCheckpoint();
}

void HemoCellFields::createCEPACfield() {
  SparseBlockStructure3D * sbStructure = lattice->getSparseBlockStructure().clone();
  ThreadAttribution * tAttribution = lattice->getMultiBlockManagement().getThreadAttribution().clone();
  plint refinement = lattice->getMultiBlockManagement().getRefinementLevel();
  lattice->getBlockCommunicator();
  CEPACfield = new MultiBlockLattice3D<T, CEPAC_DESCRIPTOR>(
      MultiBlockManagement3D(*sbStructure, tAttribution, envelopeSize, refinement),
      plb::defaultMultiBlockPolicy3D().getBlockCommunicator(), plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
      plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T, CEPAC_DESCRIPTOR>(),
      new plb::AdvectionDiffusionBGKdynamics<T, CEPAC_DESCRIPTOR>(param::tau_CEPAC));

  CEPACfield->periodicity().toggle(0, lattice->periodicity().get(0));
  CEPACfield->periodicity().toggle(1, lattice->periodicity().get(1));
  CEPACfield->periodicity().toggle(2, lattice->periodicity().get(2));

  CEPACfield->toggleInternalStatistics(false);

  integrateProcessingFunctional( // instead of integrateProcessingFunctional
      new LatticeToPassiveAdvDiff3D<T, DESCRIPTOR, CEPAC_DESCRIPTOR>(1), lattice->getBoundingBox(), *lattice,
      *CEPACfield, 1);
}

HemoCellField * HemoCellFields::addCellType(std::string name_, int constructType) {
  HemoCellField * cf = new HemoCellField(*this, name_, cellFields.size(), constructType);
  cellFields.push_back(cf);
  return cf;
}

unsigned int HemoCellFields::size() {
  return this->cellFields.size();
}

HemoCellField * HemoCellFields::operator[](unsigned int index) {
  if (index >= cellFields.size()) {
    fatal("Cell type index ", index, " requested, but there are only ", cellFields.size(), " cell types.");
  }
  return cellFields[index];
}

HemoCellField * HemoCellFields::operator[](string name) {
  for (unsigned int i = 0; i < cellFields.size(); i++) {
    if (cellFields[i]->name == name) { return cellFields[i]; }
  }
  fatal("Cell type \"", name, "\" requested, but it was not added with addCellType().");
  return NULL;
}

//SAVING FUNCTIONS
/* ******************* copyXMLreader2XMLwriter ***************************************** */
void HemoCellFields::copyXMLreader2XMLwriter(plb::XMLreader const & reader, plb::XMLwriter & writer) {
  std::string text = reader.getFirstText();
  if (!text.empty()) { writer[reader.getName()].setString(text); }
  std::vector<XMLreader *> const & children = reader.getChildren(reader.getFirstId());
  for (pluint iNode = 0; iNode < children.size(); ++iNode) {
    copyXMLreader2XMLwriter(*(children[iNode]), writer[reader.getName()]);
  }
}

void HemoCellFields::copyXMLreader2XMLwriter(plb::XMLreaderProxy readerProxy, plb::XMLwriter & writer) {
  std::string text;
  readerProxy.read(text);
  if (!text.empty()) { writer[readerProxy.getName()].setString(text); }
  std::vector<XMLreader *> const & children = readerProxy.getChildren();
  for (pluint iNode = 0; iNode < children.size(); ++iNode) {
    copyXMLreader2XMLwriter(*(children[iNode]), writer[readerProxy.getName()]);
  }
}

/*
 * Initialize variables that need to be loaded after a checkpoint AS WELL
 */
void HemoCellFields::InitAfterLoadCheckpoint() {
  if (!number_of_cells_from_checkpoint) { number_of_cells = getTotalNumberOfCells(*this); }
  std::vector<plint> const & blocks = immersedParticles->getLocalInfo().getBlocks();
  for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
    SmartBulk3D bulk(immersedParticles->getMultiBlockManagement(), blocks[iBlock]);
    Box3D blk = bulk.getBulk();
    immersedParticles->getComponent(blocks[iBlock]).setlocalDomain(blk);
    immersedParticles->getComponent(blocks[iBlock]).cellFields = this;
    immersedParticles->getComponent(blocks[iBlock]).atomicBlockId = blocks[iBlock];
    immersedParticles->getComponent(blocks[iBlock]).atomicLattice = &lattice->getComponent(blocks[iBlock]);
    if (global.enableCEPACfield && CEPACfield) {
      immersedParticles->getComponent(blocks[iBlock]).CEPAClattice = &CEPACfield->getComponent(blocks[iBlock]);
    }
    immersedParticles->getComponent(blocks[iBlock]).envelopeSize = envelopeSize;

    BlockLattice3D<T, DESCRIPTOR> * fluid = immersedParticles->getComponent(blocks[iBlock]).atomicLattice;
    immersedParticles->getComponent(blocks[iBlock]).nFluidCells = 0;
    for (unsigned int x = 0; x < fluid->getNx(); x++) {
      for (unsigned int y = 0; y < fluid->getNy(); y++) {
        for (unsigned int z = 0; z < fluid->getNz(); z++) {
          if (!fluid->get(x, y, z).getDynamics().isBoundary()) {
            immersedParticles->getComponent(blocks[iBlock]).nFluidCells++;
          }
        }
      }
    }

    //Calculate neighbours
    immersedParticles->getSparseBlockStructure().findNeighbors(
        blocks[iBlock], envelopeSize, immersedParticles->getComponent(blocks[iBlock]).neighbours);

    if (immersedParticles->getComponent(blocks[iBlock]).neighbours.size() > max_neighbours) {
      max_neighbours = immersedParticles->getComponent(blocks[iBlock]).neighbours.size();
    }
  }
}

void HemoCellFields::load(XMLreader * documentXML, unsigned int & iter, Config * cfg) {
  std::string firstField = (*(documentXML->getChildren(documentXML->getFirstId())[0])).getName();
  bool isCheckpointed = (firstField == "Checkpoint");
  if (isCheckpointed) {
    (*documentXML)["Checkpoint"]["General"]["Iteration"].read(iter);
    std::string outDir;
    (*documentXML)["Checkpoint"]["General"]["OutDirectory"].read(outDir);
    plb::global::directories().setOutputDir(outDir);
    loadDirectories(cfg, false);

    checkCheckPointFingerprint(*documentXML);

    // Restore state that cannot be re-derived from the input files. Older
    // checkpoints lack these entries (Palabos then throws); fall back to the
    // previous behaviour for those.
    try {
      int checkpointNumberOfCells = 0;
      (*documentXML)["Checkpoint"]["General"]["NumberOfCells"].read(checkpointNumberOfCells);
      number_of_cells = checkpointNumberOfCells;
      number_of_cells_from_checkpoint = true;
    } catch (plb::PlbIOException & e) {
      hlog
          << "(HemoCell) (CellFields) WARNING: checkpoint does not store NumberOfCells, deriving it from the .pos files; they must be unchanged since the original run."
          << endl;
    }
    try {
      unsigned int lastOutputIteration = 0;
      (*documentXML)["Checkpoint"]["General"]["LastOutputIteration"].read(lastOutputIteration);
      hemocell.setLastOutputIteration(lastOutputIteration);
    } catch (plb::PlbIOException & e) { hemocell.setLastOutputIteration(iter); }

    std::string & chkDir = hemo::global.checkpointDirectory;

    if (hemocell.preInlet) {
      plb::parallelIO::load(chkDir + "PRE_lattice", *hemocell.preinlet_lattice, true);
      plb::parallelIO::load(chkDir + "PRE_particleField", *preinlet_immersedParticles, true);
    }
    plb::parallelIO::load(chkDir + "lattice", *hemocell.domain_lattice, true);
    plb::parallelIO::load(chkDir + "particleField", *domain_immersedParticles, true);
    loadCEPACfield(chkDir);
  } else {
    pcout << "(HemoCell) (CellFields) loading checkpoint from non-checkpoint Config" << endl;
    std::string & chkDir = hemo::global.checkpointDirectory;
    if (hemocell.preInlet) {
      plb::parallelIO::load(chkDir + "PRE_lattice", *hemocell.preinlet_lattice, true);
      plb::parallelIO::load(chkDir + "PRE_particleField", *preinlet_immersedParticles, true);
    }
    plb::parallelIO::load(chkDir + "lattice", *hemocell.domain_lattice, true);
    plb::parallelIO::load(chkDir + "particleField", *domain_immersedParticles, true);
    loadCEPACfield(chkDir);
  }

  InitAfterLoadCheckpoint();
  syncEnvelopes();
  deleteIncompleteCells();
}

void HemoCellFields::save(XMLreader * xmlr, unsigned int iter, Config * cfg) {
  std::string & outDir = hemo::global.checkpointDirectory;

  mkpath(outDir.c_str(), 0777);


  /* Rename files, for safety reasons */
  if (global::mpi().isMainProcessor()) {
    renameFileToDotOld(outDir + "lattice.dat");
    renameFileToDotOld(outDir + "lattice.plb");
    renameFileToDotOld(outDir + "particleField.dat");
    renameFileToDotOld(outDir + "particleField.plb");
    renameFileToDotOld(outDir + "checkpoint.xml");
    if (global.enableCEPACfield && CEPACfield) {
      renameFileToDotOld(outDir + "CEPAC.dat");
      renameFileToDotOld(outDir + "CEPAC.plb");
    }
    if (hemocell.preInlet) {
      renameFileToDotOld(outDir + "PRE_lattice.dat");
      renameFileToDotOld(outDir + "PRE_lattice.plb");
      renameFileToDotOld(outDir + "PRE_particleField.dat");
      renameFileToDotOld(outDir + "PRE_particleField.plb");
    }
  }

  global::mpi().barrier();

  if (hemocell.preInlet) {
    plb::parallelIO::save(*hemocell.preinlet_lattice, outDir + "PRE_lattice", true);
    plb::parallelIO::save(*preinlet_immersedParticles, outDir + "PRE_particleField", true);
  }

  plb::parallelIO::save(*hemocell.domain_lattice, outDir + "lattice", true);
  plb::parallelIO::save(*domain_immersedParticles, outDir + "particleField", true);
  if (global.enableCEPACfield && CEPACfield) { plb::parallelIO::save(*CEPACfield, outDir + "CEPAC", true); }
}

namespace {
// MPI datatype of plint (ptrdiff_t: long or long long).
MPI_Datatype mpiPlintType() {
  static_assert(sizeof(plint) == sizeof(long) || sizeof(plint) == sizeof(long long), "no MPI datatype for plint");
  return sizeof(plint) == sizeof(long) ? MPI_LONG : MPI_LONG_LONG;
}

// Compile-time choices that change the raw particle records of checkpoints.
std::string particleLayout() {
  std::ostringstream layout;
  // No spaces: Palabos' XML reader reads strings up to the first space.
  layout << "record=" << sizeof(HemoCellParticle::serializeValues_t) << ",T=" << sizeof(T) << ",plint=" << sizeof(plint)
         << ",MATERIAL_INTEGRATION=" << HEMOCELL_MATERIAL_INTEGRATION;
#ifdef SOLIDIFY_MECHANICS
  layout << ",SOLIDIFY_MECHANICS";
#endif
  return layout.str();
}

std::string trimmed(const std::string & text) {
  const std::string::size_type first = text.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) { return ""; }
  return text.substr(first, text.find_last_not_of(" \t\n\r") - first + 1);
}

// Registered cell types in type-index order, "name:vertices" separated by commas.
std::string cellTypeList(HemoCellFields & fields) {
  std::ostringstream list;
  for (unsigned int i = 0; i < fields.size(); i++) {
    list << (i ? "," : "") << fields[i]->getIdentifier() << ":" << fields[i]->numVertex;
  }
  return list.str();
}
} // namespace

void HemoCellFields::saveCheckPointDescription(XMLreader * xmlr, unsigned int iter) {
  XMLwriter xmlw;
  std::string firstField = (*(xmlr->getChildren(xmlr->getFirstId())[0])).getName();
  if (firstField == "Checkpoint") {
    copyXMLreader2XMLwriter((*xmlr)["Checkpoint"]["hemocell"], xmlw["Checkpoint"]);
  } else {
    copyXMLreader2XMLwriter((*xmlr)["hemocell"], xmlw["Checkpoint"]);
  }
  xmlw["Checkpoint"]["General"]["Iteration"].set(iter);
  xmlw["Checkpoint"]["General"]["OutDirectory"].set(plb::global::directories().getOutputDir());
  // State that cannot be re-derived from the input files at restart.
  xmlw["Checkpoint"]["General"]["NumberOfCells"].set(number_of_cells);
  xmlw["Checkpoint"]["General"]["LastOutputIteration"].set(hemocell.getLastOutputIteration());
  // Build fingerprint, checked by load().
  xmlw["Checkpoint"]["General"]["ParticleLayout"].set(particleLayout());
  xmlw["Checkpoint"]["General"]["CellTypes"].set(cellTypeList(*this));
  xmlw["Checkpoint"]["General"]["dx"].set(param::dx, 17);
  xmlw["Checkpoint"]["General"]["dt"].set(param::dt, 17);

  const std::string file = hemo::global.checkpointDirectory + "checkpoint.xml";
  global::mpi().barrier(); // all data files are complete
  xmlw.print(file + ".tmp");
  if (global::mpi().isMainProcessor() && std::rename((file + ".tmp").c_str(), file.c_str()) != 0) {
    fatal("(Checkpoint) Cannot rename ", file, ".tmp to ", file, ".");
  }
  global::mpi().barrier();
}

void HemoCellFields::checkCheckPointFingerprint(XMLreader & documentXML) {
  std::string layout, types;
  T dx = 0, dt = 0;
  try {
    documentXML["Checkpoint"]["General"]["ParticleLayout"].read(layout);
    documentXML["Checkpoint"]["General"]["CellTypes"].read(types);
    documentXML["Checkpoint"]["General"]["dx"].read(dx);
    documentXML["Checkpoint"]["General"]["dt"].read(dt);
    layout = trimmed(layout);
    types = trimmed(types);
  } catch (plb::PlbIOException & e) {
    hlog << "(HemoCell) (CellFields) WARNING: the checkpoint has no build fingerprint (written by an older version); "
            "it cannot be checked that it matches this build and cell types."
         << endl;
    return;
  }
  if (layout != particleLayout()) {
    fatal("(Checkpoint) The checkpoint was written by a build with particle records \"", layout,
          "\", this build has \"", particleLayout(), "\"; use a library variant with the same compile-time features.");
  }
  if (types != cellTypeList(*this)) {
    fatal("(Checkpoint) The checkpoint has the cell types \"", types, "\" (name:vertices, in addCellType() order), ",
          "this run registers \"", cellTypeList(*this), "\".");
  }
  if (std::fabs(dx - param::dx) > 1e-12 * param::dx || std::fabs(dt - param::dt) > 1e-12 * param::dt) {
    fatal("(Checkpoint) The checkpoint has dx = ", dx, ", dt = ", dt, ", this run dx = ", param::dx,
          ", dt = ", param::dt, ".");
  }
}

void HemoCellFields::loadCEPACfield(const std::string & checkpointDirectory) {
  if (!(global.enableCEPACfield && CEPACfield)) { return; }
  // Checkpoints written before the CEPAC field was saved do not contain it.
  if (file_exists(checkpointDirectory + "CEPAC.plb")) {
    plb::parallelIO::load(checkpointDirectory + "CEPAC", *CEPACfield, true);
  } else {
    hlog << "(HemoCell) (CellFields) WARNING: checkpoint contains no CEPAC field, it keeps its initial state." << endl;
  }
}

void readPositionsCellFields(std::string particlePosFile) {}


//void HemoCellFields::

void HemoCellFields::HemoFindInternalParticleGridPoints::processGenericBlocks(Box3D domain,
                                                                              std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->findInternalParticleGridPoints(domain);
}

void HemoCellFields::findInternalParticleGridPoints() {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoFindInternalParticleGridPoints(), immersedParticles->getBoundingBox(), wrapper);
}


void HemoCellFields::HemoInternalGridPointsMembrane::processGenericBlocks(Box3D domain,
                                                                          std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->internalGridPointsMembrane(domain);
}

void HemoCellFields::internalGridPointsMembrane() {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoInternalGridPointsMembrane(), immersedParticles->getBoundingBox(), wrapper);
}


void HemoCellFields::HemoInterpolateFluidVelocity::processGenericBlocks(Box3D domain,
                                                                        std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->interpolateFluidVelocity(domain);
}
void HemoCellFields::interpolateFluidVelocity() {
  global.statistics.getCurrent()["interpolateFluidVelocity"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoInterpolateFluidVelocity(), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::calculateCommunicationStructure() {
  MultiBlockManagement3D management_temp(immersedParticles->getMultiBlockManagement());
  ParallelBlockCommunicator3D * communicator =
      dynamic_cast<ParallelBlockCommunicator3D const *>(&immersedParticles->getBlockCommunicator())->clone();
  communicator->duplicateOverlaps(management_temp, immersedParticles->periodicity());
  if (large_communicator) { delete large_communicator; }
  large_communicator = new CommunicationStructure3D(*communicator->communication);
  immersedParticles->getMultiBlockManagement().changeEnvelopeWidth(3);
  immersedParticles->signalPeriodicity();
  immersedParticles->getBlockCommunicator().duplicateOverlaps(*immersedParticles, modif::hemocell_no_comm);
  delete communicator;
}

// This replaces Palabos' duplicateOverlaps(modif::hemocell) on the particle
// field, which unpacked the messages of the partner ranks in arrival order.
// The unpack order decides the order of the particles in the receiving
// blocks, and with it the floating-point summation order of force spreading
// and repulsion, so multi-rank runs differed from run to run at round-off
// level. Here the messages are unpacked in ascending source rank.
void HemoCellFields::exchangeEnvelopes() {
  const ParallelBlockCommunicator3D & communicator =
      dynamic_cast<const ParallelBlockCommunicator3D &>(immersedParticles->getBlockCommunicator());
  // Palabos rebuilds its cached communication structure after the overlaps
  // changed (periodicity, envelope width); hemocell_no_comm moves no data.
  if (communicator.overlapsModified || !communicator.communication) {
    communicator.duplicateOverlaps(*immersedParticles, modif::hemocell_no_comm);
  }
  const CommunicationStructure3D & comms = *communicator.communication;
  // Tags on particleComm, distinct from those of the large-envelope exchange.
  const int tagSizes = 25;
  const int tagData = 43;
  const unsigned int particleSize = sizeof(HemoCellParticle::serializeValues_t);

  envelopeSendRanks.clear();
  for (const CommunicationInfo3D & info : comms.sendPackage) { envelopeSendRanks.push_back(info.toProcessId); }
  std::sort(envelopeSendRanks.begin(), envelopeSendRanks.end());
  envelopeSendRanks.erase(std::unique(envelopeSendRanks.begin(), envelopeSendRanks.end()), envelopeSendRanks.end());
  envelopeRecvRanks.clear();
  for (const CommunicationInfo3D & info : comms.recvPackage) { envelopeRecvRanks.push_back(info.fromProcessId); }
  std::sort(envelopeRecvRanks.begin(), envelopeRecvRanks.end());
  envelopeRecvRanks.erase(std::unique(envelopeRecvRanks.begin(), envelopeRecvRanks.end()), envelopeRecvRanks.end());
  auto slot = [](const vector<int> & ranks, int rank) {
    return unsigned(std::lower_bound(ranks.begin(), ranks.end(), rank) - ranks.begin());
  };
  const unsigned int nSend = envelopeSendRanks.size();
  const unsigned int nRecv = envelopeRecvRanks.size();

  // 1. Receive the message sizes. A partner sends one message per overlap,
  //    in the order of its send package, which matches our receive package.
  envelopeRecvSizes.resize(nRecv);
  envelopeRecvData.resize(nRecv);
  envelopeRecvInfos.resize(nRecv);
  for (unsigned int r = 0; r < nRecv; r++) {
    envelopeRecvSizes[r].clear();
    envelopeRecvInfos[r].clear();
  }
  for (const CommunicationInfo3D & info : comms.recvPackage) {
    const unsigned int r = slot(envelopeRecvRanks, info.fromProcessId);
    envelopeRecvInfos[r].push_back(&info);
    envelopeRecvSizes[r].push_back(0);
  }
  vector<MPI_Request> sizeRequests(nRecv), dataRequests(nRecv, MPI_REQUEST_NULL), sendRequests;
  for (unsigned int r = 0; r < nRecv; r++) {
    MPI_Irecv(envelopeRecvSizes[r].data(), int(envelopeRecvSizes[r].size()), MPI_INT, envelopeRecvRanks[r], tagSizes,
              particleComm, &sizeRequests[r]);
  }

  // 2. Send the particles in the bulk part of every overlap: per partner the
  //    message sizes, then all messages in one buffer.
  envelopeSendSizes.resize(nSend);
  envelopeSendData.resize(nSend);
  for (unsigned int s = 0; s < nSend; s++) {
    envelopeSendSizes[s].clear();
    envelopeSendData[s].clear();
  }
  vector<HemoCellParticle *> found;
  for (const CommunicationInfo3D & info : comms.sendPackage) {
    const unsigned int s = slot(envelopeSendRanks, info.toProcessId);
    immersedParticles->getComponent(info.fromBlockId).findParticles(info.fromDomain, found);
    vector<NoInitChar> & data = envelopeSendData[s];
    size_t offset = data.size();
    data.resize(offset + found.size() * particleSize);
    for (const HemoCellParticle * particle : found) {
      *((HemoCellParticle::serializeValues_t *)&data[offset]) = particle->sv;
      offset += particleSize;
    }
    envelopeSendSizes[s].push_back(int(found.size() * particleSize));
  }
  sendRequests.reserve(2 * nSend);
  for (unsigned int s = 0; s < nSend; s++) {
    sendRequests.emplace_back();
    MPI_Isend(envelopeSendSizes[s].data(), int(envelopeSendSizes[s].size()), MPI_INT, envelopeSendRanks[s], tagSizes,
              particleComm, &sendRequests.back());
    if (!envelopeSendData[s].empty()) {
      sendRequests.emplace_back();
      MPI_Isend(envelopeSendData[s].data(), int(envelopeSendData[s].size()), MPI_CHAR, envelopeSendRanks[s], tagData,
                particleComm, &sendRequests.back());
    }
  }

  // 3. Copies between blocks of this rank.
  for (const CommunicationInfo3D & info : comms.sendRecvPackage) {
    const AtomicBlock3D & fromBlock = immersedParticles->getComponent(info.fromBlockId);
    AtomicBlock3D & toBlock = immersedParticles->getComponent(info.toBlockId);
    toBlock.getDataTransfer().attribute(info.toDomain, info.fromDomain.x0 - info.toDomain.x0,
                                        info.fromDomain.y0 - info.toDomain.y0, info.fromDomain.z0 - info.toDomain.z0,
                                        fromBlock, modif::hemocell, info.absoluteOffset);
  }

  // 4. Post the data receives as the sizes arrive (in any order), ...
  for (unsigned int i = 0; i < nRecv; i++) {
    int r;
    MPI_Waitany(int(nRecv), sizeRequests.data(), &r, MPI_STATUS_IGNORE);
    const vector<int> & sizes = envelopeRecvSizes[r];
    const size_t total = std::accumulate(sizes.begin(), sizes.end(), size_t(0));
    envelopeRecvData[r].resize(total);
    if (total > 0) {
      MPI_Irecv(envelopeRecvData[r].data(), int(total), MPI_CHAR, envelopeRecvRanks[r], tagData, particleComm,
                &dataRequests[r]);
    }
  }
  // ... but unpack them in ascending source rank, the messages of a partner
  // in the order of the receive package.
  for (unsigned int r = 0; r < nRecv; r++) {
    MPI_Wait(&dataRequests[r], MPI_STATUS_IGNORE);
    char * data = reinterpret_cast<char *>(envelopeRecvData[r].data());
    size_t position = 0;
    for (unsigned int m = 0; m < envelopeRecvInfos[r].size(); m++) {
      const CommunicationInfo3D & info = *envelopeRecvInfos[r][m];
      const int length = envelopeRecvSizes[r][m];
      immersedParticles->getComponent(info.toBlockId)
          .getDataTransfer()
          .receive(info.toDomain, data + position, length, modif::hemocell, info.absoluteOffset);
      position += length;
    }
  }

  MPI_Waitall(int(sendRequests.size()), sendRequests.data(), MPI_STATUSES_IGNORE);
}

void HemoCellFields::HemoSyncEnvelopes::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->syncEnvelopes();
}
void HemoCellFields::syncEnvelopes() {
  global.statistics.getCurrent()["syncEnvelopes"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
    pf.removeParticles_inverse(pf.localDomain);
  }
  exchangeEnvelopes();

  if (large_communicator) {
    CommunicationStructure3D * comms = large_communicator;
    std::set<int> recv_procs, send_procs;
    std::map<int, vector<CommunicationInfo3D const *>> recv_infos, send_infos;
    for (CommunicationInfo3D const & info : comms->sendPackage) {
      send_procs.insert(info.toProcessId);
      send_infos[info.toProcessId].push_back(&info);
    }
    for (CommunicationInfo3D const & info : comms->recvPackage) {
      recv_procs.insert(info.fromProcessId);
      recv_infos[info.fromProcessId].push_back(&info);
    }

    set<plint> locals;
    for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = immersedParticles->getComponent(lbid);
      for (HemoCellParticle & particle : pf.particles) { locals.insert(particle.sv.cellId); }
    }
    vector<plint> locals_v;
    locals_v.insert(locals_v.end(), locals.begin(), locals.end());

    vector<int> recv_procs_v;
    recv_procs_v.insert(recv_procs_v.end(), recv_procs.begin(), recv_procs.end());
    vector<MPI_Request> reqs(recv_procs.size());

    // Message protocol on the private communicator `particleComm`:
    //  tag 24: list of cell IDs present on the requesting rank,
    //  tag 42: reply with all particles of those cells found in our blocks.
    // Requests are received per *known* source rank (not MPI_ANY_SOURCE): a
    // rank that finished this exchange early may already have sent its request
    // for the next call, and MPI's non-overtaking rule only orders messages
    // from the same source. Probing a specific source therefore guarantees we
    // consume the request belonging to this call. A reply (tag 42) can only be
    // sent in response to our own request, so replies cannot cross calls and
    // are received in arrival order.
    for (unsigned int i = 0; i < recv_procs_v.size(); i++) {
      MPI_Isend(locals_v.data(), locals_v.size(), mpiPlintType(), recv_procs_v[i], 24, particleComm, &reqs[i]);
    }
    sendBuffers.resize(send_procs.size());
    unsigned int send_index = 0;
    for (const int requesting_proc : send_procs) {
      MPI_Status status;
      MPI_Probe(requesting_proc, 24, particleComm, &status);
      int count;
      MPI_Get_count(&status, mpiPlintType(), &count);
      vector<plint> requested_ids(count);
      vector<NoInitChar> & sendBuffer = sendBuffers[send_index++];
      sendBuffer.clear();
      size_t offset = 0;
      MPI_Recv(requested_ids.data(), count, mpiPlintType(), requesting_proc, 24, particleComm, MPI_STATUS_IGNORE);
      for (CommunicationInfo3D const * info : send_infos[requesting_proc]) {
        HemoCellParticleField & pf = immersedParticles->getComponent(info->fromBlockId);
        // The requester knows the cells by the IDs of its periodic images.
        const plint offset_p = pf.getDataTransfer().getOffset(info->absoluteOffset);
        const map<plint, vector<int>> & ppc = pf.get_particles_per_cell();

        for (const plint requested : requested_ids) {
          const plint id = requested - offset_p;
          if (ppc.find(id) == ppc.end()) { continue; }
          for (int pid : ppc.at(id)) {
            if (pid <= -1) { continue; }
            if (pid >= (int)pf.particles.size()) { continue; }
            sendBuffer.resize(sendBuffer.size() + sizeof(HemoCellParticle::serializeValues_t));
            *((HemoCellParticle::serializeValues_t *)&sendBuffer[offset]) = pf.particles[pid].sv;
            offset += sizeof(HemoCellParticle::serializeValues_t);
          }
        }
      }
      reqs.emplace_back();
      MPI_Isend(sendBuffer.data(), sendBuffer.size(), MPI_CHAR, requesting_proc, 42, particleComm, &reqs.back());
    }

    // Receives are posted in arrival order (each source sends exactly one
    // reply), but unpacked in ascending source rank (recv_procs_v is sorted):
    // the unpack order decides the order of particles in the receiving blocks
    // and with it the floating-point summation order of force spreading and
    // repulsion. Arrival order made multi-rank runs differ from run to run at
    // round-off level.
    vector<MPI_Request> recv_reqs(recv_procs.size());
    recvBuffers.resize(recv_procs.size());
    for (unsigned int i = 0; i < recv_procs.size(); i++) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, 42, particleComm, &status);
      int count;
      MPI_Get_count(&status, MPI_CHAR, &count);
      const unsigned int slot =
          std::lower_bound(recv_procs_v.begin(), recv_procs_v.end(), status.MPI_SOURCE) - recv_procs_v.begin();
      vector<NoInitChar> & recv_buffer = recvBuffers[slot];
      recv_buffer.resize(count);
      MPI_Irecv(recv_buffer.data(), count, MPI_CHAR, status.MPI_SOURCE, 42, particleComm, &recv_reqs[slot]);
    }

    for (unsigned int i = 0; i < recv_procs_v.size(); i++) {
      MPI_Status status;
      if (MPI_SUCCESS != MPI_Wait(&recv_reqs[i], &status)) {
        fatal("(syncEnvelopes) MPI_Wait failed with error ", status.MPI_ERROR, ".");
      }
      //Get Offsets and Destinations
      for (CommunicationInfo3D const * info : recv_infos[recv_procs_v[i]]) {
        HemoCellParticleField & toBlock = immersedParticles->getComponent(info->toBlockId);
        toBlock.getDataTransfer().receive(info->toDomain, recvBuffers[i], info->absoluteOffset);
      }
    }

    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    // 3. Local copies which require no communication.
    for (unsigned iSendRecv = 0; iSendRecv < comms->sendRecvPackage.size(); ++iSendRecv) {
      CommunicationInfo3D const & info = comms->sendRecvPackage[iSendRecv];
      AtomicBlock3D const & fromBlock = immersedParticles->getComponent(info.fromBlockId);
      AtomicBlock3D & toBlock = immersedParticles->getComponent(info.toBlockId);

      toBlock.getDataTransfer().attribute(info.toDomain, 0, 0, 0, fromBlock, modif::hemocell, info.absoluteOffset);
    }
  }
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoAdvanceParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->advanceParticles();
}
void HemoCellFields::advanceParticles() {
  global.statistics.getCurrent()["advanceParticles"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoAdvanceParticles(), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSpreadParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->spreadParticleForce(domain);
}
void HemoCellFields::spreadParticleForce() {
  global.statistics.getCurrent()["spreadParticleForce"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSpreadParticleForce(), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoApplyConstitutiveModel::processGenericBlocks(Box3D domain,
                                                                      std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->applyConstitutiveModel(forced);
}
void HemoCellFields::applyConstitutiveModel(bool forced) {
  global.statistics.getCurrent()["applyConstitutiveModel"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoApplyConstitutiveModel * fnct = new HemoApplyConstitutiveModel();
  fnct->forced = forced;
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoUnifyForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->unifyForceVectors();
}
void HemoCellFields::unify_force_vectors() {
  global.statistics.getCurrent()["unifyForceVectors"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoUnifyForceVectors(), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoRepulsionForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->applyRepulsionForce();
}
void HemoCellFields::applyRepulsionForce() {
  global.statistics.getCurrent()["repulsionForce"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoRepulsionForce * fnct = new HemoRepulsionForce();
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoBoundaryRepulsionForce::processGenericBlocks(Box3D domain,
                                                                      std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->applyBoundaryRepulsionForce(resetFirst);
}
void HemoCellFields::applyBoundaryRepulsionForce(bool resetFirst) {
  global.statistics.getCurrent()["boundaryRepulsionForce"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoBoundaryRepulsionForce * fnct = new HemoBoundaryRepulsionForce();
  fnct->resetFirst = resetFirst;
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoPopulateBoundaryParticles::processGenericBlocks(Box3D domain,
                                                                         std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->populateBoundaryParticles();
}
void HemoCellFields::populateBoundaryParticles() {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoPopulateBoundaryParticles * fnct = new HemoPopulateBoundaryParticles();
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);
}

void HemoCellFields::HemoPopulateBindingSites::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->populateBindingSites(domain);
}
void HemoCellFields::populateBindingSites(plb::Box3D * box) {
  //Initialize bindingField before entering functional
  bindingFieldHelper::get(*this);

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  Box3D domain;
  if (box) {
    domain = *box;
  } else {
    domain = immersedParticles->getBoundingBox();
  }
  HemoPopulateBindingSites * fnct = new HemoPopulateBindingSites();
  applyProcessingFunctional(fnct, domain, wrapper);
}
void HemoCellFields::HemoupdateResidenceTime::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->updateResidenceTime(rtime);
}
void HemoCellFields::updateResidenceTime(unsigned int rtime) {
  global.statistics.getCurrent()["updateResidenceTime"].start();
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoupdateResidenceTime * fnct = new HemoupdateResidenceTime();
  fnct->rtime = rtime;
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);
  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSeperateForceVectors::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->separateForceVectors();
}
void HemoCellFields::separate_force_vectors() {
  global.statistics.getCurrent()["separateForceVectors"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSeperateForceVectors(), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}
void HemoCellFields::HemoDeleteIncompleteCells::processGenericBlocks(Box3D domain,
                                                                     std::vector<AtomicBlock3D *> blocks) {
  dynamic_cast<HemoCellParticleField *>(blocks[0])->deleteIncompleteCells(verbose);
}
void HemoCellFields::deleteIncompleteCells(bool verbose) {
  global.statistics.getCurrent()["deleteIncompleteCells"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  HemoDeleteIncompleteCells * fnct = new HemoDeleteIncompleteCells();
  fnct->verbose = verbose;
  applyProcessingFunctional(fnct, immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}
void HemoCellFields::deleteBrokenCells(plint envelope, bool verbose) {
  global.statistics.getCurrent()["deleteBrokenCells"].start();
  for (plint lbid : immersedParticles->getLocalInfo().getBlocks()) {
    immersedParticles->getComponent(lbid).deleteBrokenCells(envelope, verbose);
  }
  global.statistics.getCurrent().stop();
}

HemoCellFields::CompleteCellsScope::CompleteCellsScope(HemoCellFields & fields_) : fields(fields_) {
  if (fields.completeCellsScopes++ > 0) { return; }
  for (plint lbid : fields.immersedParticles->getLocalInfo().getBlocks()) {
    fields.immersedParticles->getComponent(lbid).snapshotParticles();
  }
  fields.syncEnvelopes();
  fields.deleteIncompleteCells(false);
}

HemoCellFields::CompleteCellsScope::~CompleteCellsScope() {
  if (--fields.completeCellsScopes > 0) { return; }
  for (plint lbid : fields.immersedParticles->getLocalInfo().getBlocks()) {
    fields.immersedParticles->getComponent(lbid).restoreParticles();
  }
}

void HemoCellFields::HemoGetParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  Box3D localDomain;
  intersect(domain, pf->localDomain, localDomain);
  pf->findParticles(localDomain, particles);
}
void HemoCellFields::getParticles(vector<HemoCellParticle *> & particles, Box3D & domain) {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoGetParticles(particles), domain, wrapper);
}
void HemoCellFields::HemoSetParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  Box3D localDomain;
  intersect(domain, pf->localDomain, localDomain);
  for (HemoCellParticle & particle : particles) { pf->addParticle(&particle); }
}
void HemoCellFields::addParticles(vector<HemoCellParticle> & particles) {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSetParticles(particles), immersedParticles->getBoundingBox(), wrapper);
}


void HemoCellFields::HemoDeleteNonLocalParticles::processGenericBlocks(Box3D domain,
                                                                       std::vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  pf->removeParticles_inverse(pf->localDomain.enlarge(envelopeSize));
}
void HemoCellFields::deleteNonLocalParticles(int envelope) {
  global.statistics.getCurrent()["deleteNonLocalParticles"].start();

  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoDeleteNonLocalParticles(envelope), immersedParticles->getBoundingBox(), wrapper);

  global.statistics.getCurrent().stop();
}

void HemoCellFields::HemoSolidifyCells::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  pf->solidifyCells();
}
void HemoCellFields::solidifyCells() {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoSolidifyCells(), immersedParticles->getBoundingBox(), wrapper);
}

void HemoCellFields::HemoPrepareSolidification::processGenericBlocks(Box3D domain,
                                                                     std::vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  pf->prepareSolidification();
}
void HemoCellFields::prepareSolidification() {
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(immersedParticles);
  applyProcessingFunctional(new HemoPrepareSolidification(), immersedParticles->getBoundingBox(), wrapper);
}

HemoCellFields::HemoInternalGridPointsMembrane * HemoCellFields::HemoInternalGridPointsMembrane::clone() const {
  return new HemoCellFields::HemoInternalGridPointsMembrane(*this);
}
HemoCellFields::HemoFindInternalParticleGridPoints * HemoCellFields::HemoFindInternalParticleGridPoints::clone() const {
  return new HemoCellFields::HemoFindInternalParticleGridPoints(*this);
}
HemoCellFields::HemoSeperateForceVectors * HemoCellFields::HemoSeperateForceVectors::clone() const {
  return new HemoCellFields::HemoSeperateForceVectors(*this);
}
HemoCellFields::HemoUnifyForceVectors * HemoCellFields::HemoUnifyForceVectors::clone() const {
  return new HemoCellFields::HemoUnifyForceVectors(*this);
}
HemoCellFields::HemoSpreadParticleForce * HemoCellFields::HemoSpreadParticleForce::clone() const {
  return new HemoCellFields::HemoSpreadParticleForce(*this);
}
HemoCellFields::HemoInterpolateFluidVelocity * HemoCellFields::HemoInterpolateFluidVelocity::clone() const {
  return new HemoCellFields::HemoInterpolateFluidVelocity(*this);
}
HemoCellFields::HemoAdvanceParticles * HemoCellFields::HemoAdvanceParticles::clone() const {
  return new HemoCellFields::HemoAdvanceParticles(*this);
}
HemoCellFields::HemoApplyConstitutiveModel * HemoCellFields::HemoApplyConstitutiveModel::clone() const {
  return new HemoCellFields::HemoApplyConstitutiveModel(*this);
}
HemoCellFields::HemoSyncEnvelopes * HemoCellFields::HemoSyncEnvelopes::clone() const {
  return new HemoCellFields::HemoSyncEnvelopes(*this);
}
HemoCellFields::HemoRepulsionForce * HemoCellFields::HemoRepulsionForce::clone() const {
  return new HemoCellFields::HemoRepulsionForce(*this);
}
HemoCellFields::HemoBoundaryRepulsionForce * HemoCellFields::HemoBoundaryRepulsionForce::clone() const {
  return new HemoCellFields::HemoBoundaryRepulsionForce(*this);
}
HemoCellFields::HemoDeleteIncompleteCells * HemoCellFields::HemoDeleteIncompleteCells::clone() const {
  return new HemoCellFields::HemoDeleteIncompleteCells(*this);
}
HemoCellFields::HemoGetParticles * HemoCellFields::HemoGetParticles::clone() const {
  return new HemoCellFields::HemoGetParticles(*this);
}
HemoCellFields::HemoSetParticles * HemoCellFields::HemoSetParticles::clone() const {
  return new HemoCellFields::HemoSetParticles(*this);
}
HemoCellFields::HemoPopulateBoundaryParticles * HemoCellFields::HemoPopulateBoundaryParticles::clone() const {
  return new HemoCellFields::HemoPopulateBoundaryParticles(*this);
}
HemoCellFields::HemoDeleteNonLocalParticles * HemoCellFields::HemoDeleteNonLocalParticles::clone() const {
  return new HemoCellFields::HemoDeleteNonLocalParticles(*this);
}
HemoCellFields::HemoSolidifyCells * HemoCellFields::HemoSolidifyCells::clone() const {
  return new HemoCellFields::HemoSolidifyCells(*this);
}
HemoCellFields::HemoPrepareSolidification * HemoCellFields::HemoPrepareSolidification::clone() const {
  return new HemoCellFields::HemoPrepareSolidification(*this);
}
HemoCellFields::HemoPopulateBindingSites * HemoCellFields::HemoPopulateBindingSites::clone() const {
  return new HemoCellFields::HemoPopulateBindingSites(*this);
}
HemoCellFields::HemoupdateResidenceTime * HemoCellFields::HemoupdateResidenceTime::clone() const {
  return new HemoCellFields::HemoupdateResidenceTime(*this);
}


void HemoCellFields::HemoSyncEnvelopes::getTypeOfModification(std::vector<modif::ModifT> & modified) const {
  for (pluint i = 0; i < modified.size(); i++) { modified[i] = modif::hemocell; }
}

MultiParticleField3D<HemoCellParticleField> & HemoCellFields::getParticleField3D() {
  return *immersedParticles;
};

}
