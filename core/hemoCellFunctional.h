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
#ifndef HEMOCELLFUNCTIONAL_H
#define HEMOCELLFUNCTIONAL_H

#include "constant_defaults.h"
#include "mpi.h"
#include <cstring>

#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "core/globalDefs.h"
#include "multiBlock/multiDataProcessorWrapper3D.hh"

namespace hemo {

class HemoCellFunctional : public plb::BoxProcessingFunctional3D {
public:
  HemoCellFunctional() {}
  void getModificationPattern(std::vector<bool> & isWritten) const {
    for (pluint i = 0; i < isWritten.size(); i++) { isWritten[i] = false; }
  }
  plb::BlockDomain::DomainT appliesTo() const { return plb::BlockDomain::bulk; }
  void getTypeOfModification(std::vector<plb::modif::ModifT> & modified) const {
    for (pluint i = 0; i < modified.size(); i++) { modified[i] = plb::modif::nothing; }
  }
};

template <class GatherType>
class HemoCellGatheringFunctional : public plb::BoxProcessingFunctional3D {
  //Make the structures byte accessible for easy mpi serialization
  union byteint {
    unsigned char b[sizeof(int)];
    int i;
  };
  struct IDandGatherType {
    int ID;
    GatherType g;
  };

public:
  //Numblocks should be larger than the maximum number of blocks per process
  //The total can be retrieved from Multiblock.getManagment.getsparseblock.getnumblocks
  HemoCellGatheringFunctional(std::map<int, GatherType> & gatherValues_) : gatherValues(gatherValues_) {}
  void getModificationPattern(std::vector<bool> & isWritten) const {
    for (pluint i = 0; i < isWritten.size(); i++) { isWritten[i] = false; }
  }
  plb::BlockDomain::DomainT appliesTo() const { return plb::BlockDomain::bulk; }
  void getTypeOfModification(std::vector<plb::modif::ModifT> & modified) const {
    for (pluint i = 0; i < modified.size(); i++) { modified[i] = plb::modif::nothing; }
  }
  /**
   * @brief Combine the maps of all ranks on every rank. Collective.
   *
   * Afterwards every rank holds the entries of all ranks; for a key present on
   * several ranks the value of the highest rank wins. The data volume is the
   * total map size on every rank, so use gatherToRoot() when only one rank
   * needs the result (e.g. for writing output).
   */
  static void gather(std::map<int, GatherType> & gatherValues) { exchange(gatherValues, -1); }

  /**
   * @brief Combine the maps of all ranks on rank `root` only. Collective.
   *
   * On `root` the result equals that of gather(); the other ranks keep their
   * own entries unchanged.
   */
  static void gatherToRoot(std::map<int, GatherType> & gatherValues, int root = 0) { exchange(gatherValues, root); }

private:
  // root < 0: all ranks receive (Allgatherv); otherwise only `root` (Gatherv).
  static void exchange(std::map<int, GatherType> & gatherValues, int root) {
    int be = 0;
    int sendsize = sizeof(int) + sizeof(IDandGatherType) * gatherValues.size();
    std::vector<unsigned char> sendbuffer(sendsize);

    byteint local_number;
    local_number.i = gatherValues.size();
    for (unsigned int i = 0; i < sizeof(int); i++) {
      sendbuffer[be] = local_number.b[i];
      be++;
    }
    for (auto const & entry : gatherValues) {
      IDandGatherType bg;
      bg.ID = entry.first;
      bg.g = entry.second;
      // The records follow a 4-byte count, so they are not aligned:
      // copy bytes instead of dereferencing a cast pointer.
      std::memcpy(&sendbuffer[be], &bg, sizeof(IDandGatherType));
      be += sizeof(IDandGatherType);
    }

    const int nRanks = plb::global::mpi().getSize();
    const bool allReceive = root < 0;
    const bool receiving = allReceive || plb::global::mpi().getRank() == root;
    std::vector<int> sendcounts(nRanks);
    if (allReceive) {
      MPI_Allgather(&sendsize, 1, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    } else {
      MPI_Gather(&sendsize, 1, MPI_INT, sendcounts.data(), 1, MPI_INT, root, MPI_COMM_WORLD);
    }

    int receivesize = 0;
    std::vector<int> displacements(nRanks, 0);
    if (receiving) {
      for (int j = 0; j < nRanks; j++) {
        displacements[j] = receivesize;
        receivesize += sendcounts[j];
      }
    }
    std::vector<unsigned char> receivebuffer(receivesize);

    if (allReceive) {
      MPI_Allgatherv(sendbuffer.data(), sendsize, MPI_BYTE, receivebuffer.data(), sendcounts.data(),
                     displacements.data(), MPI_BYTE, MPI_COMM_WORLD);
    } else {
      MPI_Gatherv(sendbuffer.data(), sendsize, MPI_BYTE, receivebuffer.data(), sendcounts.data(), displacements.data(),
                  MPI_BYTE, root, MPI_COMM_WORLD);
    }
    if (!receiving) { return; }

    for (int j = 0; j < nRanks; j++) {
      be = displacements[j];
      for (unsigned int i = 0; i < sizeof(int); i++) {
        local_number.b[i] = receivebuffer[be];
        be++;
      }
      for (int i = 0; i < local_number.i; i++) {
        IDandGatherType bg;
        std::memcpy(&bg, &receivebuffer[be], sizeof(IDandGatherType));
        be += sizeof(IDandGatherType);
        gatherValues[bg.ID] = bg.g;
      }
    }
  }

public:
  //This map should be set in the processingGenericBlocks function and is local to the mpi processor;
  std::map<int, GatherType> & gatherValues;
};

}
#endif /* HEMOCELLFUNCTIONAL_H */
