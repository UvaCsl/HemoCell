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

#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "core/globalDefs.h"
#include "multiBlock/multiDataProcessorWrapper3D.hh"

namespace hemo {

class HemoCellFunctional : public plb::BoxProcessingFunctional3D {
public:
    HemoCellFunctional(){}
    void getModificationPattern(std::vector<bool>& isWritten) const {
        for (pluint i = 0; i < isWritten.size(); i++) {
            isWritten[i] = false;       
        }
    }
    plb::BlockDomain::DomainT appliesTo() const { return plb::BlockDomain::bulk; }
    void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        for (pluint i = 0; i < modified.size(); i++) {
            modified[i] = plb::modif::nothing;       
        }
        
    }
};

template<class GatherType>
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
    union byteGatherType {
        unsigned char b[sizeof(IDandGatherType)];
        IDandGatherType g;
    };
public:
    //Numblocks should be larger than the maximum number of blocks per process
    //The total can be retrieved from Multiblock.getManagment.getsparseblock.getnumblocks
    HemoCellGatheringFunctional(std::map<int,GatherType> & gatherValues_) : 
      gatherValues(gatherValues_) {}
    void getModificationPattern(std::vector<bool>& isWritten) const {
        for (pluint i = 0; i < isWritten.size(); i++) {
            isWritten[i] = false;       
        }
    }
    plb::BlockDomain::DomainT appliesTo() const { return plb::BlockDomain::bulk; }
    void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        for (pluint i = 0; i < modified.size(); i++) {
            modified[i] = plb::modif::nothing;       
        }
        
    }
    static void gather(std::map<int,GatherType> & gatherValues) {
        int be = 0;
        int sendsize = sizeof(int) + sizeof(IDandGatherType)*gatherValues.size();
        std::vector<unsigned char> sendbuffer(sendsize);

        byteint local_number;
        local_number.i = gatherValues.size();
        for (unsigned int i = 0; i < sizeof(int) ; i++) {
            sendbuffer[be] = local_number.b[i];
            be++;
        }
        for (auto const & entry : gatherValues) {
            byteGatherType bg;
            bg.g.ID = entry.first;
            bg.g.g = entry.second;
            for (unsigned int i = 0; i < sizeof(IDandGatherType) ; i++) {
                sendbuffer[be] = bg.b[i];
                be++;
            }
        }
        
        std::vector<int> sendcounts(plb::global::mpi().getSize());
        MPI_Allgather(&sendsize,sizeof(int),MPI_BYTE,&sendcounts[0],sizeof(int),MPI_BYTE,MPI_COMM_WORLD);
              
        int receivesize = 0;
        std::vector<int> displacements(1,0);
        for (int & size : sendcounts) {
          receivesize += size;
          displacements.push_back(displacements.back() + size);
        }
        displacements.pop_back();
        std::vector<unsigned char> receivebuffer(receivesize);        
        
        MPI_Allgatherv(&sendbuffer[0],sendsize,MPI_BYTE,&receivebuffer[0],&sendcounts[0],&displacements[0],MPI_BYTE,MPI_COMM_WORLD);
        be = 0;
        
        for (int j = 0 ; j < plb::global::mpi().getSize() ; j++) {
            be = displacements[j];
            for (unsigned int i = 0; i < sizeof(int) ; i++) {
                local_number.b[i] = receivebuffer[be];
                be++;
            }
            for (int i = 0 ; i < local_number.i ; i++) {
                byteGatherType bg;
                for (unsigned int k = 0; k < sizeof(byteGatherType) ; k++) {
                    bg.b[k] = receivebuffer[be];
                    be++;
                }
                gatherValues[bg.g.ID] = bg.g.g;
            }
        }
        
    }
        
    //This map should be set in the processingGenericBlocks function and is local to the mpi processor;
    std::map<int,GatherType> & gatherValues; 
};

}
#endif /* HEMOCELLFUNCTIONAL_H */

