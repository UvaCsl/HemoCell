/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HemoCellFunctional.h
 * Author: vikko
 *
 * Created on 10 March 2017, 15:34
 */

#ifndef HEMOCELLFUNCTIONAL_H
#define HEMOCELLFUNCTIONAL_H

#include "hemocell_internal.h"
#include "mpi.h"

class HemoCellFunctional : public plb::BoxProcessingFunctional3D {
public:
    HemoCellFunctional(){}
    void getModificationPattern(std::vector<bool>& isWritten) const {
        for (pluint i = 0; i < isWritten.size(); i++) {
            isWritten[i] = false;       
        }
    }
    BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i = 0; i < modified.size(); i++) {
            modified[i] = modif::nothing;       
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
    HemoCellGatheringFunctional(map<int,GatherType> & gatherValues_) : 
      gatherValues(gatherValues_) {}
    void getModificationPattern(std::vector<bool>& isWritten) const {
        for (pluint i = 0; i < isWritten.size(); i++) {
            isWritten[i] = false;       
        }
    }
    BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i = 0; i < modified.size(); i++) {
            modified[i] = modif::nothing;       
        }
        
    }
    static void gather(map<int,GatherType> & gatherValues, int numblocks) {
        int be = 0;
        int sendsize = sizeof(int) + sizeof(IDandGatherType)*numblocks;
        int receivesize = sendsize*global::mpi().getSize();
        vector<unsigned char> sendbuffer;
        sendbuffer.resize(sendsize);
        vector<unsigned char> receivebuffer;
        receivebuffer.resize(receivesize);
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
            for (unsigned int i = 0; i < sizeof(byteGatherType) ; i++) {
                sendbuffer[be] = bg.b[i];
                be++;
            }
        }
        
        MPI_Allgather(&sendbuffer[0],sendsize,MPI_BYTE,&receivebuffer[0],sendsize,MPI_BYTE,MPI_COMM_WORLD);
        be = 0;
        
        for (int j = 0 ; j < global::mpi().getSize() ; j++) {
            be = j*sendsize;
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
    map<int,GatherType> & gatherValues; 
};

#endif /* HEMOCELLFUNCTIONAL_H */

