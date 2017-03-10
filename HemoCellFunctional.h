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

#include "palabos3D.h"

class HemoCellFunctional : public plb::BoxProcessingFunctional3D {
public:
    HemoCellFunctional(){}
    void getModificationPattern(std::vector<bool>& isWritten) const {
        for (pluint i = 0; i < isWritten.size(); i++) {
            isWritten[i] = true;       
        }
    }
    BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i = 0; i < modified.size(); i++) {
            modified[i] = modif::dynamicVariables;       
        }
        
    }
};

#endif /* HEMOCELLFUNCTIONAL_H */

