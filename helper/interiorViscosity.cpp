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
#include "interiorViscosity.h"
#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"

namespace hemo {
  InteriorViscosityHelper::InteriorViscosityHelper(HemoCellFields & cellFields_) : cellFields(cellFields_) {
    //Create viscosity field with same properties as fluid field underlying the particleField.
    if(cellFields.hemocell.preInlet){
      preinlet_multiInteriorViscosityField = new plb::MultiScalarField3D<T>(
            MultiBlockManagement3D (
                *cellFields.hemocell.preinlet_lattice->getSparseBlockStructure().clone(),
                cellFields.hemocell.preinlet_lattice->getMultiBlockManagement().getThreadAttribution().clone(),
                cellFields.hemocell.preinlet_lattice->getMultiBlockManagement().getEnvelopeWidth(),
                cellFields.hemocell.preinlet_lattice->getMultiBlockManagement().getRefinementLevel()),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),                
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>(),
                0);
      preinlet_multiInteriorViscosityField->periodicity().toggle(0,cellFields.hemocell.preinlet_lattice->periodicity().get(0));
      preinlet_multiInteriorViscosityField->periodicity().toggle(1,cellFields.hemocell.preinlet_lattice->periodicity().get(1));
      preinlet_multiInteriorViscosityField->periodicity().toggle(2,cellFields.hemocell.preinlet_lattice->periodicity().get(2));

      preinlet_multiInteriorViscosityField->initialize();
      
      //Make sure each particleField has access to its local scalarField
      for (const plint & bId : preinlet_multiInteriorViscosityField->getLocalInfo().getBlocks()) {
        HemoCellParticleField & pf = cellFields.preinlet_immersedParticles->getComponent(bId);
        pf.interiorViscosityField = &preinlet_multiInteriorViscosityField->getComponent(bId);
      }
    }
    domain_multiInteriorViscosityField = new plb::MultiScalarField3D<T>(
            MultiBlockManagement3D (
                *cellFields.hemocell.domain_lattice->getSparseBlockStructure().clone(),
                cellFields.hemocell.domain_lattice->getMultiBlockManagement().getThreadAttribution().clone(),
                cellFields.hemocell.domain_lattice->getMultiBlockManagement().getEnvelopeWidth(),
                cellFields.hemocell.domain_lattice->getMultiBlockManagement().getRefinementLevel()),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),                
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>(),
                0);
    domain_multiInteriorViscosityField->periodicity().toggle(0,cellFields.hemocell.domain_lattice->periodicity().get(0));
    domain_multiInteriorViscosityField->periodicity().toggle(1,cellFields.hemocell.domain_lattice->periodicity().get(1));
    domain_multiInteriorViscosityField->periodicity().toggle(2,cellFields.hemocell.domain_lattice->periodicity().get(2));

    domain_multiInteriorViscosityField->initialize();
    
    //Make sure each particleField has access to its local scalarField
    for (const plint & bId : domain_multiInteriorViscosityField->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = cellFields.domain_immersedParticles->getComponent(bId);
      pf.interiorViscosityField = &domain_multiInteriorViscosityField->getComponent(bId);
    }

    if(cellFields.hemocell.partOfpreInlet){
      multiInteriorViscosityField = preinlet_multiInteriorViscosityField;
    }
    else{
      multiInteriorViscosityField = domain_multiInteriorViscosityField;
    }
    
  }
  
  InteriorViscosityHelper::~InteriorViscosityHelper() {
    delete multiInteriorViscosityField;
  }
    
  void InteriorViscosityHelper::checkpoint() {
    if (!global.enableInteriorViscosity) { 
      pcout << "(interiorViscosityField) Checkpoint called while global.enableSolidifyMechanics is not enabled, still checkpointing but this should not happen" << endl;
      return;
    }
    std::string & outDir = hemo::global.checkpointDirectory;
    mkpath(outDir.c_str(), 0777);
    
    if (global::mpi().isMainProcessor()) {
        renameFileToDotOld(outDir + "internalViscosity.dat");
        renameFileToDotOld(outDir + "internalViscosity.plb");
        if(cellFields.hemocell.preInlet){
          renameFileToDotOld(outDir + "PRE_internalViscosity.dat");
          renameFileToDotOld(outDir + "PRE_internalViscosity.plb");
        }
    }
    if(cellFields.hemocell.preInlet){
      plb::parallelIO::save(*preinlet_multiInteriorViscosityField, outDir + "PRE_internalViscosity", true);
    }
    plb::parallelIO::save(*domain_multiInteriorViscosityField, outDir + "internalViscosity", true);
  }
  
  void InteriorViscosityHelper::restore(HemoCellFields & cellFields) {
    if (!global.enableInteriorViscosity) { 
      pcout << "(internalViscosityField) Restore called while global.enableInteriorViscosity is not enabled, not restoring" << endl;
      return;
    }
    std::string & outDir = hemo::global.checkpointDirectory;
    std::string file_dat = outDir + "internalViscosity.dat";
    std::string file_plb = outDir + "internalViscosity.plb";
    if(!(file_exists(file_dat) && file_exists(file_plb))) {
      pcout << "(internalViscosityField) Error restoring internalViscosity fields from checkpoint, they do not seem to exist" << endl;
      exit(1);
    }
    if(cellFields.hemocell.preInlet){
      std::string file_dat = outDir + "PRE_internalViscosity.dat";
      std::string file_plb = outDir + "PRE_internalViscosity.plb";
      if(!(file_exists(file_dat) && file_exists(file_plb))) {
        pcout << "(PRE_internalViscosityField) Error restoring PRE_internalViscosity fields from checkpoint, they do not seem to exist" << endl;
        exit(1);
      }
    }
    plb::parallelIO::load(outDir + "internalViscosity",*get(cellFields).domain_multiInteriorViscosityField,true);
    if(cellFields.hemocell.preInlet){
      plb::parallelIO::load(outDir + "PRE_internalViscosity",*get(cellFields).preinlet_multiInteriorViscosityField,true);
    }
    get(cellFields).refillBindingSites();
  }  
  
  void InteriorViscosityHelper::add(HemoCellParticleField & pf, const Dot3D & internalPoint, T tau) {
      pf.internalPoints.insert(internalPoint);
      pf.interiorViscosityField->get(internalPoint.x,internalPoint.y,internalPoint.z) = tau;
  }
  
  void InteriorViscosityHelper::add(HemoCellParticleField & pf, const vector<Dot3D> & internalPoints, T tau) {
    for (const Dot3D & internalPoint : internalPoints) {
      add(pf,internalPoint,tau);
    }
  }
  
  void InteriorViscosityHelper::remove(HemoCellParticleField & pf, const Dot3D & internalPoint) {
    pf.internalPoints.erase(internalPoint);
    pf.interiorViscosityField->get(internalPoint.x,internalPoint.y,internalPoint.z) = 0;
  }
  
  void InteriorViscosityHelper::remove(HemoCellParticleField & pf, const vector<Dot3D> & internalPoints) {
    for (const Dot3D & internalPoint: internalPoints) {
      remove(pf,internalPoint);
    }
  } 
  
  void InteriorViscosityHelper::empty(HemoCellParticleField & pf) {
    for (const Dot3D & internalPoint: pf.internalPoints) {
      pf.interiorViscosityField->get(internalPoint.x,internalPoint.y,internalPoint.z) = 0;
    }
    pf.internalPoints.clear();
  }
  
  void InteriorViscosityHelper::refillBindingSites() {
    if(cellFields.hemocell.preInlet){
      for (const plint & bId : cellFields.preinlet_immersedParticles->getLocalInfo().getBlocks()) {
        HemoCellParticleField & pf = cellFields.preinlet_immersedParticles->getComponent(bId);
        ScalarField3D<T> & bf = *pf.interiorViscosityField;
        Box3D domain = bf.getBoundingBox();
        for (int x = domain.x0; x <= domain.x1 ; x++) {
          for (int y = domain.y0; y <= domain.y1; y++) {
            for (int z = domain.z0; z <= domain.z1; z++) {
              if(bf.get(x,y,z)) {
                pf.internalPoints.insert({x,y,z});
                
                //WARNING this _can_ memory leak, so you should be ok if this is only called one time (from checkpointing)
                plb::Dynamics<T,DESCRIPTOR>* dynamic = cellFields.hemocell.preinlet_lattice->getBackgroundDynamics().clone();
                dynamic->setOmega(1.0/bf.get(x,y,z));
                pf.atomicLattice->get(x,y,z).attributeDynamics(dynamic);
              }
            }
          }
        }
      }
    }
    for (const plint & bId : cellFields.domain_immersedParticles->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = cellFields.domain_immersedParticles->getComponent(bId);
      ScalarField3D<T> & bf = *pf.interiorViscosityField;
      Box3D domain = bf.getBoundingBox();
      for (int x = domain.x0; x <= domain.x1 ; x++) {
        for (int y = domain.y0; y <= domain.y1; y++) {
          for (int z = domain.z0; z <= domain.z1; z++) {
            if(bf.get(x,y,z)) {
              pf.internalPoints.insert({x,y,z});
              
              //WARNING this _can_ memory leak, so you should be ok if this is only called one time (from checkpointing)
              plb::Dynamics<T,DESCRIPTOR>* dynamic = cellFields.hemocell.domain_lattice->getBackgroundDynamics().clone();
              dynamic->setOmega(1.0/bf.get(x,y,z));
              pf.atomicLattice->get(x,y,z).attributeDynamics(dynamic);
            }
          }
        }
      }
    }
  }
}