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
#include "bindingField.h"
#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"

namespace hemo {
  bindingFieldHelper::bindingFieldHelper(HemoCellFields * cellFields_) : cellFields(*cellFields_) {
    //Check if function is not called too early
    if (!cellFields_) {
      pcout << "(BindingField) ERROR: Bindingfield is requested while cellfields is not initialized. Perhaps populateBindingSites is called before hemocell.initializeCellFields?" << endl;
      exit(1);
    }
    
    //Create bindingfield with same properties as fluid field underlying the particleField.
    multiBindingField = new plb::MultiScalarField3D<bool>(
            MultiBlockManagement3D (
                *cellFields.hemocell.lattice->getSparseBlockStructure().clone(),
                cellFields.hemocell.lattice->getMultiBlockManagement().getThreadAttribution().clone(),
                cellFields.hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth(),
                cellFields.hemocell.lattice->getMultiBlockManagement().getRefinementLevel()),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),                
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<bool>(),
                0);
    multiBindingField->periodicity().toggle(0,cellFields.hemocell.lattice->periodicity().get(0));
    multiBindingField->periodicity().toggle(1,cellFields.hemocell.lattice->periodicity().get(1));
    multiBindingField->periodicity().toggle(2,cellFields.hemocell.lattice->periodicity().get(2));

    multiBindingField->initialize();
    
    //Make sure each particleField has access to its local scalarField
    for (const plint & bId : multiBindingField->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = cellFields.immersedParticles->getComponent(bId);
      pf.bindingField = &multiBindingField->getComponent(bId);
    }
  }
  
  bindingFieldHelper::~bindingFieldHelper() {
    delete multiBindingField;
  }
    
  void bindingFieldHelper::checkpoint() {
    if (!global.enableSolidifyMechanics) { 
      pcout << "(BindingField) Checkpoint called while global.enableSolidifyMechanics is not enabled, still checkpointing but this should not happen" << endl;
      return;
    }
    std::string & outDir = hemo::global.checkpointDirectory;
    mkpath(outDir.c_str(), 0777);
    
    if (global::mpi().isMainProcessor()) {
        renameFileToDotOld(outDir + "bindingSites.dat");
        renameFileToDotOld(outDir + "bindingSites.plb");
    }
    
    plb::parallelIO::save(*multiBindingField, outDir + "bindingSites", true);
  }
  
  void bindingFieldHelper::restore(HemoCellFields & cellFields) {
    if (!global.enableSolidifyMechanics) { 
      pcout << "(BindingField) Restore called while global.enableSolidifyMechanics is not enabled, not restoring" << endl;
      return;
    }

    std::string & outDir = hemo::global.checkpointDirectory;
    std::string file_dat = outDir + "bindingSites.dat";
    std::string file_plb = outDir + "bindingSites.plb";
    if(!(file_exists(file_dat) && file_exists(file_plb))) {
      pcout << "(BindingField) Error restoring bindingSites fields from checkpoint, they do not seem to exist" << endl;
      exit(1);
    }
    
    plb::parallelIO::load(outDir + "bindingSites",*get(cellFields).multiBindingField,true);
    get(cellFields).refillBindingSites();
  }  
  
  void bindingFieldHelper::add(HemoCellParticleField & pf, const Dot3D & bindingSite) {
      pf.bindingSites.insert(bindingSite);
      pf.bindingField->get(bindingSite.x,bindingSite.y,bindingSite.z) = true;
  }
  
  void bindingFieldHelper::add(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites) {
    for (const Dot3D & bindingSite : bindingSites) {
      add(pf,bindingSite);
    }
  }
  
  void bindingFieldHelper::remove(HemoCellParticleField & pf, const Dot3D & bindingSite) {
    pf.bindingSites.erase(bindingSite);
    pf.bindingField->get(bindingSite.x,bindingSite.y,bindingSite.z) = false;
  }
  
  void bindingFieldHelper::remove(HemoCellParticleField & pf, const vector<Dot3D> & bindingSites) {
    for (const Dot3D & bindingSite: bindingSites) {
      remove(pf,bindingSite);
    }
  } 
  
  void bindingFieldHelper::refillBindingSites() {
    for (const plint & bId : cellFields.immersedParticles->getLocalInfo().getBlocks()) {
      HemoCellParticleField & pf = cellFields.immersedParticles->getComponent(bId);
      ScalarField3D<bool> & bf = *pf.bindingField;
      Box3D domain = bf.getBoundingBox();
      for (int x = domain.x0; x <= domain.x1 ; x++) {
        for (int y = domain.y0; y <= domain.y1; y++) {
          for (int z = domain.z0; z <= domain.z1; z++) {
            if(bf.get(x,y,z)) {
              pf.bindingSites.insert({x,y,z});
            }
          }
        }
      }
    }
  }
}