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
#include "hemoCellFunctional.h"
#include "multiBlock/multiBlockOperations3D.hh"

void applyTimedProcessingFunctional(BoxProcessingFunctional3D* functional,
                               Box3D domain, std::vector<MultiBlock3D*> multiBlocks) {
    executeTimedDataProcessor( BoxProcessorGenerator3D(functional, domain),
                          multiBlocks );
 
};

void executeTimedDataProcessor( DataProcessorGenerator3D const& generator,
                           std::vector<MultiBlock3D*> multiBlocks )
{
    MultiProcessing3D<DataProcessorGenerator3D const, DataProcessorGenerator3D >
        multiProcessing(generator, multiBlocks);
    std::vector<DataProcessorGenerator3D*> const& retainedGenerators = multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const& atomicBlockNumbers = multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D*> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock=0; iBlock<extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock]
                = &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
				// Start Timers
				for (AtomicBlock3D* block : extractedAtomicBlocks) { block->timer.start(); }
        // Delegate to the "AtomicBlock version" of executeDataProcessor.
        plb::executeDataProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks);
				// Stop Timers
				for (AtomicBlock3D* block : extractedAtomicBlocks) { block->timer.stop(); }
    }
    // In the "executeProcessor" version, envelopes are updated right here, because the processor
    //   has already been executed. This behavior is unlike the behavior of the "addInternalProcessor" version,
    //   where envelopes are updated from within the multi-block, after the execution of internal processors.
    multiProcessing.updateEnvelopesWhereRequired();
}