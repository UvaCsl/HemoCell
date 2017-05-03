#ifndef HEMOCELL_FUNCTIONAL_CPP
#define HEMOCELL_FUNCTIONAL_CPP

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


#endif
