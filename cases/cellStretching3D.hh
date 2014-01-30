#ifndef CELL_STRETCHING_3D_HH
#define CELL_STRETCHING_3D_HH

#include "cellStretchingForces3D.h"



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
        CellStretching3D<T,Descriptor,ParticleFieldT>::CellStretching3D(
                TriangleBoundary3D<T> const& Cells_,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles_,
                plint numParticlesPerSide_, plint flowType_,
                T dx_, T dt_, T dNewton_,
                std::map<plint, Particle3D<T,Descriptor>*> * tagToParticle3D_,
                bool checkpointed_,
                T stretchForceScalarLU_, plint timesToStretch_,
                plint firstPlane_, plint secondPlane_) :
        Cells(Cells_),
        particles(&particles_), numParticlesPerSide(numParticlesPerSide_), flowType(flowType_),
        dx(dx_), dt(dt_), dNewton(dNewton_), tagToParticle3D(tagToParticle3D_),
        checkpointed(checkpointed_),
        stretchForceScalar(stretchForceScalarLU_), timesToStretch(timesToStretch_),
        firstPlane(firstPlane_), secondPlane(secondPlane_),
        convergeX((T)1,(T)100,1.0e-5), convergeY((T)1,(T)100,1.0e-5)
{
/*
 *  [3] Cell Stretching,
 *  [4] Cell Stretching Analysis and release,
 *  [5] Fast Cell Stretching Analysis and release,
 *  [8] Fast Cell Stretching release
 */

    //PLB_PRECONDITION( npar == 1 && MPI::COMM_WORLD.Get_size() == 1 );
    lateralCellParticleTags.push_back(&outerLeftTags);
    lateralCellParticleTags.push_back(&outerRightTags);
    lateralCellParticleTags.push_back(&outerFrontTags);
    lateralCellParticleTags.push_back(&outerBackTags);

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles);
    applyProcessingFunctional (
        new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerLeftTags, &outerRightTags, firstPlane),
        particles->getBoundingBox(), particleArg );
    applyProcessingFunctional (
        new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerFrontTags, &outerBackTags, secondPlane),
        particles->getBoundingBox(), particleArg );

    dStretchingForce = (200.0e-12/dNewton) / timesToStretch; // 200pN
    if (flowType == 8) {
        dStretchingForce = stretchForceScalar; // Usually 7pN in LU
        timesToStretch = 1;
        convergeX.setEpsilon(1e-7);
        convergeY.setEpsilon(1e-7);
    }
    checkInterval = 10;
    stretchReleased = 0;
    if ( (flowType == 4) or (flowType == 5) ) {
        stretchForceScalar = 0.0;
    }

    std::ostream::openmode mode = std::ostream::out;
    if (not checkpointed) { mode = mode | std::ostream::trunc; }
    else { mode = mode | std::ostream::app; }


    for (pluint ipa = 0; ipa < outerLeftTags.size(); ++ipa) {
        pcout << ipa << " : " <<outerLeftTags[ipa] << ", " <<outerRightTags[ipa] << std::endl;
    }
    if ((flowType == 3) || (flowType == 4) || (flowType == 5) || (flowType == 8) ) {
        stretchLogFile.open((plb::global::directories().getLogOutDir() + "plbStretchDeformation.log").c_str(), mode);
        stretchResultFile.open((plb::global::directories().getLogOutDir() + "plbStretchResults.log").c_str(), mode);
        stretchReleasedFile.open((plb::global::directories().getLogOutDir() + "plbStretchRelease.log").c_str(), mode);
    }
    if (not checkpointed) {
        stretchLogFile << setprecision(20) << "# t [sec]; Force[N]; D_A [m]; D_T [m]; Mean Edge Distance [LU]; Max Edge Distance [LU]; " << std::endl;
    }

}



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
CellStretching3D<T,Descriptor,ParticleFieldT>::~CellStretching3D() {};



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellStretching3D<T,Descriptor,ParticleFieldT>::setStretchScalarForce(T scalarForce) {
        stretchForceScalar = scalarForce;
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
T CellStretching3D<T,Descriptor,ParticleFieldT>::getStretchScalarForce() {
        return stretchForceScalar;
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellStretching3D<T,Descriptor,ParticleFieldT>::applyForce(plint iter, T cellDensity)
{
    if ((flowType == 3) || (flowType == 4) || (flowType == 5)  || (flowType == 8)  ) {
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        applyProcessingFunctional ( // compute force applied on the some particles by the stretching force
                new ApplyStretchingForce3D<T,Descriptor>(outerLeftTags, outerRightTags, Array<T,3>(stretchForceScalar,0,0), cellDensity, tagToParticle3D),
                particles->getBoundingBox(), particleArg );
    }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint CellStretching3D<T,Descriptor,ParticleFieldT>::hasConverged(plint iter)
{
    /*
     * Return types:
     *      0 -- It has NOT converged
     *      1 -- It HAS converged and simulation is over
     *      2 -- Previous force HAD converged, but moved to the next one.
     *      3 -- Cell has been released but has to recover the initial shape
     */
    plint converged=0;
    std::string delim = " ; ";
    if (( (flowType == 4) or (flowType == 5)  || (flowType == 8) ) && (iter % checkInterval == 0)) {
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        std::vector<Array<T,3> > angles;
        applyProcessingFunctional (
            new MeasureCellStretchDeformation3D<T,Descriptor>(lateralCellParticleTags, &stretchingDeformations, &angles, tagToParticle3D),
            particles->getBoundingBox(), particleArg );
        convergeX.takeValue(stretchingDeformations[0],false);
        convergeY.takeValue(stretchingDeformations[1],false);

        if (convergeX.hasConverged() && convergeY.hasConverged() ) {
            converged = 2;
            writeConverged(iter, converged);
            convergeX.resetValues(); convergeY.resetValues();
            stretchForceScalar += dStretchingForce;
            if ( (stretchForceScalar>timesToStretch*dStretchingForce ) and stretchReleased) {
                converged=1;
            }
            else if ( (stretchForceScalar>timesToStretch*dStretchingForce ) and not stretchReleased) {
                converged=3;
                stretchReleased = 1;
                stretchForceScalar = 0;

                stretchReleasedFile << setprecision(20) << iter*dt
                        << delim << stretchForceScalar*dNewton
                        << delim << stretchingDeformations[0]*dx
                        << delim << stretchingDeformations[1]*dx
                        << delim << std::endl;
                pcout << "[FICSION]: Releasing RBC." << std::endl;

            }
        }
    }
    return converged;
}



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellStretching3D<T,Descriptor,ParticleFieldT>::writeConverged(plint iter, plint converged)
{
    /* Checks the convergence and prints output. */
    std::string delim = " ; ";
    pcout << "# StretchForce: " << stretchForceScalar*dNewton*1.0e12 << "pN converged in iteration " << iter << "." << std::endl;
    stretchResultFile << setprecision(20) << iter*dt
            << delim << stretchForceScalar*dNewton
            << delim << stretchingDeformations[0]*dx
            << delim << stretchingDeformations[1]*dx
            << std::endl;
    std::string stretchNameId = stretchReleased?"stretchR.":"stretch.";
    writeImmersedSurfaceVTK (
        Cells, *particles,
        global::directories().getOutputDir()+createFileName(stretchNameId, plint(stretchForceScalar*dNewton*1.0e12), 3)+".pN.vtk");

}



template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellStretching3D<T,Descriptor,ParticleFieldT>::write(plint iter, T meanEdgeDistanceLU, T maxEdgeDistanceLU)
{
    if ((flowType == 3) || (flowType == 4) || (flowType == 5) || (flowType == 8)  ) {
            std::string delim = " ; ";
            std::vector<MultiBlock3D*> particleArg;
            particleArg.push_back(particles);
            std::vector<Array<T,3> > angles;
            applyProcessingFunctional (
                new MeasureCellStretchDeformation3D<T,Descriptor>(lateralCellParticleTags, &stretchingDeformations, &angles, tagToParticle3D),
                particles->getBoundingBox(), particleArg );
            pcout << "StrDeform [m] [";
            for (pluint iDirection = 0; iDirection < stretchingDeformations.size(); ++iDirection) {
                pcout << "(" << iDirection << ", " << stretchingDeformations[iDirection]*dx << "), " ;
           } pcout << "]" << std::endl;
           stretchLogFile << setprecision(20) << iter*dt
                          << delim << stretchForceScalar*dNewton
                          << delim << stretchingDeformations[0]*dx
                          << delim << stretchingDeformations[1]*dx
                          << delim << meanEdgeDistanceLU
                          << delim << maxEdgeDistanceLU
                                                            << std::endl;
           if (stretchReleased) {
               stretchReleasedFile << setprecision(20) << iter*dt
                       << delim << stretchForceScalar*dNewton
                       << delim << stretchingDeformations[0]*dx
                       << delim << stretchingDeformations[1]*dx
                       << delim<< meanEdgeDistanceLU
                       << delim  << maxEdgeDistanceLU
                                                           << std::endl;
           }
    }
}


#endif  // CELL_STRETCHING_FORCES_3D_HG

