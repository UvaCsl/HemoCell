#ifndef CELL_STRETCHING_3D_H
#define CELL_STRETCHING_3D_H

#include "cellStretchingForces3D.h"
#include "immersedCellParticleVtk3D.h"


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
class CellStretching3D {
public:
        CellStretching3D(
                TriangleBoundary3D<T> const& Cells_,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles_,
                plint numParticlesPerSide_, plint flowType_,
                T dx_, T dt_, T dNewton_,
                std::map<plint, Particle3D<T,Descriptor>*> * tagToParticle3D_,
                bool checkpointed_=false,
                T stretchForceScalarLU_=0, plint timesToStretch_=40,
                plint firstPlane_=TFL_DIRECTION_X, plint secondPlane_=TFL_DIRECTION_Z);
        virtual ~CellStretching3D();
        void applyForce(plint iter, T cellDensity) ;

/*   hasConverged Return types:
 *      0 -- It has NOT converged
 *      1 -- It HAS converged and simulation is over
 *      2 -- Previous force HAD converged, but moved to the next one.
 *      3 -- Cell has been released but has to recover the initial shape
 */
        plint hasConverged(plint iter);

        void writeConverged(plint iter, plint converged) ;
        void write(plint iter, T meanEdgeDistanceLU, T maxEdgeDistanceLU);
public:
        void setStretchScalarForce(T stretchForceScalar_);
        T getStretchScalarForce();
private:
        std::vector<plint> outerLeftTags, outerRightTags;
        std::vector<plint> outerFrontTags, outerBackTags;
        std::vector<T> stretchingDeformations;
        std::vector<std::vector<plint>*> lateralCellParticleTags;

        TriangleBoundary3D<T> const& Cells;
        MultiParticleField3D<ParticleFieldT<T,Descriptor> > * particles;
        plint numParticlesPerSide;
        plint flowType;
        T dx, dt, dNewton;
        std::map<plint, Particle3D<T,Descriptor>*> * tagToParticle3D;
        bool checkpointed;
        T stretchForceScalar;
        plint timesToStretch;
        plint firstPlane, secondPlane;
        plb_ofstream stretchLogFile, stretchResultFile, stretchReleasedFile;
        util::ValueTracer<T> convergeX, convergeY;

        T dStretchingForce;
        plint stretchReleased;
        plint checkInterval;

};





#include "cellStretching3D.hh"
#endif  // CELL_STRETCHING_3D_H

