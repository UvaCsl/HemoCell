#ifndef HEMOCELL_H
#define HEMOCELL_H

#define VERSION_MAJOR 0
#define VERSION_MINOR 1

//Load Constants
#include "constant_defaults.h"
#include "hemocell_internal.h"
#include "config.h"

/* CORE libs */
#include "hemoCellFunctional.h"
#include "hemoCellParticle.h"
#include "hemoCellFields.h"
#include "hemoCellParticleType.h"

/* IO */
#include "ParticleHdf5IO.h"
#include "FluidHdf5IO.h"

/* HELPERS */
#include "genericTools.h"
#include "meshMetrics.h"
#include "voxelizeDomain.h"

/* MECHANICS */
#include "cellMechanics.h"
#include "constantConversion.h"

/* EXTERNALS */
#include "diagonalize.hpp"
#include "external/readPositionsBloodCells.h"

class HemoCell {
  HemoCell(string configFileName);
};

#endif // HEMOCELL_H
