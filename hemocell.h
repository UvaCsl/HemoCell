#ifndef HEMOCELL_H
#define HEMOCELL_H

//Load Constants
#include "constant_defaults.h"

#include "palabos3D.h"
#include "palabos3D.hh"

#define VERSION_MAJOR 0
#define VERSION_MINOR 1

using namespace plb;

#include <algorithm>
#include <cmath>
#include <stack>
#include <stdio.h>
#include <string>
#include <vector>
#include <limits>
#include <map>

using namespace std;

//#include "fcnGenericFunctions.h"
//#include "simpleProfiler.cpp"
//#include "config.cpp"
#include "config.h"

//#include "ficsionInit.h"
//#include "ficsionInit.hh"
/* CORE libs */

#include "HemoCellFunctional.h"
#include "meshMetrics.h"

#include "surfaceParticle3D.h"
#include "cell3D.h"
//#include "helperFunctionals.h"
//#include "checkPoint.h"
#include "genericTools.h"

//#include "cellCellForces3D.h"
//#include "proximityDynamics3D.h"
//#include "cellForceChecking3D.h"
//#include "cellReductionTypes.h"
//#include "computeCellForces3D.h"
//#include "immersedBoundaryMethod3D.h"
//#include "meshGeneratingFunctions.h"

/* IO */ //IO is FUBAR
//#include "hdf5IO.h"
//#include "CellHdf5IO.h"
//#include "ParticleField3DHdf5IO.h"
//#include "config.h"


/* MODELS */
//#include "shellModel3D.h"
//#include "cellModel3D.h"
//#include "shapeMemoryModel3D.h"
//#include "restModel3D.h"
//#include "intermediateModel3D.h"

#include "cellMechanics.h"

#include "cellFields3D.h"
#include "hemoCellField.h"

#include "ParticleHdf5IO.h"
#include "FluidHdf5IO.h"


/* EXTENSIONS */
//#include "cellStretching3D.h"
//#include "cellInShearFlow3D.h"
//#include "cellStretchingForces3D.h"
//#include "rbcDisaggregation.h"

/* HELPERS */
//#include "immersedCellParticleFunctional3D.h"


/* EXTERNALS */
#include "diagonalize.hpp"

/* Particle position initialisation method */
//#include "external/orderedPositionsMultipleCells.h"
//#include "external/randomPositionsMultipleCells.h"
#include "external/readPositionsBloodCells.h"

#include "constantConversion.h"

#endif // HEMOCELL_H
