#ifndef HEMOCELL_H
#define HEMOCELL_H

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
//#include <array>
#include <limits>
#include <map>

using namespace std;

// ============== Compile time options - Set these

/*
Choose kernel. 
Phi1 [1], Phi2 [2] - Default, Phi3 [3], Phi4 [4],  Phi4c [5]
*/
#define HEMOCELL_KERNEL 2

/*
Choose material model.
1 - Dao/Suresh model (+Fedosov 2010 improvements)
2 - HO model (Under testing)
*/
#define HEMOCELL_MATERIAL_MODEL 2

/*
Choose material integration method.
Euler [1], Adams-Bashforth [2]
*/
#define HEMOCELL_MATERIAL_INTEGRATION 1

/*
Choose collision operator for LBM.
[1] - BGK <- use this for dt < 0 settings to have tau = 1 suppress oscillations
[2] - MRT <- use this in every other case
*/
#define HEMOCELL_CFD_DYNAMICS 1


// ===================================

#include "fcnGenericFunctions.h"
#include "simpleProfiler.cpp"

#include "ficsionInit.h"
#include "ficsionInit.hh"
/* CORE libs */
#include "cellField3D.h"
#include "cellParticle3D.h"
#include "surfaceParticle3D.h"
#include "helperFunctionals.h"
#include "checkPoint.h"
#include "genericTools.h"

#include "cellCellForces3D.h"
#include "proximityDynamics3D.h"
#include "cellForceChecking3D.h"
#include "cellReductionTypes.h"
#include "computeCellForces3D.h"
#include "immersedBoundaryMethod3D.h"
#include "meshGeneratingFunctions.h"

/* IO */
#include "hdf5IO.h"
#include "ParticleHdf5IO.h"
#include "CellHdf5IO.h"
#include "ParticleField3DHdf5IO.h"


/* MODELS */
#include "shellModel3D.h"
#include "cellModel3D.h"
#include "shapeMemoryModel3D.h"
#include "restModel3D.h"
#include "intermediateModel3D.h"

/* EXTENSIONS */
#include "cellStretching3D.h"
//#include "cellInShearFlow3D.h"
//#include "cellStretchingForces3D.h"
//#include "rbcDisaggregation.h"

/* HELPERS */
#include "immersedCellParticleFunctional3D.h"
#include "meshMetrics.h"


/* EXTERNALS */
#include "diagonalize.hpp"

/* Particle position initialisation method */
#include "external/orderedPositionsMultipleCells.h"
#include "external/randomPositionsMultipleCells.h"
#include "external/readPositionsBloodCells.h"



#endif // HEMOCELL_H
