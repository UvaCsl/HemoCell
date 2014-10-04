#ifndef FICSION_H
#define FICSION_H

#include "palabos3D.h"
#include "palabos3D.hh"
using namespace plb;
using namespace std;



#include <algorithm>
#include <cmath>
#include <stack>
#include <stdio.h>
#include <string>
#include <vector>
#include <limits>
#include <map>

#include "ficsionInit.h"
#include "ficsionInit.hh"
/* CORE libs */
#include "cellField3D.h"
#include "reductionParticle3D.h"
#include "immersedCellParticle3D.h"

#include "cellCellForces3D.h"
#include "cellForceChecking3D.h"
#include "cellReductionTypes.h"
#include "computeCellForces3D.h"
#include "immersedBoundaryMethod3D.h"
#include "meshGeneratingFunctions.h"
#include "hdf5IO.h"
#include "ParticleHdf5IO.h"


/* MODELS */
#include "shellModel3D.h"
#include "cellModel3D.h"
#include "shapeMemoryModel3D.h"


/* CASES */
//#include "cellStretching3D.h"
//#include "cellInShearFlow3D.h"
//#include "cellStretchingForces3D.h"
//#include "rbcDisaggregation.h"

/* HELPERS */
#include "immersedCellParticleFunctional3D.h"
#include "meshMetrics.h"


/* EXTERNALS */
#include "diagonalize.hpp"



#endif // FICSION_H
