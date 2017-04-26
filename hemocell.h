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
#include "meshGeneratingFunctions.h"

/* MECHANICS */
#include "cellMechanics.h"
#include "constantConversion.h"

/* EXTERNALS */
#include "diagonalize.hpp"
#include "external/readPositionsBloodCells.h"

class HemoCell {
  public:
  //Unfortunately, due to palabos regulations, it is required to pass the
  //commandline arguments
  HemoCell(char * configFileName, int argc, char* argv[]);

  //Set all the fluid nodes to these values
  void latticeEquilibrium(double rho, Array<double, 3> vel);

  //Initialice the cellfields structure (and thus also the particlefield)
  void initializeCellfield();

 /* Add a celltype
  * valid options for constructType are:
  * RBC_FROM_SPHERE <- RBC
  * ELLIPSOID_FROM_SPHERE <- platelet
  * use as addCelltype<RbcHO>("RBC", RBC_FROM_SPHERE) for example
  * Since it is a template, it must be in the header class, maybe move to .hh
  * file for readability ...
  */
  template<class Mechanics>
  void addCellType(string name, int constructType) {
    double aspectRatio = 0.3;
    if (constructType == ELLIPSOID_FROM_SPHERE) {
      aspectRatio = (*cfg)["MaterialModel"][name]["aspectRatio"].read<double>();
    }
    TriangleBoundary3D<double> * boundaryElement = new TriangleBoundary3D<double>(constructMeshElement(constructType, 
                           (*cfg)["ibm"]["radius"].read<double>()/param::dx, 
                           (*cfg)["MaterialModel"][name]["minNumTriangles"].read<double>(), param::dx, 
                           string(""), Array<double,3>(0.,0.,0.),aspectRatio));
    TriangularSurfaceMesh<double>  *meshElement = new TriangularSurfaceMesh<double>(boundaryElement->getMesh());

    //TODO correctly give hematocrit
    HemoCellField * cellfield = cellfields->addCellType(*meshElement, 0., name);
    Mechanics * mechanics = new Mechanics((*cfg), *cellfield);
    cellfield->mechanics = mechanics;
  }

  //Set the output of a celltype
  void setOutputs(string name, vector<int> outputs);
  //Set the output of the fluid field
  void setFluidOutputs(vector<int> outputs);

  void loadParticles(string name);

  void loadCheckPoint();
  void saveCheckPoint();

  void writeOutput();
  void iterate();

  MultiBlockLattice3D<double, DESCRIPTOR> * lattice;
	Config * cfg;
  HemoCellFields * cellfields;
  unsigned int iter = 0;
  private:
  XMLreader * documentXML; //Needed for legacy checkpoint reading TODO fix
};

#endif // HEMOCELL_H
