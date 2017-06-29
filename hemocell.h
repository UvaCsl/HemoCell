#ifndef HEMOCELL_H
#define HEMOCELL_H

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
#include "writeCellInfoCSV.h"
#include "readPositionsBloodCells.h"

/* HELPERS */
#include "genericTools.h"
#include "meshMetrics.h"
#include "voxelizeDomain.h"
#include "meshGeneratingFunctions.h"
#include "loadBalancer.h"

/* MECHANICS */
#include "cellMechanics.h"
#include "constantConversion.h"

/* EXTERNALS */
#include "diagonalize.hpp"  // TODO: Do we need this file?


class HemoCell {
  public:
  /**
   * Creates an hemocell object
   * 
   * @param configFileName the location of the main config file
   * 
   * Unfortunately, due to palabos regulations, it is required to pass the
   * commandline arguments
   */
  
  HemoCell(char * configFileName, int argc, char* argv[]);

  /**
   *  Set all the fluid nodes to these values
   * 
   *  @param rho the desired density in lbm units
   *  @param vel the desired macroscopic velocity of each node
   */
  void latticeEquilibrium(double rho, Array<double, 3> vel);

  /**
   * Initialice the cellfields structure (and thus also the particlefield)
   */
  void initializeCellfield();

 /**
  * Add a celltype
  * valid options for constructType are:
  * RBC_FROM_SPHERE <- RBC
  * ELLIPSOID_FROM_SPHERE <- platelet
  * STRING_FROM_VERTEXES ->von willibrand factor
  * use as addCelltype<RbcHO>("RBC", RBC_FROM_SPHERE) for example
  * Since it is a template, it must be in the header class, maybe move to .hh
  * file for readability ...
  */
  template<class Mechanics>
  void addCellType(string name, int constructType) {
    string materialXML = name + ".xml";
    Config *materialCfg = new Config(materialXML.c_str());
    TriangularSurfaceMesh<double> * meshElement;
    
    double aspectRatio = 0.3;
    if (constructType == ELLIPSOID_FROM_SPHERE) {
      aspectRatio = (*materialCfg)["MaterialModel"]["aspectRatio"].read<double>();
    }
    
    if(constructType == STRING_FROM_VERTEXES) {
      meshElement = constructStringMeshFromConfig(*materialCfg);
    } else {
      TriangleBoundary3D<double> * boundaryElement = new TriangleBoundary3D<double>(constructMeshElement(constructType, 
                           (*materialCfg)["MaterialModel"]["radius"].read<double>()/param::dx, 
                           (*materialCfg)["MaterialModel"]["minNumTriangles"].read<double>(), param::dx, 
                           string(""), Array<double,3>(0.,0.,0.), aspectRatio));
      meshElement = new TriangularSurfaceMesh<double>(boundaryElement->getMesh());
    }

    HemoCellField * cellfield = cellfields->addCellType(*meshElement, name);
    Mechanics * mechanics = new Mechanics((*materialCfg), *cellfield);
    cellfield->mechanics = mechanics;
    mechanics->statistics();
  }

  /**
   * Set the output of a celltype
   * the outputs string should contain constants like VELOCITY_OUTPUT defined
   * in the constants_defaults file
   * 
   * @param outputs a vector of constants that define the desired output
   * @param name the name of the CellType ("RBC", "PLT")
   */
  void setOutputs(string name, vector<int> outputs);
  
  //Sets the repulsion constant and cutoff distance, also enables repulsion
  bool repulsionEnabled = false;
  void setRepulsion(double repulsionConstant, double repulsionCutoff);

  //Set the timescale separation of the particles of a particle type
  void setMaterialTimeScaleSeparation(string name, unsigned int separation);
  
  //Set the separation of when velocity is interpolated to the particle
  void setParticlePositionUpdateTimeScaleSeparation(unsigned int separation);

  //Set the timescale separation of the repulsion force for all particles
  void setRepulsionTimeScaleSeperation(unsigned int separation);
  
  //Set the minimum distance of the particles of a type to the solid, must be called BEFORE loadparticles
  void setMinimumDistanceFromSolid(string name, double distance);
  
  //Set the output of the fluid field
  void setFluidOutputs(vector<int> outputs);

  //Load the particles
private:
  bool loadParticlesIsCalled = false;
public:
  void loadParticles();

  void loadCheckPoint();
  void saveCheckPoint();

  bool outputInSiUnits = false;
  void writeOutput();
  void iterate();

  //Load balancing library functions
  double calculateFractionalLoadImbalance();
  void doLoadBalance();
  LoadBalancer * loadBalancer;

  MultiBlockLattice3D<double, DESCRIPTOR> * lattice;
	Config * cfg;
  HemoCellFields * cellfields;
  unsigned int iter = 0;

  private:
  XMLreader * documentXML; //Needed for legacy checkpoint reading TODO fix
  int lastOutputAt;  // Store the last time (iteration) output occured
};

#endif // HEMOCELL_H
