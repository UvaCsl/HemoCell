#ifndef CELLFIELDS3D_H
#define CELLFIELDS3D_H

class CellFields3D;
#include "palabos3D.h"
#include "constant_defaults.h"
#include "palabos3D.hh"
#include "cell3D.h"
#include "cellModel3D.h"
#include "cellMechanics.h"
#include "hemoCellParticleField3D.h"
#include "shapeMemoryModel3D.h"
#include "cellReductionTypes.h"
#include "cellFieldFunctionals3D.h"
#include "cellCellForces3D.h"
#include "fcnGenericFunctions.h"
#include "HemoCellFunctional.h"
#include "rbcHO.h"
#include <sys/stat.h>
#include <unistd.h>
#include <set>
#include <string>
#include <algorithm>    // std::sort


using namespace std;
using namespace plb;

class HemoCellField;

/*
 * This class can contain many cellfields, it keeps track of all the particles
 * in all the cellfields. The option exists to get exclusive access to a single
 * cellfield if necessary.`
 * 
 * A cellField should have an local update and a between cell update function
 *
 * The between cell interaction might be different for different type of
 * particles, they can be set in a special map
 *
 *
 * TODO: light and heavy (tracking cellfield arrays to particles or not)
 *
 */
class CellFields3D
{
public:
	CellFields3D(MultiBlockLattice3D<double, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth);
    MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & getParticleField3D() { return *immersedParticles; };
    virtual void advanceParticles();
    virtual void interpolateFluidVelocity();
    virtual void spreadParticleForce();

    HemoCellField * addCellType(TriangularSurfaceMesh<double> & meshElement, double hematocrit,  std::string name_);
	  void setParticleUpdateScheme (double _cellTimeStep=1.0); //For decoupled update schemes 

    void readPositionsCellFields(std::string particlePosFile);
    /*Checkpoint functions*/
    void copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer);
    void copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer);
    void load(XMLreader * documentXML, plint & iter);
    void save(XMLreader * documentXML, plint iter);
    void InitAfterLoadCheckpoint();

    double getMaximumForce_Global() {return 0;}

    unsigned int size();

    HemoCellField * operator[](unsigned int index);

    void setFluidExternalForce(double poiseuilleForce) {}

	MultiBlockLattice3D<double, DESCRIPTOR> & lattice;
  vector<HemoCellField *> cellFields;
  pluint envelopeSize;
	MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * immersedParticles;
	//MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * reductionParticles;
  double cellTimeStep;
  void synchronizeCellQuantities(SyncRequirements _dummy) {}
  void separate_force_vectors();
  void unify_force_vectors();
  class HemoSeperateForceVectors: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSeperateForceVectors * clone() const { return new HemoSeperateForceVectors(*this);}
  };
  class HemoUnifyForceVectors: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoUnifyForceVectors * clone() const { return new HemoUnifyForceVectors(*this);}
  };
  class HemoSpreadParticleForce: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSpreadParticleForce * clone() const { return new HemoSpreadParticleForce(*this);}
  };
  class HemoInterpolateFluidVelocity: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoInterpolateFluidVelocity * clone() const { return new HemoInterpolateFluidVelocity(*this);}
  };
  class HemoAdvanceParticles: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoAdvanceParticles * clone() const { return new HemoAdvanceParticles(*this);}
  };
  class HemoApplyConstitutiveModel: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoApplyConstitutiveModel * clone() const { return new HemoApplyConstitutiveModel(*this);}
  };
  virtual void applyConstitutiveModel();
  void syncEnvelopes();
  class HemoSyncEnvelopes: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSyncEnvelopes * clone() const { return new HemoSyncEnvelopes(*this);}
   void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
     for (pluint i = 0; i < modified.size(); i++) {
       modified[i] = modif::dynamicVariables;
   } }

  };
  bool hemocellfunction = false; //true if we should allow things to communicate (under our sight, not palabos);
};

/*contains information about one particular cellfield, structlike*/
class HemoCellField{
  static vector<int> default_output;
  public:

  HemoCellField(CellFields3D& cellFields_, Cell3D<double,DESCRIPTOR> cell3D_, TriangularSurfaceMesh<double>& meshElement_);
  double getVolumeFraction() { return hematocrit;}
  double hematocrit;
  ShellModel3D<double> * model;
  TriangularSurfaceMesh<double> & getMesh() { return meshElement;}
  std::string name;
  int ctype;
  int numVertex;
  bool outputTriangles = false;
  CellFields3D & cellFields;
  vector<int> desiredOutputVariables;
  Cell3D<double,DESCRIPTOR> & cell3D;
  vector<Array<plint,3>> triangle_list;
  TriangularSurfaceMesh<double> & meshElement;
  void(*kernelMethod)(BlockLattice3D<double,DESCRIPTOR> const&,SurfaceParticle3D*) = interpolationCoefficientsPhi2;
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleField3D() {return cellFields.immersedParticles;}
  MultiBlockLattice3D<double,DESCRIPTOR> * getFluidField3D() {return &(cellFields.lattice);}
  int getNumberOfCells_Global() {return 0;}
  std::string getIdentifier() {return name;}
  Box3D getBoundingBox() { return cellFields.immersedParticles->getBoundingBox(); }
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleArg() { return cellFields.immersedParticles; }
  std::map<plint, Cell3D<double,DESCRIPTOR>* > getCellIdToCell3D() { std::map<plint,Cell3D<double,DESCRIPTOR>* > tmp; return tmp ;}
  void synchronizeSyncRequirements(SyncRequirements _dummy) {}
  void setOutputVariables(const vector<int> &);
  CellMechanics * mechanics;
  void statistics();
  MeshMetrics<double> * meshmetric;
};

vector<int> HemoCellField::default_output ({OUTPUT_POSITION});
#endif
