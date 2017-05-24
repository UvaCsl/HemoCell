#ifndef CELLFIELDS3D_H
#define CELLFIELDS3D_H

class HemoCellFields;
#include "hemocell_internal.h"
#include "hemoCellParticleField.h"
#include "genericFunctions.h"
#include "hemoCellFunctional.h"
#include "hemoCellParticleType.h"
#include <unistd.h>

class HemoCell;
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
class HemoCellFields
{
public:
    HemoCellFields(MultiBlockLattice3D<double, DESCRIPTOR> & lattice_, unsigned int particleEnvelopeWidth,HemoCell &);
    MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> & getParticleField3D();
    void createParticleField();
    ~HemoCellFields();
    virtual void advanceParticles();
    virtual void interpolateFluidVelocity();
    virtual void spreadParticleForce();

    HemoCellField * addCellType(TriangularSurfaceMesh<double> & meshElement, std::string name_);
	  void setParticleUpdateScheme (double _cellTimeStep=1.0); //For decoupled update schemes 

    void readPositionsCellFields(std::string particlePosFile);
    /*Checkpoint functions*/
    void copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer);
    void copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer);
    void load(XMLreader * documentXML, unsigned int & iter);
    void save(XMLreader * documentXML, unsigned int iter);
    void InitAfterLoadCheckpoint();

    //double getMaximumForce_Global() {return 0;}

    unsigned int size();

    HemoCellField * operator[](unsigned int index);
    HemoCellField * operator[](string name);


    //void setFluidExternalForce(double poiseuilleForce);

	MultiBlockLattice3D<double, DESCRIPTOR> * lattice;
  vector<int> desiredFluidOutputVariables;
  HemoCell & hemocell;
  vector<HemoCellField *> cellFields;
  pluint envelopeSize;
	MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * immersedParticles;
	//MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * reductionParticles;
  double cellTimeStep;
  //void synchronizeCellQuantities(SyncRequirements _dummy) {}
  void separate_force_vectors();
  void unify_force_vectors();
  class HemoSeperateForceVectors: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSeperateForceVectors * clone() const;
  };
  class HemoUnifyForceVectors: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoUnifyForceVectors * clone() const;
  };
  class HemoSpreadParticleForce: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSpreadParticleForce * clone() const;
  };
  class HemoInterpolateFluidVelocity: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoInterpolateFluidVelocity * clone() const;
  };
  class HemoAdvanceParticles: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoAdvanceParticles * clone() const;
  };
  class HemoApplyConstitutiveModel: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoApplyConstitutiveModel * clone() const;
  };
  class HemoRepulsionForce: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoRepulsionForce * clone() const;
  };
  void calculateRepulsionForce();
  class HemoDeleteIncompleteCells: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoDeleteIncompleteCells * clone() const;
  };
  void deleteIncompleteCells();
  virtual void applyConstitutiveModel();
  void syncEnvelopes();
  class HemoSyncEnvelopes: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   HemoSyncEnvelopes * clone() const;
   void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  };
  bool hemocellfunction = false; //true if we should allow things to communicate (under our sight, not palabos);
};

#endif
