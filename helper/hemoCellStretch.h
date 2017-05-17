#ifndef HEMOCELLSTRETCH_H
#define HEMOCELLSTRETCH_H

#include "hemocell_internal.h"
#include "hemoCellFields.h"
#include "hemoCellFunctional.h"

class HemoCellStretch {
  public:

  HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, double external_force_);
  void applyForce();
  
  class FindForcedLsps: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   FindForcedLsps * clone() const;
  };
  
  class ForceForcedLsps: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   ForceForcedLsps * clone() const;
  };
  
  //Vertex id of lsps
  static vector<plint> lower_lsps;
  static vector<plint> upper_lsps;
  
  HemoCellField & cellfield;
  static unsigned int n_forced_lsps;
  static double external_force;
};

#endif
