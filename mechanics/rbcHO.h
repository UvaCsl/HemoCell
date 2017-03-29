#ifndef HEMO_RbcHO
#define HEMO_RbcHO
#include "highOrderForces.h"
#include "constantConversion.h"
#include "config.h"

class RbcHO : public HighOrderForces {
  public:
  RbcHO(Config & cfg,HemoCellField & cellField_);
};
#endif
