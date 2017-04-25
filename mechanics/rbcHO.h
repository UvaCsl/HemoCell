#ifndef HEMO_RbcHO
#define HEMO_RbcHO
#include "highOrderForces.h"
#include "constantConversion.h"
#include "config.h"

class RbcHO : public HighOrderForces {
  public:
  RbcHO(Config & cfg,HemoCellField & cellField_) :
                        HighOrderForces(cellField_,
                                        param::calculate_kVolume(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kArea(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kInPlane(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kBend(cfg,"RBC")) {};
};
#endif
