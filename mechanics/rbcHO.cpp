#ifndef HEMO_RbcHO_cpp
#define HEMO_RbcHO_cpp
#include "rbcHO.h"

RbcHO::RbcHO(Config & cfg,HemoCellField & cellField_) :
                        HighOrderForces(cellField_,
                                        param::calculate_kVolume(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kArea(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kInPlane(cfg,"RBC",*cellField_.meshmetric),
                                        param::calculate_kBend(cfg,"RBC")) {};


#endif
