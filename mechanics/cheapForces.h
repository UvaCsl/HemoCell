#ifndef HEMO_CHEAPFORCES_H
#define HEMO_CHEAPFORCES_H
#include "constantConversion.h"
#include "config.h"

class CheapForces : public CellMechanics {

  public:
  //Variables
  const CommonCellConstants cellConstants;
  HemoCellField & cellField;
  const double k_volume;
  const double k_area;
  const double k_inPlane;
  const double k_bend;


  //Constructor
  public:
  CheapForces(HemoCellField & cellField_, double k_volume_, double k_area_, double k_inPlane_, double k_bend_) : CellMechanics(),
                  cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellField_)),
                  cellField(cellField_), k_volume(k_volume_), k_area(k_area_), k_inPlane(k_inPlane_), k_bend(k_bend_)
  { }

  inline void ParticleMechanics(map<int,vector<HemoCellParticle *>> particles_per_cell, map<int,bool> lpc, pluint ctype) {
		//Do something cheap
  }

  inline void statistics() {
    pcout << "Cheap forces for " << cellField.name << " cellfield" << std::endl;
    pcout << "k_volume: " << k_volume << std::endl; 
    pcout << "k_area:   " << k_area << std::endl; 
    pcout << "k_inPlane:" << k_inPlane << std::endl; 
    pcout << "k_bend: : " << k_bend << std::endl; 
  }
};

class CheapForcesXML : public CheapForces {
  public:
  CheapForcesXML(Config & cfg,HemoCellField & cellField_, string name) :
                        CheapForces(cellField_,
                                        param::calculate_kVolume(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kArea(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kInPlane(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kBend(cfg,name)) {};
};
#endif
