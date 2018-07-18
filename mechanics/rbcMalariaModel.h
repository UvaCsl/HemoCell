#ifndef HEMOCELL_RBCMALARIAMODEL_H
#define HEMOCELL_RBCMALARIAMODEL_H

#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"
namespace hemo {
class RbcMalariaModel : public CellMechanics {

	public:
	//variables
	HemoCellField & cellField;
	const T k_volume;
	const T k_area;
	const T k_link;
	const T k_bend;
        const T k_inner_link;
	const T eta_m;

	public:
	RbcMalariaModel(Config & modelCfg_, HemoCellField & cellField_);
	
	void ParticleMechanics(map<int,vector<HemoCellParticle *> > & particles_per_cell, const map<int, bool> &lpc, size_t ctype);
	
	void statistics();
	
      	static T calculate_kInnerLink(Config &cfg, MeshMetrics<T> &);
};

}
#endif
