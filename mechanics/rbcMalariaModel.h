#ifndef HEMOCELL_RBCMALARIAMODEL_H
#define HEMOCELL_RBCMALARIAMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"

class RbcMalariaModel : public CellMechanics {

	public:
	//variables
	HemoCellField & cellField;
	const double k_volume;
	const double k_area;
	const double k_link;
	const double k_bend;
        const double k_inner_link;
	const double eta_m;
	const double eta_v;

	public:
	RbcMalariaModel(Config & modelCfg_, HemoCellField & cellField_);
	
	void ParticleMechanics(map<int,vector<HemoCellParticle *> > & particles_per_cell, const map<int, bool> &lpc, size_t ctype);
	
	void statistics();
	
	static double calculate_kBend(Config & cfg, MeshMetrics<double> &);
	static double calculate_kVolume(Config &cfg, MeshMetrics<double> &);
	static double calculate_kArea(Config &cfg, MeshMetrics<double> &);
	static double calculate_kLink(Config &cfg, MeshMetrics<double> &);
      	static double calculate_kInnerLink(Config &cfg, MeshMetrics<double> &);
	static double calculate_etaM(Config &cfg);
	static double calculate_etaV(Config &cfg);
};

#endif
