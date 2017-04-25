#ifndef READ_POSISIONS_OF_BLOOD_CELLS_H
#define READ_POSISIONS_OF_BLOOD_CELLS_H

#include "hemocell_internal.h"
#include "hemoCellFields.h"
#include "config.h"

void readPositionsBloodCellField3D(HemoCellFields & cellFields, double dx, const char* positionsFileName, Config & cfg);

void getReadPositionsBloodCellsVector(Box3D realDomain,
                                           std::vector<TriangularSurfaceMesh<double>* > & meshes,
                                           std::vector<plint> & Np,
                                           std::vector<std::vector<Array<double,3> > > & positions,
                                           std::vector<std::vector<plint> > & cellIds,
                                           std::vector<std::vector<Array<double,3> > > & randomAngles,
                                           const char* positionsFileName, Config & cfg);

class ReadPositionsBloodCellField3D : public BoxProcessingFunctional3D
{
public:
    ReadPositionsBloodCellField3D (HemoCellFields & cellFields_, double dx_, const char* positionsFileName_, Config & cfg_)
            : cellFields(cellFields_), cfg(cfg_) {dx = dx_; positionsFileName=positionsFileName_;}
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ReadPositionsBloodCellField3D* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    const char* positionsFileName;
    double dx;
    HemoCellFields & cellFields;
    Config & cfg;
};

#endif
