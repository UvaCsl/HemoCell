/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef VOXELIZEDOMAIN_H
#define VOXELIZEDOMAIN_H

#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "offLattice/triangleBoundary3D.hh"


#include "array.h"
#include "config.h"

#include <memory>

namespace hemo {
class CopyFromNeighbor : public plb::BoxProcessingFunctional3D_S<int> {
public:
  CopyFromNeighbor(hemo::Array<plb::plint, 3> offset_) : offset(offset_){};

  virtual void process(plb::Box3D domain, plb::ScalarField3D<int> & field1);

  virtual CopyFromNeighbor * clone() const;

  virtual void getTypeOfModification(std::vector<plb::modif::ModifT> & modified) const;

  virtual plb::BlockDomain::DomainT appliesTo() const;

private:
  hemo::Array<plb::plint, 3> offset;
};

/**
 * @brief Voxelise a closed STL surface into a flag matrix (1 fluid, 0 solid).
 *
 * The geometry is scaled so that it spans @p refDirLength lattice nodes in
 * direction @p refDir (0 x, 1 y, 2 z), and opened at both x ends (inflow and
 * outflow). Collective. The caller owns the created objects.
 * @param[in] meshFileName STL file, relative to the working directory
 * @param[in] extendedEnvelopeWidth envelope of the voxelised blocks [lattice nodes]
 * @param[in] refDirLength extent of the geometry along @p refDir [lattice nodes]
 * @param[in] refDir reference direction (0 x, 1 y, 2 z)
 * @param[out] voxelizedDomain the voxelised domain (its block management)
 * @param[out] flagMatrix 1 at fluid nodes, 0 elsewhere
 * @param[in] blockSize atomic block size [lattice nodes]; <= 0 for one block per rank
 * @param[in] particleEnvelope only used for advice in the log (0: 25)
 */
void getFlagMatrixFromSTL(std::string meshFileName, plb::plint extendedEnvelopeWidth, plb::plint refDirLength,
                          plb::plint refDir, plb::VoxelizedDomain3D<T> *& voxelizedDomain,
                          plb::MultiScalarField3D<int> *& flagMatrix, plint blockSize, int particleEnvelope = 0);
/// As above, with the created objects owned by @p voxelizedDomain and
/// @p flagMatrix (formerly std::auto_ptr, which C++17 removed).
void getFlagMatrixFromSTL(std::string meshFileName, plb::plint extendedEnvelopeWidth, plb::plint refDirLength,
                          plb::plint refDir, std::unique_ptr<plb::VoxelizedDomain3D<T>> & voxelizedDomain,
                          std::unique_ptr<plb::MultiScalarField3D<int>> & flagMatrix, plint blockSize,
                          int particleEnvelope = 0);

/**
 * @brief Voxelise the geometry of a case configuration (see above).
 *
 * Reads from `<domain>`: `geometry` (the STL file, relative to the working
 * directory), `fluidEnvelope`, `refDirN`, `refDir`, `blockSize` and
 * `particleEnvelope`. Replaces the call with these six values that most
 * cases repeat. Collective.
 */
void getFlagMatrixFromSTL(Config & cfg, std::unique_ptr<plb::VoxelizedDomain3D<T>> & voxelizedDomain,
                          std::unique_ptr<plb::MultiScalarField3D<int>> & flagMatrix);
/// As above, with raw pointers owned by the caller (e.g. for a PreInlet,
/// which keeps a pointer to the flag matrix).
void getFlagMatrixFromSTL(Config & cfg, plb::VoxelizedDomain3D<T> *& voxelizedDomain,
                          plb::MultiScalarField3D<int> *& flagMatrix);
}
#endif
