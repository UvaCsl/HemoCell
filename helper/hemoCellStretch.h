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
#ifndef HEMOCELLSTRETCH_H
#define HEMOCELLSTRETCH_H

#include "hemoCellFields.h"
#include "hemoCellFunctional.h"

namespace hemo {

/**
 * @brief Stretch a single cell by pulling vertices at both ends along x, as
 * in an optical-tweezers experiment.
 *
 * The data are per object, so several simulations (or stretches) in one
 * process do not influence each other. The object must live as long as
 * applyForce() is called.
 */
class HemoCellStretch {
public:
  /**
   * @brief Select the vertices to pull.
   *
   * Picks the @p n_forced_lsps_ vertices of the cell with the lowest and the
   * highest x position. Synchronises the particle envelopes; all MPI ranks
   * must call it.
   * @param[in] cellfield_ Cell type with exactly one cell (cell ID 0).
   * @param[in] n_forced_lsps_ Number of vertices pulled at each end.
   * @param[in] external_force_ Total force per end [lattice force units],
   *            divided equally over its vertices.
   */
  HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, T external_force_);
  /// Add the pulling force (times #scale) to the vertices; call after every
  /// iterate() (needs a material update interval of 1). All MPI ranks must call it.
  void applyForce();

  class FindForcedLsps : public HemoCellFunctional {
  public:
    explicit FindForcedLsps(HemoCellStretch & stretch_) : stretch(stretch_) {}

  private:
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>) override;
    FindForcedLsps * clone() const override;
    HemoCellStretch & stretch;
  };

  class ForceForcedLsps : public HemoCellFunctional {
  public:
    explicit ForceForcedLsps(const HemoCellStretch & stretch_) : stretch(stretch_) {}

  private:
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D *>) override;
    ForceForcedLsps * clone() const override;
    const HemoCellStretch & stretch;
  };

  /// Vertex IDs of the vertices pulled towards -x (lower) and +x (upper).
  std::vector<plint> lower_lsps;
  std::vector<plint> upper_lsps;

  HemoCellField & cellfield;
  unsigned int n_forced_lsps;
  /// Force per pulled vertex [lattice force units].
  T external_force;
  /// Factor applied to the force at the next applyForce() calls (default 1).
  T scale = 1.0;
};
}
#endif
