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
#include "readPositionsBloodCells.h"
#include <fstream>
#include <string>
#include <sstream>
#include "logfile.h"
#include "hemocell.h"
#include "preInlet.h"
#include "hemoCellFunctional.h"

#include <map>
#include <set>
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/triangularSurfaceMesh.hh"

using namespace std;
using namespace plb;


namespace hemo {

/// Rotates a triangular surface mesh inplace. The rotation is applied in order
/// X, Y, Z. NOTE: plb::TriangularSurfaceMesh does provide an method to rotate
/// the surface mesh, however, that routine is defined in ZYX order.
template <typename T>
void rotateTriangularMeshXYZ(plb::TriangularSurfaceMesh<T> * mesh, T alpha, T beta, T gamma) {
  // alpha, beta, gamma are arbitrary angles in radians (any sign). The range
  // assertion copied from Palabos' rotate() referred to undefined variables
  // and broke debug (PLB_DEBUG) builds; it does not apply here.

  // Rotation matrix around x axis (column-first)
  T a[3][3];
  a[0][0] = (T)1.0;
  a[0][1] = (T)0.0;
  a[0][2] = (T)0.0;
  a[1][0] = (T)0.0;
  a[1][1] = std::cos(alpha);
  a[1][2] = std::sin(alpha);
  a[2][0] = (T)0.0;
  a[2][1] = -std::sin(alpha);
  a[2][2] = std::cos(alpha);

  // Rotation matrix around y axis (column-first)
  T b[3][3];
  b[0][0] = std::cos(beta);
  b[0][1] = (T)0.0;
  b[0][2] = -std::sin(beta);
  b[1][0] = (T)0.0;
  b[1][1] = (T)1.0;
  b[1][2] = (T)0.0;
  b[2][0] = std::sin(beta);
  b[2][1] = (T)0.0;
  b[2][2] = std::cos(beta);

  // Ry * Rx
  T c[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      c[i][j] = (T)0.0;
      for (int k = 0; k < 3; k++) { c[i][j] += a[k][j] * b[i][k]; }
    }
  }

  // Rotation matrix around z axis (column-first)
  b[0][0] = std::cos(gamma);
  b[0][1] = std::sin(gamma);
  b[0][2] = (T)0.0;
  b[1][0] = -std::sin(gamma);
  b[1][1] = std::cos(gamma);
  b[1][2] = (T)0.0;
  b[2][0] = (T)0.0;
  b[2][1] = (T)0.0;
  b[2][2] = (T)1.0;

  // Rz * [Ry*Rx]
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a[i][j] = (T)0.0;
      for (int k = 0; k < 3; k++) { a[i][j] += c[k][j] * b[i][k]; }
    }
  }

  for (auto iVertex = 0; iVertex < mesh->getNumVertices(); iVertex++) {
    Array<T, 3> x = mesh->getVertex(iVertex);
    for (int i = 0; i < 3; i++) {
      mesh->getVertex(iVertex)[i] = (T)0.0;
      for (int j = 0; j < 3; j++) { mesh->getVertex(iVertex)[i] += a[i][j] * x[j]; }
    }
  }
}

inline void meshRotation(TriangularSurfaceMesh<T> * mesh, hemo::Array<T, 3> rotationAngles) {
  plb::Array<T, 2> xRange, yRange, zRange;
  mesh->computeBoundingBox(xRange, yRange, zRange);
  plb::Array<T, 3> meshCenter =
      plb::Array<T, 3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * (T)0.5;
  mesh->translate((T)-1.0 * meshCenter);
  rotateTriangularMeshXYZ(mesh, rotationAngles[0], rotationAngles[1], rotationAngles[2]);
  mesh->translate(meshCenter);
}

// Whether a vertex is clear of walls as far as this block's lattice can tell:
// its node and the nodes within denyLayerSize of it are not boundary nodes
// (nodes outside the block's lattice are not checked).
inline bool vertexClearOfWalls(BlockLattice3D<T, DESCRIPTOR> & fluid, const hemo::Array<T, 3> & vertex,
                               int denyLayerSize) {
  const Box3D fluidbb = fluid.getBoundingBox();
  //Custom contains function for generic blocklattice3D structure
  Dot3D absfluidloc(int(vertex[0] + 0.5), int(vertex[1] + 0.5), int(vertex[2] + 0.5));
  Dot3D relfluidloc = absfluidloc - fluid.getLocation();

  if (!(relfluidloc.x < fluidbb.x0 || relfluidloc.x > fluidbb.x1 || relfluidloc.y < fluidbb.y0 ||
        relfluidloc.y > fluidbb.y1 || relfluidloc.z < fluidbb.z0 || relfluidloc.z > fluidbb.z1)) {
    // Test when particle is inside if in a boundary -> dont add this particle
    if (fluid.get(relfluidloc.x, relfluidloc.y, relfluidloc.z).getDynamics().isBoundary()) { return false; }
  }

  for (int px = -denyLayerSize; px <= denyLayerSize; ++px) {
    for (int py = -denyLayerSize; py <= denyLayerSize; ++py) {
      for (int pz = -denyLayerSize; pz <= denyLayerSize; ++pz) {
        const Dot3D rrelfluidloc = {relfluidloc.x + px, relfluidloc.y + py, relfluidloc.z + pz};
        if (!(rrelfluidloc.x < fluidbb.x0 || rrelfluidloc.x > fluidbb.x1 || rrelfluidloc.y < fluidbb.y0 ||
              rrelfluidloc.y > fluidbb.y1 || rrelfluidloc.z < fluidbb.z0 || rrelfluidloc.z > fluidbb.z1)) {
          if (fluid.get(rrelfluidloc.x, rrelfluidloc.y, rrelfluidloc.z).getDynamics().isBoundary()) { return false; }
        }
      }
    }
  }
  return true;
}

// Deny particles that are in the outer most layer, aka. the "shear layer"
inline int denyLayerSize(HemoCellParticleField & particleField, pluint celltype) {
  return ((*particleField.cellFields)[celltype]->minimumDistanceFromSolid * 1e-6) / param::dx;
}

inline void positionCellInParticleField(HemoCellParticleField & particleField, BlockLattice3D<T, DESCRIPTOR> & fluid,
                                        TriangularSurfaceMesh<T> * mesh, hemo::Array<T, 3> startingPoint, plint cellId,
                                        pluint celltype) {
  plint nVertices = mesh->getNumVertices();
  const int deny = denyLayerSize(particleField, celltype);

  for (plint iVertex = 0; iVertex < nVertices; ++iVertex) {
    hemo::Array<T, 3> vertex = startingPoint + mesh->getVertex(iVertex);

    //If we cannot place it in the particle field continue
    if (!particleField.isContainedABS(vertex, particleField.getBoundingBox())) { continue; }
    if (!vertexClearOfWalls(fluid, vertex, deny)) { continue; }

    particleField.addParticle(HemoCellParticle(vertex, cellId, iVertex, celltype).sv);
  }
}

int getTotalNumberOfCells(HemoCellFields & cellFields) {
  int nCells, totalCells = 0;
  for (pluint j = 0; j < cellFields.size(); j++) {
    fstream fIn;
    fIn.open(cellFields[j]->name + ".pos", fstream::in);
    fIn >> nCells;
    totalCells += nCells;
    fIn.close();
  }
  return totalCells;
}

namespace {

/// Parse the text of one .pos file: the cell count, then one row per cell with
/// the centre (x y z, micrometres) and the rotation about X, Y and Z (degrees).
/// Angles are converted to radians with the sign flipped (right- to
/// left-handed frame), positions are shifted to the pre-inlet origin when one
/// is active, and cell ids continue from nextCellId.
std::vector<InitialCellPlacement> parsePositions(const std::string & text, const std::string & name,
                                                 HemoCellFields & cellFields, T posRatio, plint & nextCellId) {
  std::istringstream in(text);
  plint declared = 0;
  in >> declared;
  hlog << "(readPositionsBloodCells) Particle count in file (" << name << "): " << declared << "." << endl;

  const bool shiftToPreInlet = cellFields.hemocell.preInlet && cellFields.hemocell.preInlet->initialized;
  std::vector<InitialCellPlacement> cells;
  if (declared > 0) { cells.reserve(declared); }
  for (plint i = 0; i < declared; i++) {
    InitialCellPlacement cell;
    in >> cell.position[0] >> cell.position[1] >> cell.position[2] >> cell.angles[0] >> cell.angles[1] >>
        cell.angles[2];
    if (!in) {
      hlog << "(readPositionsBloodCells) WARNING: " << name << ".pos declares " << declared << " cells but only " << i
           << " complete rows could be read; using those." << endl;
      break;
    }
    for (int d = 0; d < 3; d++) {
      cell.angles[d] *= PI / 180.0; // Deg to Rad
      cell.angles[d] *= -1.0;       // Right- to left-handed coordinate system
    }
    if (shiftToPreInlet) { // Translate system to preInlet Location (set 0,0,0 point)
      cell.position[0] += cellFields.hemocell.preInlet->location.x0 / posRatio;
      cell.position[1] += cellFields.hemocell.preInlet->location.y0 / posRatio;
      cell.position[2] += cellFields.hemocell.preInlet->location.z0 / posRatio;
    }
    cell.cellId = nextCellId++;
    cells.push_back(cell);
  }
  return cells;
}

} // namespace


/* ******** ReadPositionMultipleCellField3D *********************************** */
void ReadPositionsBloodCellField3D::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  int numberOfCellFields = blocks.size() - 1;
  BlockLattice3D<T, DESCRIPTOR> & fluid = *dynamic_cast<BlockLattice3D<T, DESCRIPTOR> *>(blocks[0]);
  std::vector<HemoCellParticleField *> particleFields(numberOfCellFields);
  std::vector<TriangularSurfaceMesh<T> *> meshes(numberOfCellFields);
  std::vector<ElementsOfTriangularSurfaceMesh<T>> emptyEoTSM(numberOfCellFields);

  Dot3D fLocation(fluid.getLocation());
  Box3D realDomain(domain.x0 + fLocation.x, domain.x1 + fLocation.x, domain.y0 + fLocation.y, domain.y1 + fLocation.y,
                   domain.z0 + fLocation.z, domain.z1 + fLocation.z);
  plb::Array<T, 2> xRange, yRange, zRange;

  for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
    TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(cellFields[iCF]->getMesh(), emptyEoTSM[iCF]);
    mesh->computeBoundingBox(xRange, yRange, zRange);
    mesh->translate(plb::Array<T, 3>(-(xRange[0] + xRange[1]) / 2.0, -(yRange[0] + yRange[1]) / 2.0,
                                     -(zRange[0] + zRange[1]) / 2.0));
    meshes[iCF] = mesh;

    particleFields[iCF] = (dynamic_cast<HemoCellParticleField *>(blocks[iCF + 1]));
    particleFields[iCF]->removeParticles(particleFields[iCF]->getBoundingBox());
  }

  // Positions are in micrometres; this converts them to lattice units.
  // Note: this method uses the center of the particles for location
  const T posRatio = 1e-6 / dx;

  // Periodic images. During a run, a cell across a periodic boundary exists
  // as images on both sides: the copy shifted by the domain length L along an
  // axis has the cell ID shifted as by a crossing (+-1, +-offset_y,
  // +-offset_z times the number of cells; HemoCellParticleDataTransfer). The
  // same images are placed here if cellFields.placeCellsAcrossPeriodicBoundaries
  // is set, so that cells crossing a periodic boundary in the .pos file are
  // placed instead of dropped. Then only cells whose centre lies in the domain
  // along every periodic axis are placed (cells centred outside were never
  // placed): a .pos file packed for a larger box would otherwise wrap its
  // outer cells onto the inner ones. Images that lie outside a block's
  // particle field are skipped by the check below. Without it, the loop is
  // the original one (only the cells themselves), and the cells centred in
  // the domain that cross a periodic boundary are counted (`crossing`).
  const bool placeImages = cellFields.placeCellsAcrossPeriodicBoundaries;
  struct PeriodicImage {
    hemo::Array<T, 3> shift; // [lattice units]
    plint idOffset;
    bool original; // the cell itself (no shift)
  };
  std::vector<PeriodicImage> images;
  const Box3D domainBox = cellFields.immersedParticles->getBoundingBox();
  bool periodic[3];
  {
    const plint length[3] = {domainBox.getNx(), domainBox.getNy(), domainBox.getNz()};
    const plint idStep[3] = {1, plint(cellFields.periodicity_limit_offset_y),
                             plint(cellFields.periodicity_limit_offset_z)};
    int range[3];
    for (int d = 0; d < 3; d++) {
      periodic[d] = cellFields.immersedParticles->periodicity().get(d);
      range[d] = periodic[d] ? 1 : 0;
    }
    for (int i = -range[0]; i <= range[0]; i++) {
      for (int j = -range[1]; j <= range[1]; j++) {
        for (int k = -range[2]; k <= range[2]; k++) {
          PeriodicImage image;
          image.shift = {T(i * length[0]), T(j * length[1]), T(k * length[2])};
          image.idOffset = (i * idStep[0] + j * idStep[1] + k * idStep[2]) * plint(cellFields.number_of_cells);
          image.original = (i == 0 && j == 0 && k == 0);
          images.push_back(image);
        }
      }
    }
  }

  for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
    for (const InitialCellPlacement & cell : cells[iCF]) {
      // The centre must lie in the domain along the periodic axes (a node's
      // cell is (x - 0.5, x + 0.5], as for isContainedABS()).
      const hemo::Array<T, 3> centre = cell.position * posRatio;
      const T low[3] = {T(domainBox.x0), T(domainBox.y0), T(domainBox.z0)};
      const T high[3] = {T(domainBox.x1), T(domainBox.y1), T(domainBox.z1)};
      bool home = true;
      for (int d = 0; d < 3; d++) {
        if (periodic[d] && !(centre[d] > low[d] - 0.5 && centre[d] <= high[d] + 0.5)) { home = false; }
      }
      if (placeImages && !home) { continue; }
      const bool count = !placeImages && home && crossing;
      for (const PeriodicImage & image : images) {
        // Without image placement, images are only checked for counting.
        if (!placeImages && !image.original && !count) { continue; }
        //Check if it actually fits (mostly) in this atomic block
        const hemo::Array<T, 3> position = cell.position * posRatio + image.shift;
        if (position[0] < realDomain.x0 - particleEnvelope || position[0] > realDomain.x1 + particleEnvelope ||
            position[1] < realDomain.y0 - particleEnvelope || position[1] > realDomain.y1 + particleEnvelope ||
            position[2] < realDomain.z0 - particleEnvelope || position[2] > realDomain.z1 + particleEnvelope) {
          continue;
        }

        ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
        TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);

        meshRotation(meshCopy, cell.angles);

        if (count) {
          // A crossing cell counts if all its images are clear of the walls,
          // as checked on every block they reach (one wall vetoes it), i.e.
          // if image placement would place it.
          const int deny = denyLayerSize(*particleFields[iCF], iCF);
          bool crosses = false;
          bool clear = true;
          for (plint iVertex = 0; iVertex < meshCopy->getNumVertices() && clear; ++iVertex) {
            const hemo::Array<T, 3> vertex = position + meshCopy->getVertex(iVertex);
            for (int d = 0; d < 3; d++) {
              if (vertex[d] > low[d] - 0.5 && vertex[d] <= high[d] + 0.5) { continue; }
              if (periodic[d]) {
                crosses = true;
              } else {
                clear = false; // outside the domain: not placed either way
              }
            }
            clear = clear && vertexClearOfWalls(fluid, vertex, deny);
          }
          if (!clear) {
            crossing->blocked.insert(cell.cellId);
          } else if (crosses && image.original) {
            crossing->candidates.insert(cell.cellId);
          }
        }

        if (placeImages || image.original) {
          positionCellInParticleField(*(particleFields[iCF]), fluid, meshCopy, position, cell.cellId + image.idOffset,
                                      iCF);
        }
        delete meshCopy;
      }
    }

    particleFields[iCF]->deleteIncompleteCells(iCF, false);
    std::vector<HemoCellParticle *> particles;
    particleFields[iCF]->findParticles(particleFields[iCF]->getBoundingBox(), particles, iCF);

    delete meshes[iCF];
  }
}


ReadPositionsBloodCellField3D * ReadPositionsBloodCellField3D::clone() const {
  return new ReadPositionsBloodCellField3D(*this);
}

void ReadPositionsBloodCellField3D::getTypeOfModification(std::vector<modif::ModifT> & modified) const {
  modified[0] = modif::nothing; // Fluid field.
  for (unsigned int iField = 1; iField < modified.size(); ++iField) {
    modified[iField] = modif::nothing; // Particle fields.
  }
}

void ReadPositionsBloodCellField3D::getModificationPattern(std::vector<bool> & isWritten) const {
  isWritten[0] = false; // Fluid field.
  for (unsigned int iField = 1; iField < isWritten.size(); ++iField) {
    isWritten[iField] = false; // Particle fields.
  }
}


BlockDomain::DomainT ReadPositionsBloodCellField3D::appliesTo() const {
  return BlockDomain::bulk;
}


void readPositionsBloodCellField3D(HemoCellFields & cellFields, T dx, Config & cfg) {
  std::vector<MultiBlock3D *> fluidAndParticleFieldsArg;

  fluidAndParticleFieldsArg.push_back(cellFields.lattice);

  for (pluint icf = 0; icf < cellFields.size(); ++icf) {
    fluidAndParticleFieldsArg.push_back(cellFields[icf]->getParticleField3D());
  }
  hlog << "(readPositionsBloodCells) Reading particle positions..." << std::endl;

  // Rank 0 reads each .pos file and broadcasts its text; every rank parses it
  // once and the functional reuses the result for all of its atomic blocks
  // (it used to re-read and re-parse every file for each block).
  const T posRatio = 1e-6 / dx;
  std::vector<std::vector<InitialCellPlacement>> cells(cellFields.size());
  plint nextCellId = 0;
  for (pluint j = 0; j < cellFields.size(); j++) {
    const std::string name = cellFields[j]->name;
    std::string text;
    int opened = 0;
    if (global::mpi().isMainProcessor()) {
      std::ifstream fIn(name + ".pos");
      if (fIn.is_open()) {
        opened = 1;
        std::ostringstream contents;
        contents << fIn.rdbuf();
        text = contents.str();
      }
    }
    global::mpi().bCast(&opened, 1);
    if (!opened) {
      fatal("(readPositionsBloodCells) The positions file ", name,
            ".pos does not exist (paths are relative to the working directory). Create it, e.g. with "
            "tools/packCells, or write a file containing 0 for no cells of this type.");
    }
    global::mpi().bCast(text);
    cells[j] = parsePositions(text, name, cellFields, posRatio, nextCellId);
    cellFields[j]->cellsInPositionFile = cells[j].size();
  }
  // Cell ids run over all types, so this is the total over all files. Set
  // before placing cells: the periodic cell-id mapping depends on it.
  cellFields.number_of_cells = nextCellId;

  CellsCrossingPeriodicBoundary crossing;
  if (cellFields.hemocell.preInlet && cellFields.hemocell.preInlet->initialized &&
      !cellFields.hemocell.partOfpreInlet) {
  } else {
    const plint particleEnvelope = cfg["domain"]["particleEnvelope"].read<int>();
    applyProcessingFunctional(new ReadPositionsBloodCellField3D(cellFields, dx, cells, particleEnvelope, &crossing),
                              cellFields.lattice->getBoundingBox(), fluidAndParticleFieldsArg);
    hlogfile << "Mpi Process: " << global::mpi().getRank() << " Completed loading particles" << std::endl;
  }
  // Cells skipped across a periodic boundary, counted per type on the main
  // rank (collective: every rank takes part, also those without cells).
  // A cell counts if some block found it clear of the walls and none vetoed it.
  std::map<plint, bool> crossingIds;
  for (const plint id : crossing.candidates) { crossingIds[id] = true; }
  std::map<plint, bool> blocked;
  for (const plint id : crossing.blocked) { blocked[id] = true; }
  HemoCellGatheringFunctional<bool>::gatherToRoot(crossingIds);
  HemoCellGatheringFunctional<bool>::gatherToRoot(blocked);
  for (const auto & entry : blocked) { crossingIds.erase(entry.first); }
  for (pluint j = 0; j < cellFields.size(); j++) {
    cellFields[j]->cellsCrossingPeriodicBoundary = 0;
    for (const InitialCellPlacement & cell : cells[j]) {
      if (crossingIds.count(cell.cellId)) { cellFields[j]->cellsCrossingPeriodicBoundary++; }
    }
  }
}

}
