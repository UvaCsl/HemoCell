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

The library is distre from all 3 mirrors simultaneously.
Reply
Collapse
ibuted in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "constant_defaults.h"
#include "voxelizeDomain.h"
#include "logfile.h"

#include "palabos3D.h"
#include "palabos3D.hh"

namespace hemo {
  using namespace std;
  using namespace plb;
// ----------------------- Copy from neighbour ------------------------------------
void CopyFromNeighbor::process(
        Box3D domain, ScalarField3D<int> &field1) {
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field1.get(iX, iY, iZ) = field1.get(iX + offset[0], iY + offset[1], iZ + offset[2]);
            }
        }
    }
}

CopyFromNeighbor *CopyFromNeighbor::clone() const {
    return new CopyFromNeighbor(*this);
}

void CopyFromNeighbor::getTypeOfModification(std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::allVariables;
}

BlockDomain::DomainT CopyFromNeighbor::appliesTo() const {
    return BlockDomain::bulk;
}


// ---------------------- Read in STL geometry ---------------------------------

void getFlagMatrixFromSTL(std::string meshFileName, plint extendedEnvelopeWidth, plint refDirLength, plint refDir,
                          VoxelizedDomain3D<T> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix, plint blockSize, int particleEnvelope) {
    plint extraLayer = 0;   // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
    plint borderWidth = 1;  // Because the Guo boundary condition acts in a one-cell layer.
    
    // Requirement: margin>=borderWidth.
    plint margin = 1;  // Extra margin of allocated cells around the obstacle.

    TriangleSet<T> *triangleSet = new TriangleSet<T>(meshFileName, DBL);

    DEFscaledMesh<T> *defMesh =
            new DEFscaledMesh<T>(*triangleSet, refDirLength, refDir, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    voxelizedDomain = new VoxelizedDomain3D<T>(
            boundary, voxelFlag::inside, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
    
    // Print out some info
    hlog << "(main) Voxelisation is done. Resulting domain parameters are: " << endl;
    hlog << getMultiBlockInfo(voxelizedDomain->getVoxelMatrix()) << std::endl;
	
    // Use particle envelope size if available
    if(!particleEnvelope) {
      hlog << "(Voxelizer) (Warning) particleEnvelopeSize not given, setting to 25" << endl;
      particleEnvelope = 25;
    }
    // Print out if parameters where optimal.
    if (blockSize > 0) {
      if (blockSize < particleEnvelope+1 || blockSize > particleEnvelope*1.25) {
        hlog << "(Voxelizer) (Warning) BlockSize is non optimal (" << blockSize << "), consider setting it to [" << to_string(particleEnvelope+1) << "," << to_string(particleEnvelope*1.25) << "]" << endl;
      }
      plint numBlocks = voxelizedDomain->getVoxelMatrix().getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
      if (numBlocks < global::mpi().getSize()) {
        hlog << "(Voxelizer) (Warning) There are fewer blocks ("<< numBlocks <<") than processors ("<<global::mpi().getSize() <<"), you are wasting CPUs!" << endl;
      } else if (numBlocks > global::mpi().getSize()) {
        hlog << "(Voxelizer) (Warning) There are more blocks (" << numBlocks <<") than processors ("<<global::mpi().getSize() <<"), you might want to consider running with " << numBlocks << " CPUs" << endl;
      }
    } else {
      VoxelizedDomain3D<T> * temp = new VoxelizedDomain3D<T>(
            boundary, voxelFlag::inside, extraLayer, borderWidth, extendedEnvelopeWidth, particleEnvelope+3);
      plint numBlocks = temp->getVoxelMatrix().getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
      if (numBlocks != global::mpi().getSize()) {
        hlog << "(Voxelizer) (Warning) Running with " << global::mpi().getSize() << " CPU's, consider running with " << numBlocks << " CPU's for optimal efficiency" << endl;
      }
      delete temp;
    }
    
    flagMatrix = new MultiScalarField3D<int>((MultiBlock3D &) voxelizedDomain->getVoxelMatrix());

    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix->getBoundingBox(), 1);
    setToConstant(*flagMatrix, voxelizedDomain->getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix->getBoundingBox(), 1);


	// Since the domain is closed, open up the two ends by copying the slice before it.
    Box3D domainBox = flagMatrix->getBoundingBox();
    plint nx = domainBox.getNx();
    plint ny = domainBox.getNy();
    plint nz = domainBox.getNz();

    Box3D domain(0, 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor(hemo::Array<plint, 3>({1, 0, 0})), domain, *flagMatrix);

    domain = Box3D(nx - 2, nx - 1, 0, ny - 1, 0, nz - 1);
    applyProcessingFunctional(new CopyFromNeighbor(hemo::Array<plint, 3>({-1, 0, 0})), domain, *flagMatrix);

}

}
