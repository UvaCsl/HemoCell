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
#include "ParticleHdf5IO.h"
#include "hdf5Output.h"
#include "hemocell.h"

#include <hdf5.h>
#include <hdf5_hl.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace plb;

namespace hemo {

namespace {
// Write rows x columns values of a 2D output as a dataset; returns the number of rows.
template <typename Stored, typename Value>
hsize_t writeRows(hid_t file, const std::string & name, hid_t type, const vector<vector<Value>> & rows) {
  hsize_t dims[2];
  dims[0] = rows.size();
  dims[1] = rows.empty() ? 0 : rows[0].size();
  hsize_t chunk[2];
  chunk[0] = std::max<hsize_t>(1, std::min<hsize_t>(1000, dims[0]));
  chunk[1] = std::max<hsize_t>(1, dims[1]);
  std::vector<Stored> data;
  data.reserve(dims[0] * dims[1]);
  for (const vector<Value> & row : rows) {
    for (const Value & v : row) { data.push_back(Stored(v)); }
  }
  writeCompressedDataset(file, name, type, 2, dims, chunk, data.data());
  return dims[0];
}
} // namespace

/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
WriteCellField3DInMultipleHDF5Files::WriteCellField3DInMultipleHDF5Files(HemoCellField & cellField3D_, plint iter_,
                                                                         std::string identifier_, T dx_, T dt_,
                                                                         std::vector<XdmfBlockFile> * written_)
    : cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_), written(written_){};

void WriteCellField3DInMultipleHDF5Files::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  PLB_PRECONDITION(blocks.size() > 0);
  int id = global::mpi().getRank();
  long int size = global::mpi().getSize();

  HemoCellParticleField & particleField = *dynamic_cast<HemoCellParticleField *>(blocks[0]);

  if (cellField3D.desiredOutputVariables.size() == 0) {
    return; //No output desired, no problem
  }

  XdmfBlockFile record;
  record.identifier = identifier;
  if (cellField3D.cellFields.hemocell.partOfpreInlet) { record.identifier += "_PRE"; }
  record.blockId = particleField.atomicBlockId;
  record.file = "hdf5/" + zeroPadNumber(iter) + '/' + record.identifier + "." + zeroPadNumber(iter) + ".p." +
                to_string(particleField.atomicBlockId) + ".h5";
  /************************************************************/
  /**            Initialise HDF5 file                        **/
  /************************************************************/
  std::string fileName = global::directories().getOutputDir() + "/" + record.file;
  hid_t file_id;
  file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  H5LTset_attribute_double(file_id, "/", "dx", &dx, 1);
  H5LTset_attribute_double(file_id, "/", "dt", &dt, 1);
  long int iterHDF5 = iter;
  H5LTset_attribute_long(file_id, "/", "iteration", &iterHDF5, 1);
  H5LTset_attribute_int(file_id, "/", "processorId", &id, 1);
  H5LTset_attribute_long(file_id, "/", "numberOfProcessors", &size, 1);

  /************************************************************/
  /**            Write output to HDF5 file                   **/
  /************************************************************/
  bool hasPosition = false;
  vector<vector<T>> output;
  for (pluint i = 0; i < cellField3D.desiredOutputVariables.size(); i++) {
    std::string vectorname = "";
    particleField.passthroughpass(cellField3D.desiredOutputVariables[i], domain, output, cellField3D.ctype, vectorname);
    if (vectorname == "") { continue; }
    const hsize_t rows = writeRows<float>(file_id, vectorname, H5T_NATIVE_FLOAT, output);
    // Blocks without vertices are left out of the XDMF file, so the
    // components of an empty output do not matter.
    XdmfAttribute attribute;
    attribute.name = vectorname;
    attribute.components = output.empty() ? 0 : int(output[0].size());
    record.attributes.push_back(attribute);

    if (cellField3D.desiredOutputVariables[i] == OUTPUT_POSITION) {
      hasPosition = true;
      record.vertices = rows;
      long int nP = rows;
      H5LTset_attribute_long(file_id, "/", "numberOfParticles", &nP, 1);
    }
  }

  if (cellField3D.outputTriangles) {
    vector<vector<plint>> triangles;
    std::string vectorname;
    particleField.outputTriangles(domain, triangles, cellField3D.ctype, vectorname);
    long int nT = writeRows<int>(file_id, vectorname, H5T_NATIVE_INT, triangles);
    record.triangles = nT;
    H5LTset_attribute_long(file_id, "/", "numberOfTriangles", &nT, 1);
  }

  if (std::find(cellField3D.desiredOutputVariables.begin(), cellField3D.desiredOutputVariables.end(),
                OUTPUT_INNER_LINKS) != cellField3D.desiredOutputVariables.end()) {
    vector<vector<plint>> links;
    std::string vectorname;
    particleField.outputInnerLinks(domain, links, cellField3D.ctype, vectorname);
    if (!links.empty()) {
      long int nT = writeRows<int>(file_id, vectorname, H5T_NATIVE_INT, links);
      record.innerLinks = nT;
      H5LTset_attribute_long(file_id, "/", "numberOfInnerLinks", &nT, 1);
    }
  }

  H5Fclose(file_id);
  // Without positions there is no geometry to describe.
  if (hasPosition) { written->push_back(record); }
}

WriteCellField3DInMultipleHDF5Files * WriteCellField3DInMultipleHDF5Files::clone() const {
  return new WriteCellField3DInMultipleHDF5Files(*this);
}

void WriteCellField3DInMultipleHDF5Files::getTypeOfModification(std::vector<modif::ModifT> & modified) const {
  for (pluint i = 0; i < modified.size(); ++i) { modified[i] = modif::nothing; }
}

BlockDomain::DomainT WriteCellField3DInMultipleHDF5Files::appliesTo() const {
  return BlockDomain::bulk;
}


void writeCellField3D_HDF5(HemoCellFields & cellFields, T dx, T dt, plint iter, std::string preString) {
  global.statistics.getCurrent()["writeCellField"].start();

  std::vector<XdmfBlockFile> written;
  for (pluint i = 0; i < cellFields.size(); i++) {
    std::string identifier = preString + cellFields[i]->getIdentifier();
    WriteCellField3DInMultipleHDF5Files * bprf =
        new WriteCellField3DInMultipleHDF5Files(*cellFields[i], iter, identifier, dx, dt, &written);
    vector<MultiBlock3D *> wrapper;
    wrapper.push_back(cellFields[i]->getParticleArg());
    applyProcessingFunctional(bprf, cellFields[i]->getParticleField3D()->getBoundingBox(), wrapper);
  }
  writeXdmf(written, iter, cellFields.hemocell.outputInSiUnits ? iter * dt : double(iter));

  global.statistics.getCurrent().stop();
}

}
