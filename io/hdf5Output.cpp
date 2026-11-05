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
#include "hdf5Output.h"

namespace hemo {

void writeCompressedDataset(hid_t file, const std::string & name, hid_t type, int rank, const hsize_t * dims,
                            const hsize_t * chunk, const void * data) {
  hid_t space = H5Screate_simple(rank, dims, NULL);
  hid_t properties = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(properties, rank, chunk);
  H5Pset_shuffle(properties);
  H5Pset_deflate(properties, outputDeflateLevel);
  hid_t dataset = H5Dcreate2(file, name.c_str(), type, space, H5P_DEFAULT, properties, H5P_DEFAULT);
  H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  H5Pclose(properties);
  H5Sclose(space);
}

} // namespace hemo
