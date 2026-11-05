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
#ifndef HEMO_HDF5_OUTPUT_H
#define HEMO_HDF5_OUTPUT_H

#include <hdf5.h>

#include <string>

namespace hemo {

/// Deflate (zlib) level of every HDF5 output dataset. The bytes are shuffled
/// before deflating, which compresses floating-point data better and faster
/// than a higher level without shuffling.
constexpr unsigned int outputDeflateLevel = 4;

/**
 * @brief Create the dataset @p name in @p file and write @p data to it,
 * chunked, byte-shuffled and deflate-compressed.
 *
 * @param[in] file Open HDF5 file.
 * @param[in] name Dataset name.
 * @param[in] type HDF5 native type of @p data, also the stored type
 *                 (e.g. `H5T_NATIVE_FLOAT`, `H5T_NATIVE_INT`).
 * @param[in] rank Number of dimensions.
 * @param[in] dims Size of each dimension, slowest varying first; may be 0.
 * @param[in] chunk Chunk size of each dimension, each at least 1.
 * @param[in] data Array of `dims[0] * ... * dims[rank-1]` values in C order.
 */
void writeCompressedDataset(hid_t file, const std::string & name, hid_t type, int rank, const hsize_t * dims,
                            const hsize_t * chunk, const void * data);

} // namespace hemo
#endif
