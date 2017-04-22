//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef QUESO_FILE_PTR
#define QUESO_FILE_PTR

#include <queso/Defines.h>

#ifdef QUESO_HAS_HDF5
#include <hdf5.h>
#endif

namespace QUESO {

/*! \struct FilePtrSetStruct
 *  \brief Struct for handling data input and output from files.
 *
 *  This struct deals with data input and output from files.
 *  It encapsulates the input/output stream class std:: fstream.
 */
struct FilePtrSetStruct {
  //! Struct constructor
  FilePtrSetStruct();

  //! Destructor
  ~FilePtrSetStruct();

  //! Provides a stream interface to write data to files.
  std::ofstream* ofsVar;

  //! Provides a stream interface to read data from files.
  std::ifstream* ifsVar;
#ifdef QUESO_HAS_HDF5
  hid_t  h5Var;
#endif
};

}  // End namespace QUESO

#endif  // QUESO_FILE_PTR
