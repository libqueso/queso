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

#ifndef UQ_ARRAY_OF_ONE_D_TABLES_H
#define UQ_ARRAY_OF_ONE_D_TABLES_H

#include <queso/Environment.h>
#include <queso/VectorSpace.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\file ArrayOfOneDTables.h
 * \brief Class to accommodate arrays of one-dimensional tables.
 *
 * \class ArrayOfOneDTables
 * \brief Class to accommodate arrays of one-dimensional tables.
 *
 * Arrays of one-dimensional tables are necessary in the calculation (storage), for
 * instance, of CDFs and MDF of vector functions (refer to BaseVectorCdf,
 * BaseVectorMdf, and derived classes) given the (array of) grid points
 * (ArrayOfOneDGrids).
 */

template<class V = GslVector, class M = GslMatrix>
class ArrayOfOneDTables
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  ArrayOfOneDTables(const char* prefix, const VectorSpace<V,M>& rowSpace);

  //! Destructor.
  ~ArrayOfOneDTables();
  //@}

  //! @name Math methods
  //@{
  //! Sets the one-dimensional table.
  /*! This methods assigns the array \c values to  position \c rowId of the one-dimensional table.*/
  void                       setOneDTable(unsigned int rowId, const std::vector<double>& values);

  //! Returns the array located at position \c rowId of the one-dimensional table.
  const std::vector<double>& oneDTable   (unsigned int rowId) const;
  //@}

  //! @name I/O method
  //@{
  //! Prints the values in this array of tables.
  /*! It prints the arrays (inner for-loop) in each position of the table (outer for-loop).*/
  void                       print       (std::ostream& os)   const;
  //@}
private:
  const BaseEnvironment&          m_env;
        std::string                      m_prefix;
  const VectorSpace<V,M>&         m_rowSpace;
  DistArray<std::vector<double>*> m_oneDTables;
};

}  // End namespace QUESO

#endif // UQ_ARRAY_OF_ONE_D_TABLES_H
