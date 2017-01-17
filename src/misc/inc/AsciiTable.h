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

#ifndef UQ_ASCII_TABLE_H
#define UQ_ASCII_TABLE_H

#include <queso/Environment.h>
#include <queso/DistArray.h>
#include <queso/Map.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\file AsciiTable.h
 * \brief Class to read ASCII values from a table.
 *
 * \class AsciiTable
 * \brief Class for reading ASCII values from a table in a file.*/

template<class V = GslVector, class M = GslMatrix>
class AsciiTable
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! This constructor reads the data from file \c fileName, checking whether the data in each
   * column of the file is or not a string, and whether the data in each row is or not valid.*/
  AsciiTable(const BaseEnvironment& env,
                          unsigned int            numRows,
                          unsigned int            numExtraCols,
                    const std::vector<bool>*      extraColIsString,
                    const std::string&            fileName);
  //! Destructor.
  ~AsciiTable();
  //@}

  //! @name Property methods
  //@{
  //! Returns the number of rows in the table.
  unsigned int                         numRows     ()                 const;

  //! Returns the number of columns in the table.
  unsigned int                         numCols     ()                 const;

  //! Returns the string stored in column \c j.
  const DistArray<std::string>& stringColumn(unsigned int j)   const;

  //! Returns the value (double) stored in column \c j.
  const V&                             doubleColumn(unsigned int j)   const;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the table.
  void                                 print       (std::ostream& os) const;
  //@}
private:
  Map* newMap(); // See template specialization

  const BaseEnvironment&               m_env;
  unsigned int                                m_numRows;
  unsigned int                                m_numCols;
  std::vector<bool>                           m_colIsString;
  std::string                                 m_fileName;

  const Map*                           m_map;
  std::vector<DistArray<std::string>*> m_stringColumns;
  std::vector<V*>                             m_doubleColumns;

  void readColumnsFromFile();
};

}  // End namespace QUESO

#endif // UQ_ASCII_TABLE_H
