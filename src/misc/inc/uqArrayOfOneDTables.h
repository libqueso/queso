//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_ARRAY_OF_ONE_D_TABLES_H__
#define __UQ_ARRAY_OF_ONE_D_TABLES_H__

#include <uqEnvironment.h>

namespace QUESO {

/*!\file uqArrayOfOneDTables
 * \brief Class to accommodate arrays of one-dimensional tables.
 * 
 * \class ArrayOfOneDTablesClass
 * \brief Class to accommodate arrays of one-dimensional tables.
 *  
 * Arrays of one-dimensional tables are necessary in the calculation (storage), for 
 * instance, of CDFs and MDF of vector functions (refer to BaseVectorCdfClass, 
 * BaseVectorMdfClass, and derived classes) given the (array of) grid points 
 * (ArrayOfOneDGridsClass).
 */

template <class V, class M>
class ArrayOfOneDTablesClass
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  ArrayOfOneDTablesClass(const char* prefix, const VectorSpaceClass<V,M>& rowSpace);

  //! Destructor.
  ~ArrayOfOneDTablesClass();
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
  const BaseEnvironmentClass&          m_env;
        std::string                      m_prefix;
  const VectorSpaceClass<V,M>&         m_rowSpace;
  DistArrayClass<std::vector<double>*> m_oneDTables;
};
// Default constructor -------------------------------------------------
template <class V, class M>
ArrayOfOneDTablesClass<V,M>::ArrayOfOneDTablesClass(
  const char* prefix, 
  const VectorSpaceClass<V,M>& rowSpace)
  :
  m_env       (rowSpace.env()),
  m_prefix    ((std::string)(prefix)+""),
  m_rowSpace  (rowSpace      ),
  m_oneDTables(m_rowSpace.map(),1)
{
  for (unsigned int i = 0; i < (unsigned int) m_oneDTables.MyLength(); ++i) {
    m_oneDTables(i,0) = NULL;
  }
}
// Destructor ----------------------------------------------------------
template <class V, class M>
ArrayOfOneDTablesClass<V,M>::~ArrayOfOneDTablesClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_oneDTables.MyLength(); ++i) {
    if (m_oneDTables(i,0)) delete m_oneDTables(i,0);
  }
}
// Math methods --------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDTablesClass<V,M>::setOneDTable(unsigned int rowId, const std::vector<double>& values)
{
  UQ_FATAL_TEST_MACRO(rowId >= (unsigned int) m_oneDTables.MyLength(),
                      m_env.worldRank(),
                      "ArrayOfOneDTablesClass<T>::setOneDTable()",
                      "rowId is out of range");

  if (m_oneDTables(rowId,0) == NULL) m_oneDTables(rowId,0) = new std::vector<double>(0);
  else                               m_oneDTables(rowId,0)->clear();

  std::vector<double>& vec = *(m_oneDTables(rowId,0));
  vec.resize(values.size(),0.);
  for (unsigned int j = 0; j < values.size(); ++j) {
    vec[j] = values[j];
  }

  return;
}
//----------------------------------------------------------------------
template <class V, class M>
const std::vector<double>&
ArrayOfOneDTablesClass<V,M>::oneDTable(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= (unsigned int) m_oneDTables.MyLength(),
                      m_env.worldRank(),
                      "ArrayOfOneDTablesClass<T>::oneDTable()",
                      "rowId is out of range");

  ArrayOfOneDTablesClass<V,M>* tmp = const_cast<ArrayOfOneDTablesClass<V,M>*>(this);

  UQ_FATAL_TEST_MACRO(tmp->m_oneDTables(rowId,0) == NULL,
                      m_env.worldRank(),
                      "ArrayOfOneDTablesClass<T>::oneDTable()",
                      "requested row is still NULL");

  return *(tmp->m_oneDTables(rowId,0));
}
// I/O methods----------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDTablesClass<V,M>::print(std::ostream& os) const
{
  ArrayOfOneDTablesClass<V,M>* tmp = const_cast<ArrayOfOneDTablesClass<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_oneDTables.MyLength(); ++i) {
    const std::vector<double>& tmpVec = *(tmp->m_oneDTables(i,0));
    os << m_prefix << i << "_values_sub" << m_env.subIdString() << " = zeros(" << tmpVec.size()
       << ","                                                                  << 1
       << ");"
       << std::endl;
    os << m_prefix << i << "_values_sub" << m_env.subIdString() << " = [";
    for (unsigned int j = 0; j < tmpVec.size(); ++j) {
      os << tmpVec[j] << " ";
    }
    os << "];"
       << std::endl;
  }

  return;
}
//----------------------------------------------------------------------
template <class V, class M>
std::ostream& operator<< (std::ostream& os, const ArrayOfOneDTablesClass<V,M>& obj)
{
  obj.print(os);
  return os;
}

}  // End namespace QUESO

#endif // __UQ_ARRAY_OF_ONE_D_TABLES_H__
