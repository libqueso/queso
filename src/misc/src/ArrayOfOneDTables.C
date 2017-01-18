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

#include <queso/ArrayOfOneDTables.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -------------------------------------------------
template <class V, class M>
ArrayOfOneDTables<V,M>::ArrayOfOneDTables(
  const char* prefix,
  const VectorSpace<V,M>& rowSpace)
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
ArrayOfOneDTables<V,M>::~ArrayOfOneDTables()
{
  for (unsigned int i = 0; i < (unsigned int) m_oneDTables.MyLength(); ++i) {
    if (m_oneDTables(i,0)) delete m_oneDTables(i,0);
  }
}

// Math methods --------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDTables<V,M>::setOneDTable(unsigned int rowId, const std::vector<double>& values)
{
  queso_require_less_msg(rowId, (unsigned int) m_oneDTables.MyLength(), "rowId is out of range");

  if (m_oneDTables(rowId,0) == NULL) m_oneDTables(rowId,0) = new std::vector<double>(0);
  else                               m_oneDTables(rowId,0)->clear();

  std::vector<double>& vec = *(m_oneDTables(rowId,0));
  vec.resize(values.size(),0.);
  for (unsigned int j = 0; j < values.size(); ++j) {
    vec[j] = values[j];
  }

  return;
}

template <class V, class M>
const std::vector<double>&
ArrayOfOneDTables<V,M>::oneDTable(unsigned int rowId) const
{
  queso_require_less_msg(rowId, (unsigned int) m_oneDTables.MyLength(), "rowId is out of range");

  ArrayOfOneDTables<V,M>* tmp = const_cast<ArrayOfOneDTables<V,M>*>(this);

  queso_require_msg(tmp->m_oneDTables(rowId,0), "requested row is still NULL");

  return *(tmp->m_oneDTables(rowId,0));
}

// I/O methods----------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDTables<V,M>::print(std::ostream& os) const
{
  ArrayOfOneDTables<V,M>* tmp = const_cast<ArrayOfOneDTables<V,M>*>(this);
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

template <class V, class M>
std::ostream& operator<< (std::ostream& os, const ArrayOfOneDTables<V,M>& obj)
{
  obj.print(os);
  return os;
}

}  // End namespace QUESO

template class QUESO::ArrayOfOneDTables<QUESO::GslVector, QUESO::GslMatrix>;
