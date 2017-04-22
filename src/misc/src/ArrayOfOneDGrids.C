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

#include <queso/OneDGrid.h>
#include <queso/UniformOneDGrid.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -------------------------------------------------
template <class V, class M>
ArrayOfOneDGrids<V,M>::ArrayOfOneDGrids(
  const char*                    prefix,
  const VectorSpace<V,M>& rowSpace)
  :
  m_env         (rowSpace.env()    ),
  m_prefix      ((std::string)(prefix)+""),
  m_rowSpace    (rowSpace          ),
  m_oneDGrids   (m_rowSpace.map(),1),
  m_sizes       (NULL),
  m_minPositions(NULL),
  m_maxPositions(NULL)
{
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    m_oneDGrids(i,0) = NULL;
  }
}

// Destructor ----------------------------------------------------------
template <class V, class M>
ArrayOfOneDGrids<V,M>::~ArrayOfOneDGrids()
{
  if (m_maxPositions) delete m_maxPositions;
  if (m_minPositions) delete m_minPositions;
  if (m_sizes       ) delete m_sizes;

  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    if (m_oneDGrids(i,0)) delete m_oneDGrids(i,0);
  }
}

// Property methods-----------------------------------------------------
template <class V, class M>
const VectorSpace<V,M>&
ArrayOfOneDGrids<V,M>::rowSpace() const
{
  return m_rowSpace;
}

template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::sizes() const
{
  queso_require_msg(m_sizes, "sizes is still NULL");

  return *m_sizes;
}

template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::minPositions() const
{
  queso_require_msg(m_minPositions, "minPositions is still NULL");

  return *m_minPositions;
}

template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::maxPositions() const
{
  queso_require_msg(m_maxPositions, "maxPositions is still NULL");

  return *m_maxPositions;
}

template <class V, class M>
void
ArrayOfOneDGrids<V,M>::setUniformGrids(
  const V& sizesVec,
  const V& minPositionsVec,
  const V& maxPositionsVec)
{
  if (m_sizes        == NULL) m_sizes        = new V(sizesVec);
  else                       *m_sizes        = sizesVec;

  if (m_minPositions == NULL) m_minPositions = new V(minPositionsVec);
  else                       *m_minPositions = minPositionsVec;

  if (m_maxPositions == NULL) m_maxPositions = new V(maxPositionsVec);
  else                       *m_maxPositions = maxPositionsVec;

  char strI[65];
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    sprintf(strI,"%u_",i);
    m_oneDGrids(i,0) = new UniformOneDGrid<double>(m_env,
                                                          (m_prefix+strI).c_str(),
                                                          (unsigned int) sizesVec[i],
                                                          minPositionsVec[i],
                                                          maxPositionsVec[i]);
  }

  return;
}

template <class V, class M>
const BaseOneDGrid<double>&
ArrayOfOneDGrids<V,M>::grid(unsigned int rowId) const
{
  queso_require_less_msg(rowId, m_rowSpace.dimLocal(), "rowId is out of range");

  ArrayOfOneDGrids<V,M>* tmp = const_cast<ArrayOfOneDGrids<V,M>*>(this);
  return *(tmp->m_oneDGrids(rowId,0));
}

// I/O methods----------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDGrids<V,M>::print(std::ostream& os) const
{
  ArrayOfOneDGrids<V,M>* tmp = const_cast<ArrayOfOneDGrids<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    os << *(tmp->m_oneDGrids(i,0))
       << std::endl;
  }

  return;
}

}  // End namespace QUESO

template class QUESO::ArrayOfOneDGrids<QUESO::GslVector, QUESO::GslMatrix>;
