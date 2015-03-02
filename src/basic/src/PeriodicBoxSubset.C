//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/VectorSpace.h>
#include <queso/PeriodicBoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default, shaped constructor
template<class V, class M>
PeriodicBoxSubset<V,M>::PeriodicBoxSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    const V& minValues,
    const V& maxValues)
  : VectorSubset<V,M>(prefix,vectorSpace,0.),
    m_minValues(minValues),
    m_maxValues(maxValues)
{
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != maxValues.sizeLocal(),
                      m_env.worldRank(),
                      "PeriodicBoxSubset<V,M>::PeriodicBoxSubset()",
                      "vectors 'minValues' and 'maxValues' should have the same size");
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != vectorSpace.dimLocal(),
                      m_env.worldRank(),
                      "PeriodicBoxSubset<V,M>::PeriodicBoxSubset()",
                      "sizes of vectors 'minValues' and 'maxValues' should be equal to dimension of the vector space");
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    UQ_FATAL_TEST_MACRO(minValues[i] > maxValues[i],
                        m_env.worldRank(),
                        "PeriodicBoxSubset<V,M>::PeriodicBoxSubset()",
                        "it should happen minValue <= maxValue for all dimensions");
  }

  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    m_volume *= (m_maxValues[i] - m_minValues[i]);
    m_periods[i] = m_maxValues[i] - m_minValues[i];
  }
}

// Destructor
template<class V, class M>
PeriodicBoxSubset<V,M>::~PeriodicBoxSubset()
{
}

// Math methods
template<class V, class M>
bool PeriodicBoxSubset<V,M>::contains(const V& vec) const
{
  return (!vec.atLeastOneComponentSmallerThan(m_minValues) &&
          !vec.atLeastOneComponentBiggerThan (m_maxValues));
}


template<class V, class M>
const V& PeriodicBoxSubset<V,M>::minValues() const
{
  return m_minValues;
}

template<class V, class M>
const V& PeriodicBoxSubset<V,M>::maxValues() const
{
  return m_maxValues;
}


template<class V, class M>
const V& PeriodicBoxSubset<V,M>::periods() const
{
  return m_periods;
}

// I/O method
template <class V, class M>
void PeriodicBoxSubset<V,M>::print(std::ostream& os) const
{
  os << "In PeriodicBoxSubset<V,M>::print()"
     << ": m_minValues = " << m_minValues
     << ", m_maxValues = " << m_maxValues
     << ", m_periods= "    << m_periods
     << ", m_volume = "    << m_volume
     << std::endl;

  return;
}

}  // End namespace QUESO

template class QUESO::PeriodicBoxSubset<QUESO::GslVector, QUESO::GslMatrix>;
