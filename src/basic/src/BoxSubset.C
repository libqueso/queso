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

#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default, shaped constructor
template<class V, class M>
BoxSubset<V,M>::BoxSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    const V& minValues,
    const V& maxValues)
  : VectorSubset<V,M>(prefix,vectorSpace,0.),
    m_minValues(minValues),
    m_maxValues(maxValues)
{
  queso_require_equal_to_msg(minValues.sizeLocal(), maxValues.sizeLocal(), "vectors 'minValues' and 'maxValues' should have the same size");
  queso_require_equal_to_msg(minValues.sizeLocal(), vectorSpace.dimLocal(), "sizes of vectors 'minValues' and 'maxValues' should be equal to dimension of the vector space");
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    queso_require_less_equal_msg(minValues[i], maxValues[i], "it should happen minValue <= maxValue for all dimensions");
  }

  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    m_volume *= (m_maxValues[i] - m_minValues[i]);
  }
}

// Destructor
template<class V, class M>
BoxSubset<V,M>::~BoxSubset()
{
}

// Math methods
template<class V, class M>
bool BoxSubset<V,M>::contains(const V& vec) const
{
  // prudenci, 2012-09-26: allow boundary values because of 'beta' realizer, which can generate a sample with boundary value '1'
  //return (!vec.atLeastOneComponentSmallerOrEqualThan(m_minValues) &&
  //        !vec.atLeastOneComponentBiggerOrEqualThan (m_maxValues));
  return (!vec.atLeastOneComponentSmallerThan(m_minValues) &&
          !vec.atLeastOneComponentBiggerThan (m_maxValues));
}


template<class V, class M>
void BoxSubset<V,M>::centroid(V& vec) const
{
  vec = m_minValues;
  vec += m_maxValues;
  vec *= 0.5;
}



template<class V, class M>
void BoxSubset<V,M>::moments(M& mat) const
{
  mat.zeroLower();
  mat.zeroUpper();
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    double length_i = (m_maxValues[i] - m_minValues[i]);
    mat(i,i) = length_i*length_i*length_i/12;
  }
}


template<class V, class M>
const V& BoxSubset<V,M>::minValues() const
{
  return m_minValues;
}

template<class V, class M>
const V& BoxSubset<V,M>::maxValues() const
{
  return m_maxValues;
}

// I/O method
template <class V, class M>
void BoxSubset<V,M>::print(std::ostream& os) const
{
  os << "In BoxSubset<V,M>::print()"
     << ": m_minValues = " << m_minValues
     << ", m_maxValues = " << m_maxValues
     << ", m_volume = "    << m_volume
     << std::endl;

  return;
}

}  // End namespace QUESO

template class QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>;
