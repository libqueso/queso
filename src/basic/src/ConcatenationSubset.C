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

#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/ConcatenationSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default, shaped constructor
template<class V, class M>
ConcatenationSubset<V,M>::ConcatenationSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    const VectorSet<V,M>& set1,
    const VectorSet<V,M>& set2)
  : VectorSubset<V,M>(prefix, vectorSpace, set1.volume()*set2.volume()),
    m_sets(2, (const VectorSet<V,M>*) NULL)
{
  m_sets[0] = &set1;
  m_sets[1] = &set2;
}

// Default, shaped constructor
template<class V, class M>
ConcatenationSubset<V,M>::ConcatenationSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    double volume,
    const std::vector<const VectorSet<V,M>* >& sets)
  : VectorSubset<V,M>(prefix, vectorSpace, volume),
    m_sets(sets.size(), (const VectorSet<V,M>*) NULL)
{
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    m_sets[i] = sets[i];
  }
}

// Destructor
template<class V, class M>
ConcatenationSubset<V,M>::~ConcatenationSubset()
{
}

// Mathematical methods
template<class V, class M>
bool ConcatenationSubset<V,M>::contains(const V& vec) const
{
  bool result = true;

  std::vector<V*> vecs(m_sets.size(),(V*) NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_sets[i]->vectorSpace().zeroVector());
  }

  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vec.cwExtract(cummulativeSize,*(vecs[i]));
    cummulativeSize += vecs[i]->sizeLocal();
  }

  queso_require_equal_to_msg(vec.sizeLocal(), cummulativeSize, "incompatible vector sizes");

  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    result = result && m_sets[i]->contains(*(vecs[i]));
  }
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    delete vecs[i];
  }

  return (result);
}

template<class V, class M>
void ConcatenationSubset<V,M>::centroid(V& vec) const
{
  unsigned int cumulativeSize = 0;
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    V subvec(m_sets[i]->vectorSpace().zeroVector());
    m_sets[i]->centroid(subvec);
    vec.cwSet(cumulativeSize,subvec);

    cumulativeSize += subvec.sizeLocal();
  }

  queso_require_equal_to_msg(vec.sizeLocal(), cumulativeSize, "incompatible vector sizes");
}

template<class V, class M>
void ConcatenationSubset<V,M>::moments(M& mat) const
{
  unsigned int cumulativeSize = 0;
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    const Map & map = m_sets[i]->vectorSpace().map();
    const unsigned int n_columns = map.NumGlobalElements();
    M submat(m_sets[i]->vectorSpace().env(),
             map, n_columns);
    m_sets[i]->moments(submat);
    mat.cwSet(cumulativeSize,cumulativeSize,submat);

    cumulativeSize += n_columns;
  }

  queso_require_equal_to_msg(mat.numCols(), cumulativeSize, "incompatible vector sizes");
}

// I/O methods
template <class V, class M>
void ConcatenationSubset<V,M>::print(std::ostream& os) const
{
  os << "In ConcatenationSubset<V,M>::print()"
     << ": m_sets.size() = " << m_sets.size()
     << std::endl;
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    os << "m_sets[" << i << "] = " << *(m_sets[i]);
    if (i < (m_sets.size()-1)) {
      os << ", ";
    }
  }
  os << std::endl;

  return;
}

}  // End namespace QUESO

template class QUESO::ConcatenationSubset<QUESO::GslVector, QUESO::GslMatrix>;
