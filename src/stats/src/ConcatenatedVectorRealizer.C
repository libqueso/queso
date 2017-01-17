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

#include <queso/ConcatenatedVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
ConcatenatedVectorRealizer<V,M>::ConcatenatedVectorRealizer(
  const char*                           prefix,
  const BaseVectorRealizer<V,M>& realizer1,
  const BaseVectorRealizer<V,M>& realizer2,
  const VectorSet<V,M>&          unifiedImageSet)
  :
  BaseVectorRealizer<V,M>( ((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::min(realizer1.subPeriod(),realizer2.subPeriod()) ), // 2011/Oct/02
  m_realizers(2,(const BaseVectorRealizer<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering ConcatenatedVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_realizers[0] = &realizer1;
  m_realizers[1] = &realizer2;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving ConcatenatedVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
ConcatenatedVectorRealizer<V,M>::ConcatenatedVectorRealizer(
  const char*                                                prefix,
  const std::vector<const BaseVectorRealizer<V,M>* >& realizers,
  unsigned int                                               minPeriod,
  const VectorSet<V,M>&                               unifiedImageSet)
  :
  BaseVectorRealizer<V,M>( ((std::string)(prefix)+"gen").c_str(),unifiedImageSet,minPeriod),
  m_realizers(realizers.size(),(const BaseVectorRealizer<V,M>*) NULL)
{
  for (unsigned int i = 0; i < m_realizers.size(); ++i) {
    m_realizers[i] = realizers[i];
  }
}
// Destructor --------------------------------------
template<class V, class M>
ConcatenatedVectorRealizer<V,M>::~ConcatenatedVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
ConcatenatedVectorRealizer<V,M>::realization(V& nextValues) const
{
  std::vector<V*> vecs(m_realizers.size(),(V*)NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_realizers[i]->unifiedImageSet().vectorSpace().zeroVector());
    //std::cout << "In ConcatenatedVectorRealizer<V,M>::realization: v[i]->sizeLocal() = " << v[i]->sizeLocal() << std::endl;
    m_realizers[i]->realization(*(vecs[i]));
  }

  //std::cout << "In ConcatenatedVectorRealizer<V,M>::realization: nextValues.sizeLocal() = " << nextValues.sizeLocal() << std::endl;
  std::vector<const V*> constVecs(m_realizers.size(),(V*)NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    constVecs[i] = vecs[i];
  }
  nextValues.cwSetConcatenated(constVecs);
  //std::cout << "In ConcatenatedVectorRealizer<V,M>::realization: succeeded" << std::endl;

  for (unsigned int i = 0; i < vecs.size(); ++i) {
    delete vecs[i];
  }

  return;
}

}  // End namespace QUESO

template class QUESO::ConcatenatedVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
