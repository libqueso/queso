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

#include <queso/SequentialVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
SequentialVectorRealizer<V,M>::SequentialVectorRealizer(
  const char*                           prefix,
  const BaseVectorSequence<V,M>& chain)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"seq").c_str(),chain.unifiedBoxPlain(),chain.subSequenceSize()),
  m_chain                 (chain),
  m_currentChainPos       (0),
  m_unifiedSampleExpVector(new V(chain.unifiedMeanPlain()          )), // IMPORTANT
  m_unifiedSampleVarVector(new V(chain.unifiedSampleVariancePlain()))  // IMPORTANT
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In SequentialVectorRealizer<V,M>::constructor()"
                            << ": m_chain.subSequenceSize() = " << m_chain.subSequenceSize()
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
SequentialVectorRealizer<V,M>::~SequentialVectorRealizer()
{
  delete m_unifiedSampleVarVector;
  delete m_unifiedSampleExpVector;
}
// Realization-related methods----------------------
template<class V, class M>
void
SequentialVectorRealizer<V,M>::realization(V& nextParamValues) const
{
  m_chain.getPositionValues(m_currentChainPos++,nextParamValues);
  if (m_currentChainPos >= m_subPeriod) m_currentChainPos = 0;

  return;
}
//--------------------------------------------------
template <class V, class M>
const V&
SequentialVectorRealizer<V,M>::unifiedSampleExpVector() const
{
  return *m_unifiedSampleExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
SequentialVectorRealizer<V,M>::unifiedSampleVarVector() const
{
  return *m_unifiedSampleVarVector;
}

}  // End namespace QUESO

template class QUESO::SequentialVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
