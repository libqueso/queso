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

#include <queso/GenericVectorRV.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V, class M>
GenericVectorRV<V,M>::GenericVectorRV(
  const char*                     prefix,
  const VectorSet <V,M>& imageSet)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"gen").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
GenericVectorRV<V,M>::GenericVectorRV(
  const char*                           prefix,
  const VectorSet         <V,M>& imageSet,
  BaseJointPdf      <V,M>& pdf,
  BaseVectorRealizer<V,M>& realizer,
  const BaseVectorCdf     <V,M>& subCdf,
  const BaseVectorCdf     <V,M>& unifiedCdf,
  const BaseVectorMdf     <V,M>& mdf)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"gen").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = &pdf;
  m_realizer   = &realizer;
  m_subCdf     = &subCdf;
  m_unifiedCdf = &unifiedCdf;
  m_mdf        = &mdf;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
GenericVectorRV<V,M>::~GenericVectorRV()
{
}
// Random variable-handling methods ----------------
template<class V, class M>
void
GenericVectorRV<V,M>::setPdf(BaseJointPdf<V,M>& pdf)
{
  m_pdf = &pdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
GenericVectorRV<V,M>::setRealizer(BaseVectorRealizer<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
GenericVectorRV<V,M>::setSubCdf(BaseVectorCdf<V,M>& subCdf)
{
  m_subCdf = &subCdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
GenericVectorRV<V,M>::setUnifiedCdf(BaseVectorCdf<V,M>& unifiedCdf)
{
  m_unifiedCdf = &unifiedCdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
GenericVectorRV<V,M>::setMdf(BaseVectorMdf<V,M>& mdf)
{
  m_mdf = &mdf;
  return;
}
//--------------------------------------------------
template <class V, class M>
void
GenericVectorRV<V,M>::print(std::ostream& os) const
{
  os << "GenericVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
