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

#include <queso/ConcatenatedVectorRV.h>
#include <queso/ConcatenatedVectorRealizer.h>
#include <queso/ConcatenatedJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <limits>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
ConcatenatedVectorRV<V,M>::ConcatenatedVectorRV(
  const char*                     prefix,
  const BaseVectorRV<V,M>& rv1,
  const BaseVectorRV<V,M>& rv2,
  const VectorSet<V,M>&    imageSet)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"concat").c_str(),imageSet),
  m_rvs                   (2,(const BaseVectorRV      <V,M>*) NULL),
  m_pdfs                  (2,(const BaseJointPdf      <V,M>*) NULL),
  m_realizers             (2,(const BaseVectorRealizer<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ConcatenatedVectorRV<V,M>::constructor(1)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_rvs[0]       = &rv1;
  m_rvs[1]       = &rv2;
  m_pdfs[0]      = &(m_rvs[0]->pdf());
  m_pdfs[1]      = &(m_rvs[1]->pdf());
  m_realizers[0] = m_rvs[0]->has_realizer() ?
                   &(m_rvs[0]->realizer()) : NULL;
  m_realizers[1] = m_rvs[1]->has_realizer() ?
                   &(m_rvs[1]->realizer()) : NULL;

  m_pdf          = new ConcatenatedJointPdf<V,M>(m_prefix.c_str(),
                                                        *(m_pdfs[0]),
                                                        *(m_pdfs[1]),
                                                        m_imageSet);

  // Iff we have all sub-realizers, we can make our own realizer
  m_realizer     = (m_rvs[0]->has_realizer() &&
                    m_rvs[1]->has_realizer()) ?
                     new ConcatenatedVectorRealizer<V,M>(m_prefix.c_str(),
                                                         *(m_realizers[0]),
                                                         *(m_realizers[1]),
                                                         m_imageSet) :
                     NULL;

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ConcatenatedVectorRV<V,M>::constructor(1)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
ConcatenatedVectorRV<V,M>::ConcatenatedVectorRV(
  const char*                                          prefix,
  const std::vector<const BaseVectorRV<V,M>* >& rvs,
  const VectorSet<V,M>&                         imageSet)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"concat").c_str(),imageSet),
  m_rvs                   (rvs.size(),(const BaseVectorRV      <V,M>*) NULL),
  m_pdfs                  (rvs.size(),(const BaseJointPdf      <V,M>*) NULL),
  m_realizers             (rvs.size(),(const BaseVectorRealizer<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ConcatenatedVectorRV<V,M>::constructor(2)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  bool have_all_subrealizers = true;  // maybe
  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    m_rvs [i]      = rvs[i];
    m_pdfs[i]      = &(m_rvs[i]->pdf());

    // If our sub-RV has no realizer, we leave our pointer to it set
    // to NULL
    if (m_rvs[i]->has_realizer())
      m_realizers[i] = &(m_rvs[i]->realizer());
    else
      have_all_subrealizers = false;
  }

  m_pdf        = new ConcatenatedJointPdf<V,M>(m_prefix.c_str(),
                                                      m_pdfs,
                                                      m_imageSet);

  unsigned int minPeriod = std::numeric_limits<unsigned int>::max();
  for (unsigned int i = 0; i < m_realizers.size(); ++i) {
    if (m_realizers[i] && minPeriod > m_realizers[i]->subPeriod()) {
      minPeriod = m_realizers[i]->subPeriod();
    }
  }

  m_realizer   = have_all_subrealizers ?
                   new ConcatenatedVectorRealizer<V,M>(m_prefix.c_str(),
                                                       m_realizers,
                                                       minPeriod,
                                                       m_imageSet) :
                   NULL;
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ConcatenatedVectorRV<V,M>::constructor(2)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
ConcatenatedVectorRV<V,M>::~ConcatenatedVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
//--------------------------------------------------
template <class V, class M>
void
ConcatenatedVectorRV<V,M>::print(std::ostream& os) const
{
  os << "ConcatenatedVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
