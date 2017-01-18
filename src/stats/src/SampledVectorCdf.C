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

#include <queso/SampledVectorCdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V,class M>
SampledVectorCdf<V,M>::SampledVectorCdf(
  const char*                          prefix,
  const ArrayOfOneDGrids <V,M>& oneDGrids,
  const ArrayOfOneDTables<V,M>& cdfValues)
  :
  BaseVectorCdf<V,M>(prefix,oneDGrids.rowSpace()),
  m_cdfs(m_pdfSupport.vectorSpace().map(),1)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SampledVectorCdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  char strI[65];
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    sprintf(strI,"%u_",i);
    m_cdfs(i,0) = new SampledScalarCdf<double>(m_env,
                                                      ((std::string)(m_prefix)+strI).c_str(),
                                                      oneDGrids.grid(i),
                                                      cdfValues.oneDTable(i));
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SampledVectorCdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
SampledVectorCdf<V,M>::~SampledVectorCdf()
{
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    if (m_cdfs(i,0)) delete m_cdfs(i,0);
  }
}
// Math methods--------------------------------------
template<class V, class M>
void
SampledVectorCdf<V,M>::values(
  const V& paramValues,
        V& cdfVec) const
{
  queso_error_msg("incomplete code");
  return;
}
// --------------------------------------------------
template<class V, class M>
const BaseScalarCdf<double>&
SampledVectorCdf<V,M>::cdf(unsigned int rowId) const
{
  queso_require_less_msg(rowId, m_pdfSupport.vectorSpace().dimLocal(), "rowId is out of range");

  SampledVectorCdf<V,M>* tmp = const_cast<SampledVectorCdf<V,M>*>(this);
  return *(tmp->m_cdfs(rowId,0));

}
// I/O methods---------------------------------------
template <class V, class M>
void
SampledVectorCdf<V,M>::print(std::ostream& os) const
{
  SampledVectorCdf<V,M>* tmp = const_cast<SampledVectorCdf<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    os << (*tmp->m_cdfs(i,0))
       << std::endl;
  }

  return;
}
//---------------------------------------------------
template<class V, class M>
void
SampledVectorCdf<V,M>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  queso_require_greater_equal_msg(m_env.subRank(), 0, "unexpected subRank");

  SampledVectorCdf<V,M>* tmp = const_cast<SampledVectorCdf<V,M>*>(this);
  char compId[16+1];
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    sprintf(compId,"%d",i);
    tmp->m_cdfs(i,0)->subWriteContents(varNamePrefix+"comp"+compId,fileName,fileType,allowedSubEnvIds);
  }

  return;
}

}  // End namespace QUESO
