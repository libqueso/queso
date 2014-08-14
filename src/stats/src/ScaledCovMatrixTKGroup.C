//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <queso/ScaledCovMatrixTKGroup.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianJointPdf.h>

namespace QUESO {

// Default constructor ------------------------------
template<class V, class M>
ScaledCovMatrixTKGroup<V,M>::ScaledCovMatrixTKGroup(
  const char*                    prefix,
  const VectorSpace<V,M>& vectorSpace, // FIX ME: vectorSubset ???
  const std::vector<double>&     scales,
  const M&                       covMatrix)
  :
  BaseTKGroup<V,M>(prefix,vectorSpace,scales),
  m_originalCovMatrix    (covMatrix)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering ScaledCovMatrixTKGroup<V,M>::constructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In ScaledCovMatrixTKGroup<V,M>::constructor()"
                           << ": m_scales.size() = "                << m_scales.size()
                           << ", m_preComputingPositions.size() = " << m_preComputingPositions.size()
                           << ", m_rvs.size() = "                   << m_rvs.size()
                           << ", m_originalCovMatrix = "            << m_originalCovMatrix
                           << std::endl;
  }

  setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving ScaledCovMatrixTKGroup<V,M>::constructor()"
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
ScaledCovMatrixTKGroup<V,M>::~ScaledCovMatrixTKGroup()
{
}
// Math/Stats methods--------------------------------
template<class V, class M>
bool
ScaledCovMatrixTKGroup<V,M>::symmetric() const
{
  return true;
}
//---------------------------------------------------
template<class V, class M>
const GaussianVectorRV<V,M>&
ScaledCovMatrixTKGroup<V,M>::rv(unsigned int stageId)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() == 0,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv1()",
                      "m_rvs.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs[0] == NULL, // Yes, '0', because that is the id used below
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv1()",
                      "m_rvs[0] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv1()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] == NULL,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv1()",
                      "m_preComputingPositions[stageId] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In ScaledCovMatrixTKGroup<V,M>::rv1()"
                            << ", stageId = " << stageId
                            << ": about to call m_rvs[0]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageId] // FIX ME: might demand parallelism
                            << std::endl;
  }
  m_rvs[0]->updateLawExpVector(*m_preComputingPositions[stageId]);

  return (*m_rvs[0]);
}
//---------------------------------------------------
template<class V, class M>
const GaussianVectorRV<V,M>&
ScaledCovMatrixTKGroup<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < stageIds.size(),
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv2()",
                      "m_rvs.size() < stageIds.size()");

  UQ_FATAL_TEST_MACRO(m_rvs[stageIds.size()-1] == NULL,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv2()",
                      "m_rvs[stageIds.size()-1] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageIds[0],
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv2()",
                      "m_preComputingPositions.size() <= stageIds[0]");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageIds[0]] == NULL,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::rv2()",
                      "m_preComputingPositions[stageIds[0]] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In ScaledCovMatrixTKGroup<V,M>::rv2()"
                            << ", stageIds.size() = " << stageIds.size()
                            << ", stageIds[0] = "     << stageIds[0]
                            << ": about to call m_rvs[stageIds.size()-1]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageIds[0]] // FIX ME: might demand parallelism
                            << std::endl;
  }
  m_rvs[stageIds.size()-1]->updateLawExpVector(*m_preComputingPositions[stageIds[0]]);

  return (*m_rvs[stageIds.size()-1]);
}
//---------------------------------------------------
template<class V, class M>
void
ScaledCovMatrixTKGroup<V,M>::updateLawCovMatrix(const M& covMatrix)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In ScaledCovMatrixTKGroup<V,M>::updateLawCovMatrix()"
                              << ", m_scales.size() = " << m_scales.size()
                              << ", i = "               << i
                              << ", m_scales[i] = "     << m_scales[i]
                              << ", factor = "          << factor
                              << ": about to call m_rvs[i]->updateLawCovMatrix()"
                              << ", covMatrix = \n" << factor*covMatrix // FIX ME: might demand parallelism
                              << std::endl;
    }
    m_rvs[i]->updateLawCovMatrix(factor*covMatrix);
  }

  return;
}

// Misc methods -------------------------------------
template<class V, class M>
bool
ScaledCovMatrixTKGroup<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering ScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  BaseTKGroup<V,M>::setPreComputingPosition(position,stageId);
  //setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In ScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ", position = "        << position
                           << ", stageId = "         << stageId
                           << ": preComputingPos = " << *m_preComputingPositions[stageId];
    if (stageId < m_scales.size()) {
      *m_env.subDisplayFile() << ", factor = " << 1./m_scales[stageId]/m_scales[stageId];
    }
    if (stageId < m_rvs.size()) {
      const GaussianJointPdf<V,M>* pdfPtr = dynamic_cast< const GaussianJointPdf<V,M>* >(&(m_rvs[stageId]->pdf()));
      *m_env.subDisplayFile() << ", rvCov = " << pdfPtr->lawCovMatrix(); // FIX ME: might demand parallelism
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving ScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  return true;
}
//---------------------------------------------------
template<class V, class M>
void
ScaledCovMatrixTKGroup<V,M>::clearPreComputingPositions()
{
  BaseTKGroup<V,M>::clearPreComputingPositions();
  return;
}


// Private methods------------------------------------
template<class V, class M>
void
ScaledCovMatrixTKGroup<V,M>::setRVsWithZeroMean()
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() == 0,
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::setRVsWithZeroMean()",
                      "m_rvs.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs.size() != m_scales.size(),
                      m_env.worldRank(),
                      "ScaledCovMatrixTKGroup<V,M>::setRVsWithZeroMean()",
                      "m_rvs.size() != m_scales.size()");

  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    UQ_FATAL_TEST_MACRO(m_rvs[i] != NULL,
                        m_env.worldRank(),
                        "ScaledCovMatrixTKGroup<V,M>::setRVsWithZeroMean()",
                        "m_rvs[i] != NULL");
    m_rvs[i] = new GaussianVectorRV<V,M>(m_prefix.c_str(),
                                                *m_vectorSpace,
                                                m_vectorSpace->zeroVector(),
                                                factor*m_originalCovMatrix);
  }

  return;
}
// I/O methods---------------------------------------
template<class V, class M>
void
ScaledCovMatrixTKGroup<V,M>::print(std::ostream& os) const
{
  BaseTKGroup<V,M>::print(os);
  return;
}

}  // End namespace QUESO

template class QUESO::ScaledCovMatrixTKGroup<QUESO::GslVector, QUESO::GslMatrix>;
