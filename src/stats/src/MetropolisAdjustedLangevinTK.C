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

#include <queso/MetropolisAdjustedLangevinTK.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianJointPdf.h>
#include <queso/BayesianJointPdf.h>

namespace QUESO {

template <class V, class M>
MetropolisAdjustedLangevinTK<V, M>::MetropolisAdjustedLangevinTK(
  const char * prefix,
  const BayesianJointPdf<V, M> & targetPdf,
  const std::vector<double> & scales,
  const M & covMatrix)
  :
  BaseTKGroup<V, M>(prefix, targetPdf.domainSet().vectorSpace(), scales),
  m_originalCovMatrix(covMatrix),
  m_targetPdf(targetPdf),
  m_time_step(1.0)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering MetropolisAdjustedLangevinTK<V, M>::constructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In MetropolisAdjustedLangevinTK<V, M>::constructor()"
                           << ": m_scales.size() = "                << m_scales.size()
                           << ", m_preComputingPositions.size() = " << m_preComputingPositions.size()
                           << ", m_rvs.size() = "                   << m_rvs.size()
                           << ", m_originalCovMatrix = "            << m_originalCovMatrix
                           << std::endl;
  }

  setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving MetropolisAdjustedLangevinTK<V, M>::constructor()"
                           << std::endl;
  }
}

template <class V, class M>
MetropolisAdjustedLangevinTK<V, M>::~MetropolisAdjustedLangevinTK()
{
}

template <class V, class M>
bool
MetropolisAdjustedLangevinTK<V, M>::symmetric() const
{
  return false;
}

template <class V, class M>
const GaussianVectorRV<V, M> &
MetropolisAdjustedLangevinTK<V, M>::rv(unsigned int stageId) const
{
  queso_require_not_equal_to(m_rvs.size(), 0);
  queso_require(m_rvs[0]);
  queso_require_greater(m_preComputingPositions.size(), stageId);
  queso_require(m_preComputingPositions[stageId]);

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In MetropolisAdjustedLangevinTK<V, M>::rv1()"
                            << ", stageId = " << stageId
                            << ": about to call m_rvs[0]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageId] // FIX ME: might demand parallelism
                            << std::endl;
  }

  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[0]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageId]);

  return (*gaussian_rv);
}

template <class V, class M>
const GaussianVectorRV<V, M> &
MetropolisAdjustedLangevinTK<V, M>::rv(const std::vector<unsigned int> & stageIds)
{
  queso_require_greater_equal(m_rvs.size(), stageIds.size());
  queso_require(m_rvs[stageIds.size()-1]);
  queso_require_greater(m_preComputingPositions.size(), stageIds[0]);
  queso_require(m_preComputingPositions[stageIds[0]]);

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In MetropolisAdjustedLangevinTK<V, M>::rv2()"
                            << ", stageIds.size() = " << stageIds.size()
                            << ", stageIds[0] = "     << stageIds[0]
                            << ": about to call m_rvs[stageIds.size()-1]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageIds[0]] // FIX ME: might demand parallelism
                            << std::endl;
  }

  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[stageIds.size()-1]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageIds[0]]);

  return (*gaussian_rv);
}

template <class V, class M>
const GaussianVectorRV<V, M> &
MetropolisAdjustedLangevinTK<V, M>::rv(const V & position) const
{
  queso_require_not_equal_to(m_rvs.size(), 0);
  queso_require(m_rvs[0]);

  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[this->m_stageId]);

  // 'position' is a position in the chain?  Hopefully?  Anyway, assume it's a
  // position in the chain and then that means we need to modify it slightly so
  // to get the transition distribution for the next state in the chain.

  V grad(this->m_targetPdf.domainSet().vectorSpace().zeroVector());

  // Get the gradient of the log-posterior.  This is so inefficient it's
  // painful.  We should be caching the gradient evaluations.
  this->m_targetPdf.lnValue(position, grad);

  // Euler time-step
  grad *= 0.5 * this->m_time_step;

  // Add on current position
  grad += position;

  // Update the mean of the transition kernel.  The vector gets copied, so
  // we're ok.
  gaussian_rv->updateLawExpVector(grad);

  return (*gaussian_rv);
}

template <class V, class M>
void
MetropolisAdjustedLangevinTK<V, M>::updateLawCovMatrix(const M & covMatrix)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In MetropolisAdjustedLangevinTK<V, M>::updateLawCovMatrix()"
                              << ", m_scales.size() = " << m_scales.size()
                              << ", i = "               << i
                              << ", m_scales[i] = "     << m_scales[i]
                              << ", factor = "          << factor
                              << ": about to call m_rvs[i]->updateLawCovMatrix()"
                              << ", covMatrix = \n" << factor*covMatrix // FIX ME: might demand parallelism
                              << std::endl;
    }
    dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[i])->updateLawCovMatrix(factor*m_time_step*covMatrix);
  }
}

template <class V, class M>
bool
MetropolisAdjustedLangevinTK<V, M>::setPreComputingPosition(const V & position, unsigned int stageId)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering MetropolisAdjustedLangevinTK<V, M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  BaseTKGroup<V, M>::setPreComputingPosition(position, stageId);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In MetropolisAdjustedLangevinTK<V, M>::setPreComputingPosition()"
                           << ", position = "        << position
                           << ", stageId = "         << stageId
                           << ": preComputingPos = " << *m_preComputingPositions[stageId];
    if (stageId < m_scales.size()) {
      *m_env.subDisplayFile() << ", factor = " << 1./m_scales[stageId]/m_scales[stageId];
    }
    if (stageId < m_rvs.size()) {
      const GaussianJointPdf<V, M>* pdfPtr = dynamic_cast< const GaussianJointPdf<V, M>* >(&(m_rvs[stageId]->pdf()));
      *m_env.subDisplayFile() << ", rvCov = " << pdfPtr->lawCovMatrix(); // FIX ME: might demand parallelism
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving MetropolisAdjustedLangevinTK<V, M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  return true;
}

template <class V, class M>
void
MetropolisAdjustedLangevinTK<V, M>::clearPreComputingPositions()
{
  BaseTKGroup<V, M>::clearPreComputingPositions();
}

template <class V, class M>
void
MetropolisAdjustedLangevinTK<V, M>::setRVsWithZeroMean()
{
  queso_require_not_equal_to(m_rvs.size(), 0);
  queso_require_equal_to(m_rvs.size(), m_scales.size());

  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    queso_require(!(m_rvs[i]));
    m_rvs[i] = new GaussianVectorRV<V, M>(m_prefix.c_str(),
                                         *m_vectorSpace,
                                         m_vectorSpace->zeroVector(),
                                         factor*m_time_step*m_originalCovMatrix);
  }
}

template <class V, class M>
void
MetropolisAdjustedLangevinTK<V, M>::print(std::ostream & os) const
{
  BaseTKGroup<V, M>::print(os);
}

template class MetropolisAdjustedLangevinTK<GslVector, GslMatrix>;

}  // End namespace QUESO
