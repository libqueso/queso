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

#include <queso/TransformedScaledCovMatrixTKGroup.h>
#include <queso/InvLogitGaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class V, class M>
TransformedScaledCovMatrixTKGroup<V,M>::TransformedScaledCovMatrixTKGroup(
    const char * prefix,
    const BoxSubset<V,M> & boxSubset,
    const std::vector<double> & scales,
    const M & covMatrix)
  : BaseTKGroup<V, M>(prefix, boxSubset.vectorSpace(), scales),
    m_boxSubset(boxSubset),
    m_originalCovMatrix(covMatrix)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering TransformedScaledCovMatrixTKGroup<V,M>::constructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In TransformedScaledCovMatrixTKGroup<V,M>::constructor()"
                           << ": m_scales.size() = "                << m_scales.size()
                           << ", m_preComputingPositions.size() = " << m_preComputingPositions.size()
                           << ", m_rvs.size() = "                   << m_rvs.size()
                           << ", m_originalCovMatrix = "            << m_originalCovMatrix
                           << std::endl;
  }

  // Set RVs to have zero mean in the Gaussian space
  setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving TransformedScaledCovMatrixTKGroup<V,M>::constructor()"
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
TransformedScaledCovMatrixTKGroup<V,M>::~TransformedScaledCovMatrixTKGroup()
{
}
// Math/Stats methods--------------------------------
template<class V, class M>
bool
TransformedScaledCovMatrixTKGroup<V,M>::symmetric() const
{
  return false;
}
//---------------------------------------------------
template<class V, class M>
const InvLogitGaussianVectorRV<V,M>&
TransformedScaledCovMatrixTKGroup<V,M>::rv(unsigned int stageId) const
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");

  queso_require_msg(m_rvs[0], "m_rvs[0] == NULL");

  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");

  queso_require_msg(m_preComputingPositions[stageId], "m_preComputingPositions[stageId] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In TransformedScaledCovMatrixTKGroup<V,M>::rv1()"
                            << ", stageId = " << stageId
                            << ": about to call m_rvs[0]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageId] // FIX ME: might demand parallelism
                            << std::endl;
  }

  InvLogitGaussianVectorRV<V, M> * invlogit_gaussian =
    dynamic_cast<InvLogitGaussianVectorRV<V, M> * >(m_rvs[0]);

  V transformedPreComputingPositions(*m_preComputingPositions[stageId]);
  transformToGaussianSpace(*m_preComputingPositions[stageId],
      transformedPreComputingPositions);

  invlogit_gaussian->updateLawExpVector(transformedPreComputingPositions);

  return (*invlogit_gaussian);
}
//---------------------------------------------------
template<class V, class M>
const InvLogitGaussianVectorRV<V,M>&
TransformedScaledCovMatrixTKGroup<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  queso_require_greater_equal_msg(m_rvs.size(), stageIds.size(), "m_rvs.size() < stageIds.size()");

  queso_require_msg(m_rvs[stageIds.size()-1], "m_rvs[stageIds.size()-1] == NULL");

  queso_require_greater_msg(m_preComputingPositions.size(), stageIds[0], "m_preComputingPositions.size() <= stageIds[0]");

  queso_require_msg(m_preComputingPositions[stageIds[0]], "m_preComputingPositions[stageIds[0]] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In TransformedScaledCovMatrixTKGroup<V,M>::rv2()"
                            << ", stageIds.size() = " << stageIds.size()
                            << ", stageIds[0] = "     << stageIds[0]
                            << ": about to call m_rvs[stageIds.size()-1]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageIds[0]] // FIX ME: might demand parallelism
                            << std::endl;
  }

  InvLogitGaussianVectorRV<V, M> * invlogit_gaussian =
    dynamic_cast<InvLogitGaussianVectorRV<V, M> * >(m_rvs[stageIds.size()-1]);

  V transformedPreComputingPositions(*m_preComputingPositions[stageIds[0]]);
  transformToGaussianSpace(*m_preComputingPositions[stageIds[0]],
      transformedPreComputingPositions);

  invlogit_gaussian->updateLawExpVector(transformedPreComputingPositions);

  return (*invlogit_gaussian);
}

template <class V, class M>
const InvLogitGaussianVectorRV<V, M> &
TransformedScaledCovMatrixTKGroup<V, M>::rv(const V & position) const
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");
  queso_require_msg(m_rvs[0], "m_rvs[0] == NULL");
  // queso_require_greater_msg(m_preComputingPositions.size(), this->m_stageId, "m_preComputingPositions.size() <= stageId");
  // queso_require_msg(m_preComputingPositions[this->m_stageId], "m_preComputingPositions[stageId] == NULL");

  InvLogitGaussianVectorRV<V, M> * invlogit_gaussian =
    dynamic_cast<InvLogitGaussianVectorRV<V, M> * >(m_rvs[this->m_stageId]);

  V transformedPreComputingPositions(position);
  transformToGaussianSpace(position, transformedPreComputingPositions);
  invlogit_gaussian->updateLawExpVector(transformedPreComputingPositions);

  return (*invlogit_gaussian);
}

//---------------------------------------------------
template<class V, class M>
void
TransformedScaledCovMatrixTKGroup<V,M>::updateLawCovMatrix(const M& covMatrix)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In TransformedScaledCovMatrixTKGroup<V,M>::updateLawCovMatrix()"
                              << ", m_scales.size() = " << m_scales.size()
                              << ", i = "               << i
                              << ", m_scales[i] = "     << m_scales[i]
                              << ", factor = "          << factor
                              << ": about to call m_rvs[i]->updateLawCovMatrix()"
                              << ", covMatrix = \n" << factor*covMatrix // FIX ME: might demand parallelism
                              << std::endl;
    }

    InvLogitGaussianVectorRV<V, M> * invlogit_gaussian =
      dynamic_cast<InvLogitGaussianVectorRV<V, M> * >(m_rvs[i]);

    invlogit_gaussian->updateLawCovMatrix(factor*covMatrix);
  }

  return;
}

// Misc methods -------------------------------------
template<class V, class M>
bool
TransformedScaledCovMatrixTKGroup<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering TransformedScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  BaseTKGroup<V,M>::setPreComputingPosition(position,stageId);
  //setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In TransformedScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ", position = "        << position
                           << ", stageId = "         << stageId
                           << ": preComputingPos = " << *m_preComputingPositions[stageId];
    if (stageId < m_scales.size()) {
      *m_env.subDisplayFile() << ", factor = " << 1./m_scales[stageId]/m_scales[stageId];
    }
    if (stageId < m_rvs.size()) {
      const InvLogitGaussianJointPdf<V,M>* pdfPtr = dynamic_cast< const InvLogitGaussianJointPdf<V,M>* >(&(m_rvs[stageId]->pdf()));
      *m_env.subDisplayFile() << ", rvCov = " << pdfPtr->lawCovMatrix(); // FIX ME: might demand parallelism
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving TransformedScaledCovMatrixTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  return true;
}
//---------------------------------------------------
template<class V, class M>
void
TransformedScaledCovMatrixTKGroup<V,M>::clearPreComputingPositions()
{
  BaseTKGroup<V,M>::clearPreComputingPositions();
  return;
}

template <class V, class M>
unsigned int
TransformedScaledCovMatrixTKGroup<V, M>::set_dr_stage(unsigned int stageId)
{
  unsigned int old_stageId = this->m_stageId;
  this->m_stageId = stageId;
  return old_stageId;
}

// Private methods------------------------------------
template<class V, class M>
void
TransformedScaledCovMatrixTKGroup<V,M>::setRVsWithZeroMean()
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");

  queso_require_equal_to_msg(m_rvs.size(), m_scales.size(), "m_rvs.size() != m_scales.size()");

  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    queso_require_msg(!(m_rvs[i]), "m_rvs[i] != NULL");
    m_rvs[i] = new InvLogitGaussianVectorRV<V,M>(m_prefix.c_str(),
        m_boxSubset, m_vectorSpace->zeroVector(),
        factor*m_originalCovMatrix);
  }

  return;
}

template<class V, class M>
void
TransformedScaledCovMatrixTKGroup<V,M>::print(std::ostream& os) const
{
  BaseTKGroup<V,M>::print(os);
  return;
}

template<class V, class M>
void
TransformedScaledCovMatrixTKGroup<V, M>::transformToGaussianSpace(
    const V & physicalPoint, V & transformedPoint) const
{
  V min_domain_bounds(this->m_boxSubset.minValues());
  V max_domain_bounds(this->m_boxSubset.maxValues());

  for (unsigned int i = 0; i < transformedPoint.sizeLocal(); i++) {
    double min_val = min_domain_bounds[i];
    double max_val = max_domain_bounds[i];

    if (queso_isfinite(min_val) &&
        queso_isfinite(max_val)) {
        // Left- and right-hand sides are finite.  Do full transform.
        transformedPoint[i] = std::log(physicalPoint[i] - min_val) -
            std::log(max_val - physicalPoint[i]);
    }
    else if (queso_isfinite(min_val) &&
             !queso_isfinite(max_val)) {
      // Left-hand side finite, but right-hand side is not.
      // Do only left-hand transform.
      transformedPoint[i] = std::log(physicalPoint[i] - min_val);
    }
    else if (!queso_isfinite(min_val) &&
             queso_isfinite(max_val)) {
      // Right-hand side is finite, but left-hand side is not.
      // Do only right-hand transform.
      transformedPoint[i] = -std::log(max_val - physicalPoint[i]);
    }
    else {
      // No transform.
      transformedPoint[i] = physicalPoint[i];
    }
  }
}

}  // End namespace QUESO

template class QUESO::TransformedScaledCovMatrixTKGroup<QUESO::GslVector, QUESO::GslMatrix>;
