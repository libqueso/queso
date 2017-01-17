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

#include <queso/HessianCovMatricesTKGroup.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor ------------------------------
template<class V, class M>
HessianCovMatricesTKGroup<V,M>::HessianCovMatricesTKGroup(
  const char*                                   prefix,
  const VectorSpace<V,M>&                vectorSpace,
  const std::vector<double>&                    scales,
  const ScalarFunctionSynchronizer<V,M>& targetPdfSynchronizer)
  :
  BaseTKGroup<V,M>(prefix,vectorSpace,scales),
  m_targetPdfSynchronizer(targetPdfSynchronizer),
  m_originalNewtonSteps  (scales.size()+1,NULL), // Yes, +1
  m_originalCovMatrices  (scales.size()+1,NULL)  // Yes, +1
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering HessianCovMatricesTKGroup<V,M>::constructor()"
                           << std::endl;
  }

  m_rvs.resize(scales.size()+1,NULL); // Yes, +1 (IMPORTANT: overwrite initialization done by BaseTKGroup<V,M>)

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::constructor()"
                           << ": m_scales.size() = "                   << m_scales.size()
                           << ", m_preComputingPositions.size() = "    << m_preComputingPositions.size()
                           << ", m_rvs.size() = "                      << m_rvs.size()
                           << ", m_originalNewtonSteps.size() = "      << m_originalNewtonSteps.size()
                           << ", m_originalCovMatrices.size() = "      << m_originalCovMatrices.size()
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving HessianCovMatricesTKGroup<V,M>::constructor()"
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
HessianCovMatricesTKGroup<V,M>::~HessianCovMatricesTKGroup()
{
}
// Math/Stats methods--------------------------------
template<class V, class M>
bool
HessianCovMatricesTKGroup<V,M>::symmetric() const
{
  return false;
}
//---------------------------------------------------
template<class V, class M>
const GaussianVectorRV<V,M>&
HessianCovMatricesTKGroup<V,M>::rv(unsigned int stageId) const
{
  queso_require_greater_msg(m_rvs.size(), stageId, "m_rvs.size() <= stageId");

  queso_require_msg(m_rvs[stageId], "m_rvs[stageId] == NULL");

  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");

  queso_require_msg(m_preComputingPositions[stageId], "m_preComputingPositions[stageId] == NULL");

  GaussianVectorRV<V, M> * gaussian_rv =
    dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[stageId]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageId] + *m_originalNewtonSteps[stageId]);

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In HessianCovMatrixTKGroup<V,M>::rv1()"
                            << ", stageId = " << stageId
                            << ": about to call m_rvs[stageId]->updateLawCovMatrix()"
                            << ", covMatrix = \n" << *m_originalCovMatrices[stageId] // FIX ME: might demand parallelism
                            << std::endl;
  }

  gaussian_rv->updateLawCovMatrix(*m_originalCovMatrices[stageId]);

  return *gaussian_rv;
}
//---------------------------------------------------
template<class V, class M>
const GaussianVectorRV<V,M>&
HessianCovMatricesTKGroup<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  queso_require_greater_msg(m_rvs.size(), stageIds[0], "m_rvs.size() <= stageIds[0]");

  queso_require_msg(m_rvs[stageIds[0]], "m_rvs[stageIds[0]] == NULL");

  queso_require_greater_msg(m_preComputingPositions.size(), stageIds[0], "m_preComputingPositions.size() <= stageIds[0]");

  queso_require_msg(m_preComputingPositions[stageIds[0]], "m_preComputingPositions[stageIds[0]] == NULL");

  double factor = 1./m_scales[stageIds.size()-1]/m_scales[stageIds.size()-1];

  GaussianVectorRV<V, M> * gaussian_rv =
    dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[stageIds[0]]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageIds[0]] + factor*(*m_originalNewtonSteps[stageIds[0]]));

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In HessianCovMatrixTKGroup<V,M>::rv2()"
                            << ", stageIds.size() = " << stageIds.size()
                            << ", stageIds[0] = "     << stageIds[0]
                            << ", factor = "          << factor
                            << ": about to call m_rvs[stageIds[0]]->updateLawCovVector()"
                            << ", covMatrix = \n" << factor*(*m_originalCovMatrices[stageIds[0]]) // FIX ME: might demand parallelism
                            << std::endl;
  }
  gaussian_rv->updateLawCovMatrix(factor*(*m_originalCovMatrices[stageIds[0]]));

  return *gaussian_rv;
}

template <class V, class M>
const GaussianVectorRV<V, M> &
HessianCovMatricesTKGroup<V, M>::rv(const V & position) const
{
  queso_require_greater_msg(m_rvs.size(), this->m_stageId, "m_rvs.size() <= stageId");
  queso_require_msg(m_rvs[this->m_stageId], "m_rvs[stageId] == NULL");
  // queso_require_greater_msg(m_preComputingPositions.size(), this->m_stageId, "m_preComputingPositions.size() <= stageId");
  // queso_require_msg(m_preComputingPositions[this->m_stageId], "m_preComputingPositions[stageId] == NULL");

  GaussianVectorRV<V, M> * gaussian_rv =
    dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[this->m_stageId]);

  gaussian_rv->updateLawExpVector(position + *m_originalNewtonSteps[this->m_stageId]);

  gaussian_rv->updateLawCovMatrix(*m_originalCovMatrices[this->m_stageId]);

  return *gaussian_rv;
}

// Misc methods--------------------------------------
template<class V, class M>
bool
HessianCovMatricesTKGroup<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  bool validPreComputingPosition = true;

  // Verify consistency of sizes
  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");

  queso_require_equal_to_msg(m_preComputingPositions.size(), m_rvs.size(), "m_preComputingPositions.size() != m_rvs.size()");

  queso_require_equal_to_msg(m_preComputingPositions.size(), m_originalNewtonSteps.size(), "m_preComputingPositions.size() != m_originalNewtonSteps.size()");

  queso_require_equal_to_msg(m_preComputingPositions.size(), m_originalCovMatrices.size(), "m_preComputingPositions.size() != m_originalCovMatrices.size()");

  // Verify data is not null
  queso_require_msg(!(m_preComputingPositions[stageId]), "m_preComputingPositions[stageId] != NULL");

  queso_require_msg(!(m_rvs[stageId]), "m_rvs[stageId] != NULL");

  queso_require_msg(!(m_originalNewtonSteps[stageId]), "m_originalNewtonSteps[stageId] != NULL");

  queso_require_msg(!(m_originalCovMatrices[stageId]), "m_originalCovMatrices[stageId] != NULL");

  BaseTKGroup<V,M>::setPreComputingPosition(position,stageId);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                           << ", position = "                          << position
                           << ", stageId = "                           << stageId
                           << ": m_originalNewtonSteps.size() = "      << m_originalNewtonSteps.size()
                           << ", m_originalCovMatrices.size() = "      << m_originalCovMatrices.size()
                           << ", m_preComputingPositions.size() = "    << m_preComputingPositions.size()
                           << ", m_rvs.size() = "                      << m_rvs.size()
                           << std::endl;
  }

  if (m_targetPdfSynchronizer.domainSet().contains(position)) {
    M* tmpHessian = m_vectorSpace->newMatrix();
    M* tmpCovMat  = m_vectorSpace->newMatrix();
    V* tmpGrad    = m_vectorSpace->newVector();

    double logPrior = 0.;
    double logLikelihood = 0.;
    double logTarget = 0.;
    logTarget = m_targetPdfSynchronizer.callFunction(&position, // Might demand parallel environment
                                                     NULL,
                                                     tmpGrad,
                                                     tmpHessian,
                                                     NULL,
                                                     &logPrior,
                                                     &logLikelihood);
    if (logTarget) {}; // just to remove compiler warning

    // IMPORTANT: covariance matrix = (Hessian)^{-1} !!!
    V unitVector(m_vectorSpace->zeroVector());
    V multVector(m_vectorSpace->zeroVector());
    for (unsigned int j = 0; j < tmpHessian->numCols(); ++j) {
      if (j > 0) unitVector[j-1] = 0.;
      unitVector[j] = 1.;
      tmpHessian->invertMultiply(unitVector, multVector);
      for (unsigned int i = 0; i < tmpHessian->numRowsLocal(); ++i) {
        (*tmpCovMat)(i,j) = multVector[i];
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
      *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                             << ", position = "  << position
                             << ", stageId = "   << stageId
                             << ":\n H = "       << *tmpHessian
                             << "\n H^{-1} = "   << *tmpCovMat
                             << "\n H*H^{-1} = " << (*tmpHessian)*(*tmpCovMat)
                             << "\n H^{-1}*H = " << (*tmpCovMat)*(*tmpHessian)
                             << std::endl;
    }

    // Force covariance matrix to be symmetric, as the Hessian (supposedly) is
    *tmpCovMat = .5*((*tmpCovMat) + tmpCovMat->transpose());

    // Test if covariance matrix is positive definite
    M lowerChol(*tmpCovMat);
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                              << ", position = "  << position
                              << ", stageId = "   << stageId
                              << ": calling lowerChol.chol()"
                              << ", lowerChol = " << lowerChol
                              << std::endl;
    }
    int iRC = lowerChol.chol();
    if (iRC) {
      std::cerr << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition(): chol failed\n";
    }
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                              << ", position = "  << position
                              << ", stageId = "   << stageId
                              << ": got lowerChol.chol() with iRC = " << iRC
                              << std::endl;
    }

    bool covIsPositiveDefinite = !iRC;

    if (covIsPositiveDefinite) {

      //                    m_env.worldRank(),
      //                    "HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()",
      //                    "stageId is too large for m_scales");
      //double factor = 1./m_scales[stageId]/m_scales[stageId];
      //*tmpCovMat *= factor;

      m_originalNewtonSteps[stageId] = new V(-1.*(*tmpCovMat)*(*tmpGrad));
      m_originalCovMatrices[stageId] = new M(*tmpCovMat);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
        *m_env.subDisplayFile() << "In HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                               << ", position = "        << position
                               << ", stageId = "         << stageId
                               << ", about to instantiate a Gaussian RV"
                               << ": tmpHessian = "      << *tmpHessian
                               << ", preComputingPos = " << *m_preComputingPositions[stageId]
                               << ", tmpCovMat = "       << *tmpCovMat
                               << ", tmpGrad = "         << *tmpGrad
                               << ", preComputedPos = "  << *m_preComputingPositions[stageId] + *m_originalNewtonSteps[stageId]
                             //<< ", factor = "          << factor
                             //<< ", rvCov = "           << factor*(*tmpCovMat)
                               << std::endl;
      }
      m_rvs[stageId] = new GaussianVectorRV<V,M>(m_prefix.c_str(),
                                                        *m_vectorSpace,
                                                        *m_preComputingPositions[stageId] + *m_originalNewtonSteps[stageId],
                                                        *m_originalCovMatrices[stageId]);
    }
    else {
      validPreComputingPosition = false;
    }

    delete tmpGrad;
    delete tmpCovMat;
    delete tmpHessian;
  }
  else {
    validPreComputingPosition = false;
  }

  if (validPreComputingPosition == false) {
    // Put "default" values on variables
    V tmpGrad  (m_vectorSpace->zeroVector());
    M tmpCovMat(tmpGrad,1.); // = identity matrix
    m_originalNewtonSteps[stageId] = new V(-1.*tmpCovMat*tmpGrad);
    m_originalCovMatrices[stageId] = new M(tmpCovMat);
    m_rvs[stageId] = new GaussianVectorRV<V,M>(m_prefix.c_str(),
                                                      *m_vectorSpace,
                                                      *m_preComputingPositions[stageId],
                                                      tmpCovMat);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving HessianCovMatricesTKGroup<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  return validPreComputingPosition;
}
//---------------------------------------------------
template<class V, class M>
void
HessianCovMatricesTKGroup<V,M>::clearPreComputingPositions()
{
  queso_require_equal_to_msg(m_preComputingPositions.size(), m_originalNewtonSteps.size(), "m_preComputingPositions.size() != m_originalNewtonSteps.size()");

  queso_require_equal_to_msg(m_preComputingPositions.size(), m_originalCovMatrices.size(), "m_preComputingPositions.size() != m_originalCovMatrices.size()");

  BaseTKGroup<V,M>::clearPreComputingPositions();

  // RVs are not deleted in base class because the cov matrices are constant in the case of scaledTK class
  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    if (m_rvs[i]) {
      delete m_rvs[i];
      m_rvs[i] = NULL;
    }
  }

  for (unsigned int i = 0; i < m_originalNewtonSteps.size(); ++i) {
    if (m_originalNewtonSteps[i]) {
      delete m_originalNewtonSteps[i];
      m_originalNewtonSteps[i] = NULL;
    }
  }

  for (unsigned int i = 0; i < m_originalCovMatrices.size(); ++i) {
    if (m_originalCovMatrices[i]) {
      delete m_originalCovMatrices[i];
      m_originalCovMatrices[i] = NULL;
    }
  }

  return;
}

template <class V, class M>
unsigned int
HessianCovMatricesTKGroup<V, M>::set_dr_stage(unsigned int stageId)
{
  unsigned int old_stageId = this->m_stageId;
  this->m_stageId = stageId;
  return old_stageId;
}


// I/O methods---------------------------------------
template<class V, class M>
void
HessianCovMatricesTKGroup<V,M>::print(std::ostream& os) const
{
  BaseTKGroup<V,M>::print(os);
  return;
}

}  // End namespace QUESO

template class QUESO::HessianCovMatricesTKGroup<QUESO::GslVector, QUESO::GslMatrix>;
