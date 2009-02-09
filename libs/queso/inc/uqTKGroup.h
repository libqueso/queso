/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_TRANSITION_KERNEL_GROUP_H__
#define __UQ_TRANSITION_KERNEL_GROUP_H__

#include <uqVectorRV.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseTKGroupClass {
public:
  uqBaseTKGroupClass();
  uqBaseTKGroupClass(const char*                    prefix,
                     const uqVectorSpaceClass<V,M>& vectorSpace,
                     const std::vector<double>&     scales);
  virtual ~uqBaseTKGroupClass();

          const uqBaseEnvironmentClass& env() const;

  virtual       bool                          symmetric                 () const = 0;
  virtual const uqGaussianVectorRVClass<V,M>& rv                        (unsigned int                     stageId ) = 0;
  virtual const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& stageIds) = 0;
          const V&                            preComputingPosition      (unsigned int stageId) const;
  virtual       void                          setPreComputingPosition   (const V& position, unsigned int stageId);
  virtual       void                          clearPreComputingPositions();
  virtual       void                          print                     (std::ostream& os) const;

protected:
  const   uqEmptyEnvironmentClass*                    m_emptyEnv;
  const   uqBaseEnvironmentClass&                     m_env;
          std::string                                 m_prefix;
  const   uqVectorSpaceClass<V,M>*                    m_vectorSpace;
          std::vector<double>                         m_scales;
          std::vector<const V*>                       m_preComputingPositions;
          std::vector<uqGaussianVectorRVClass<V,M>* > m_rvs; // Gaussian, not Base... And nothing const...
};

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass()
  :
  m_emptyEnv             (new uqEmptyEnvironmentClass()),
  m_env                  (*m_emptyEnv),
  m_prefix               (""),
  m_vectorSpace          (NULL),
  m_scales               (0),
  m_preComputingPositions(NULL),
  m_rvs                  (0)
{
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const std::vector<double>&     scales)
  :
  m_emptyEnv             (NULL),
  m_env                  (vectorSpace.env()),
  m_prefix               (prefix),
  m_vectorSpace          (&vectorSpace),
  m_scales               (scales.size(),1.),
  m_preComputingPositions(scales.size()+1,NULL), // Yes, +1
  m_rvs                  (scales.size(),NULL) // IMPORTANT: it stays like this for scaledTK, but it will be overwritten to '+1' by hessianTK constructor
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    m_scales[i] = scales[i];
  }
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::~uqBaseTKGroupClass()
{
  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    if (m_rvs[i]) delete m_rvs[i];
  }
  for (unsigned int i = 0; i < m_preComputingPositions.size(); ++i) {
    if (m_preComputingPositions[i]) delete m_preComputingPositions[i];
  }
  if (m_emptyEnv) delete m_emptyEnv;
}

template<class V, class M>
const uqBaseEnvironmentClass&
uqBaseTKGroupClass<V,M>::env() const
{
  return m_env;
}

template<class V, class M>
const V&
uqBaseTKGroupClass<V,M>::preComputingPosition(unsigned int stageId) const
{
  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.rank(),
                      "uqBaseTKGroupClass<V,M>::preComputingPosition()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] == NULL,
                      m_env.rank(),
                      "uqBaseTKGroupClass<V,M>::preComputingPosition()",
                      "m_preComputingPositions[stageId] == NULL");

  return *m_preComputingPositions[stageId];
}

template<class V, class M>
void
uqBaseTKGroupClass<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.rank(),
                      "uqBaseTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] != NULL,
                      m_env.rank(),
                      "uqBaseTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions[stageId] != NULL");

  m_preComputingPositions[stageId] = new V(position);

  return;
}

template<class V, class M>
void
uqBaseTKGroupClass<V,M>::clearPreComputingPositions()
{
  for (unsigned int i = 0; i < m_preComputingPositions.size(); ++i) {
    if (m_preComputingPositions[i]) {
      delete m_preComputingPositions[i];
      m_preComputingPositions[i] = NULL;
    }
  }

  return;
}

template<class V, class M>
void
uqBaseTKGroupClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// TK with scaled cov matrix
//*****************************************************
template<class V, class M>
class uqScaledCovMatrixTKGroupClass : public uqBaseTKGroupClass<V,M> {
public:
  uqScaledCovMatrixTKGroupClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& vectorSpace,
                                const std::vector<double>&     scales,
                                const M&                       covMatrix);
 ~uqScaledCovMatrixTKGroupClass();

        bool                          symmetric                 () const;
  const uqGaussianVectorRVClass<V,M>& rv                        (unsigned int                     stageId );
  const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& stageIds);
        void                          setPreComputingPosition   (const V& position, unsigned int stageId);
        void                          clearPreComputingPositions();
        void                          updateCovMatrix           (const M& covMatrix);
        void                          print                     (std::ostream& os) const;

private:
        void                          setRVsWithZeroMean        ();
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_vectorSpace;
  using uqBaseTKGroupClass<V,M>::m_scales;
  using uqBaseTKGroupClass<V,M>::m_preComputingPositions;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  M m_covMatrix;
};

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::uqScaledCovMatrixTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace, // FIX ME: vectorSubset ???
  const std::vector<double>&     scales,
  const M&                       covMatrix)
  :
  uqBaseTKGroupClass<V,M>(prefix,vectorSpace,scales),
  m_covMatrix            (covMatrix)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqScaledCovMatrixTKGroupClass<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqScaledCovMatrixTKGroupClass<V,M>::constructor()"
              << ": m_scales.size() = "                << m_scales.size()
              << ", m_preComputingPositions.size() = " << m_preComputingPositions.size()
              << ", m_rvs.size() = "                   << m_rvs.size()
              << ", m_covMatrix = "                    << m_covMatrix
              << std::endl;
  }

  setRVsWithZeroMean();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqScaledCovMatrixTKGroupClass<V,M>::constructor()"
              << std::endl;
  }
}

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::~uqScaledCovMatrixTKGroupClass()
{
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::setRVsWithZeroMean()
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() == 0,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::setRVsWithZeroMean()",
                      "m_rvs.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs.size() != m_scales.size(),
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::setRVsWithZeroMean()",
                      "m_rvs.size() != m_scales.size()");

  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    UQ_FATAL_TEST_MACRO(m_rvs[i] != NULL,
                        m_env.rank(),
                        "uqScaledCovMatrixTKGroupClass<V,M>::setRVsWithZeroMean()",
                        "m_rvs[i] != NULL");
    m_rvs[i] = new uqGaussianVectorRVClass<V,M>(m_prefix.c_str(),
                                                *m_vectorSpace,
                                                m_vectorSpace->zeroVector(),
                                                factor*m_covMatrix);
  }

  return;
}

template<class V, class M>
bool
uqScaledCovMatrixTKGroupClass<V,M>::symmetric() const
{
  return true;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKGroupClass<V,M>::rv(unsigned int stageId)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() == 0,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv1()",
                      "m_rvs.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs[0] == NULL, // Yes, '0', because that is the id used below
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv1()",
                      "m_rvs[0] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv1()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] == NULL,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv1()",
                      "m_preComputingPositions[stageId] == NULL");

  m_rvs[0]->updateExpVector(*m_preComputingPositions[stageId]);

  return (*m_rvs[0]);
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKGroupClass<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < stageIds.size(),
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv2()",
                      "m_rvs.size() < stageIds.size()");

  UQ_FATAL_TEST_MACRO(m_rvs[stageIds.size()-1] == NULL,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv2()",
                      "m_rvs[stageIds.size()-1] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageIds[0],
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv2()",
                      "m_preComputingPositions.size() <= stageIds[0]");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageIds[0]] == NULL,
                      m_env.rank(),
                      "uqScaledCovMatrixTKGroupClass<V,M>::rv2()",
                      "m_preComputingPositions[stageIds[0]] == NULL");

  m_rvs[stageIds.size()-1]->updateExpVector(*m_preComputingPositions[stageIds[0]]);

  return (*m_rvs[stageIds.size()-1]);
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqScaledCovMatrixTKGroupClass<V,M>::setPreComputingPosition()"
              << ": position = " << position
              << ", stageId = "  << stageId
              << std::endl;
  }

  uqBaseTKGroupClass<V,M>::setPreComputingPosition(position,stageId);
  //setRVsWithZeroMean();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqScaledCovMatrixTKGroupClass<V,M>::setPreComputingPosition()"
              << ", position = "        << position
              << ", stageId = "         << stageId
              << ": preComputingPos = " << *m_preComputingPositions[stageId];
    if (stageId < m_rvs.size()) {
      const uqGaussianVectorPdfClass<V,M>* pdfPtr = dynamic_cast< const uqGaussianVectorPdfClass<V,M>* >(&(m_rvs[stageId]->pdf()));
      std::cout << ", factor = " << 1./m_scales[stageId]/m_scales[stageId]
                << ", rvCov = "  << pdfPtr->covMatrix();
    }
    std::cout << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqScaledCovMatrixTKGroupClass<V,M>::setPreComputingPosition()"
              << ": position = " << position
              << ", stageId = "  << stageId
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::clearPreComputingPositions()
{
  uqBaseTKGroupClass<V,M>::clearPreComputingPositions();
  return;
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::updateCovMatrix(const M& covMatrix)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    m_rvs[i]->updateCovMatrix(factor*covMatrix);
  }

  return;
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::print(std::ostream& os) const
{
  uqBaseTKGroupClass<V,M>::print(os);
  return;
}

//*****************************************************
// TK with Hessians
//*****************************************************
template<class V, class M>
class uqHessianCovMatricesTKGroupClass : public uqBaseTKGroupClass<V,M> {
public:
  uqHessianCovMatricesTKGroupClass(const char*                      prefix,
                                   const uqVectorSpaceClass<V,M>&   vectorSpace,
                                   const std::vector<double>&       scales,
                                   const uqBaseVectorPdfClass<V,M>& targetPdf);
 ~uqHessianCovMatricesTKGroupClass();

        bool                          symmetric                 () const;
  const uqGaussianVectorRVClass<V,M>& rv                        (unsigned int                     stageId );
  const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& stageIds);
        void                          setPreComputingPosition   (const V& position, unsigned int stageId);
        void                          clearPreComputingPositions();
        void                          print                     (std::ostream& os) const;

private:
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_vectorSpace;
  using uqBaseTKGroupClass<V,M>::m_scales;
  using uqBaseTKGroupClass<V,M>::m_preComputingPositions;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  const uqBaseVectorPdfClass<V,M>& m_targetPdf;
  std::vector<V*>                  m_preComputedPosPlusNewton;
//std::vector<const M*>            m_preComputedHessians; // Hessians are stored inside the rvs
};

template<class V, class M>
uqHessianCovMatricesTKGroupClass<V,M>::uqHessianCovMatricesTKGroupClass(
  const char*                      prefix,
  const uqVectorSpaceClass<V,M>&   vectorSpace,
  const std::vector<double>&       scales,
  const uqBaseVectorPdfClass<V,M>& targetPdf)
  :
  uqBaseTKGroupClass<V,M>   (prefix,vectorSpace,scales),
  m_targetPdf               (targetPdf),
  m_preComputedPosPlusNewton(scales.size()+1,NULL) // Yes, +1
//m_preComputedHessians     (scales.size()+1,NULL) // Yes, +1
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqHessianCovMatricesTKGroupClass<V,M>::constructor()"
              << std::endl;
  }

  m_rvs.resize(scales.size()+1,NULL); // Yes, +1 (IMPORTANT: overwrite initialization done by uqBaseTKGroupClass<V,M>)

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::constructor()"
              << ": m_scales.size() = "                   << m_scales.size()
              << ", m_preComputingPositions.size() = "    << m_preComputingPositions.size()
              << ", m_rvs.size() = "                      << m_rvs.size()
              << ", m_preComputedPosPlusNewton.size() = " << m_preComputedPosPlusNewton.size()
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqHessianCovMatricesTKGroupClass<V,M>::constructor()"
              << std::endl;
  }
}

template<class V, class M>
uqHessianCovMatricesTKGroupClass<V,M>::~uqHessianCovMatricesTKGroupClass()
{
}

template<class V, class M>
bool
uqHessianCovMatricesTKGroupClass<V,M>::symmetric() const
{
  return false;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqHessianCovMatricesTKGroupClass<V,M>::rv(unsigned int stageId)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() <= stageId,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv1()",
                      "m_rvs.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_rvs[stageId] == NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv1()",
                      "m_rvs[stageId] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv1()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] == NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv1()",
                      "m_preComputingPositions[stageId] == NULL");

  m_rvs[stageId]->updateExpVector(*m_preComputedPosPlusNewton[stageId]);

  return *m_rvs[stageId];
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqHessianCovMatricesTKGroupClass<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() <= stageIds[0],
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv2()",
                      "m_rvs.size() <= stageIds[0]");

  UQ_FATAL_TEST_MACRO(m_rvs[stageIds[0]] == NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv2()",
                      "m_rvs[stageIds[0]] == NULL");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageIds[0],
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv2()",
                      "m_preComputingPositions.size() <= stageIds[0]");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageIds[0]] == NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::rv2()",
                      "m_preComputingPositions[stageIds[0]] == NULL");

  m_rvs[stageIds[0]]->updateExpVector(*m_preComputedPosPlusNewton[stageIds[0]]);

  return *m_rvs[stageIds[0]];
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()"
              << ": position = " << position
              << ", stageId = "  << stageId
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() != m_preComputedPosPlusNewton.size(),
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions.size() != m_preComputedPosPlusNewton.size()");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() <= stageId,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] != NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions[stageId] != NULL");

  UQ_FATAL_TEST_MACRO(m_rvs.size() <= stageId,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_rvs.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_rvs[stageId] != NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_rvs[stageId] != NULL");

  UQ_FATAL_TEST_MACRO(m_preComputedPosPlusNewton.size() <= stageId,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputedPosPlusNewton.size() <= stageId");

  UQ_FATAL_TEST_MACRO(m_preComputedPosPlusNewton[stageId] != NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputedPosPlusNewton[stageId] != NULL");

  uqBaseTKGroupClass<V,M>::setPreComputingPosition(position,stageId);

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()"
              << ", position = "                          << position
              << ", stageId = "                           << stageId
              << ": m_preComputedPosPlusNewton.size() = " << m_preComputedPosPlusNewton.size()
              << ", m_preComputingPositions.size() = "    << m_preComputingPositions.size()
              << ", m_rvs.size() = "                      << m_rvs.size()
              << std::endl;
  }

  M* tmpHessian = m_vectorSpace->newMatrix();
  M* tmpCovMat  = m_vectorSpace->newMatrix();
  V* tmpGrad    = m_vectorSpace->newVector();
#if 1
  m_targetPdf.minus2LnValue(position,
                            NULL,
                            tmpGrad,
                            tmpHessian,
                            NULL);

  // IMPORTANT: covariance matrix = (Hessian)^{-1} !!!
  V unitVector(m_vectorSpace->zeroVector());
  V multVector(m_vectorSpace->zeroVector());
  for (unsigned int j = 0; j < tmpHessian->numCols(); ++j) {
    if (j > 0) unitVector[j-1] = 0.;
    unitVector[j] = 1.;
    tmpHessian->invertMultiply(unitVector, multVector);
    for (unsigned int i = 0; i < tmpHessian->numRows(); ++i) {
      (*tmpCovMat)(i,j) = multVector[i];
    }
  }
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()"
              << ", position = "  << position
              << ", stageId = "   << stageId
              << ":\n H^{-1} = "  << (*tmpCovMat)
              << "\n H*H^{-1} = " << (*tmpHessian)*(*tmpCovMat)
              << "\n H^{-1}*H = " << (*tmpCovMat)*(*tmpHessian)
              << std::endl;
  }

  double factor = 1./m_scales[stageId]/m_scales[stageId];
  //*tmpCovMat *= factor;

  m_preComputedPosPlusNewton[stageId] = m_vectorSpace->newVector();

  // Check existence of vector to be used
  UQ_FATAL_TEST_MACRO(m_preComputingPositions[stageId] == NULL,
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()",
                      "m_preComputingPositions[stageId] == NULL");
  *(m_preComputedPosPlusNewton[stageId]) = *m_preComputingPositions[stageId] - (*tmpCovMat)*(*tmpGrad); // or tmpHessian->invertMultiply(*tmpGrad);

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()"
              << ", position = "        << position
              << ", stageId = "         << stageId
              << ", about to instantiate a Gaussian rv"
              << ": tmpHessian = "      << *tmpHessian
              << ", preComputingPos = " << *m_preComputingPositions[stageId]
              << ", tmpCovMat = "       << *tmpCovMat
              << ", tmpGrad = "         << *tmpGrad
              << ", preComputedPos = "  << *(m_preComputedPosPlusNewton[stageId])
              << ", factor = "          << factor
              << ", rvCov = "           << factor*(*tmpCovMat)
              << std::endl;
  }
  m_rvs[stageId] = new uqGaussianVectorRVClass<V,M>(m_prefix.c_str(),
                                                    *m_vectorSpace,
                                                    *m_preComputedPosPlusNewton[stageId],
                                                    factor*(*tmpCovMat));
#else
  for (unsigned int i = 0; i < m_preComputedPosPlusNewton.size(); ++i) {
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition(): pos 002, i = " << i << std::endl;
    m_targetPdf.minus2LnValue(position,NULL,tmpGrad,tmpCovMat,NULL);
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition(): pos 003, i = " << i << std::endl;
    //double factor = 1./m_scales[i]/m_scales[i];
    //*tmpCovMat *= factor;

    m_preComputedPosPlusNewton[i] = m_vectorSpace->newVector();
    *(m_preComputedPosPlusNewton[i]) = *m_preComputingPositions[i] - tmpCovMat->invertMultiply(*tmpGrad);

    m_rvs[i] = new uqGaussianVectorRVClass<V,M>(m_prefix.c_str(),
                                                *m_vectorSpace,
                                                *m_preComputedPosPlusNewton[i],
                                                *tmpCovMat);
    std::cout << "In uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition(): pos 004, i = " << i << std::endl;
  }
#endif
  delete tmpGrad;
  delete tmpCovMat;
  delete tmpHessian;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqHessianCovMatricesTKGroupClass<V,M>::setPreComputingPosition()"
              << ": position = " << position
              << ", stageId = "  << stageId
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::clearPreComputingPositions()
{
  UQ_FATAL_TEST_MACRO(m_preComputingPositions.size() != m_preComputedPosPlusNewton.size(),
                      m_env.rank(),
                      "uqHessianCovMatricesTKGroupClass<V,M>::clearPreComputingPositions()",
                      "m_preComputingPositions.size() != m_preComputedPosPlusNewton.size()");

  uqBaseTKGroupClass<V,M>::clearPreComputingPositions();

  // RVs are not deleted in base class because the cov matrices are constant in the case of scaledTK class
  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    if (m_rvs[i]) {
      delete m_rvs[i];
      m_rvs[i] = NULL;
    }
  }

  for (unsigned int i = 0; i < m_preComputedPosPlusNewton.size(); ++i) {
    if (m_preComputedPosPlusNewton[i]) {
      delete m_preComputedPosPlusNewton[i];
      m_preComputedPosPlusNewton[i] = NULL;
    }
  }
  return;
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::print(std::ostream& os) const
{
  uqBaseTKGroupClass<V,M>::print(os);
  return;
}
#endif // __UQ_TRANSITION_KERNEL_GROUP_H__
