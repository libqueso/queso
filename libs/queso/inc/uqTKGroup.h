/* uq/libs/queso/inc/uqTK.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

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
  m_rvs                  (scales.size(),NULL)
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
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    m_rvs[i] = new uqGaussianVectorRVClass<V,M>(m_prefix.c_str(),
                                                *m_vectorSpace,
                                                m_vectorSpace->zeroVector(),
                                                factor*covMatrix);
  }
}

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::~uqScaledCovMatrixTKGroupClass()
{
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
  uqBaseTKGroupClass<V,M>::setPreComputingPosition(position,stageId);
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

  M* tmpHessian = m_vectorSpace->newMatrix();
  V* tmpGrad    = m_vectorSpace->newVector();
  for (unsigned int i = 0; i < m_preComputedPosPlusNewton.size(); ++i) {
    m_targetPdf.minus2LnValue(position,NULL,tmpGrad,tmpHessian,NULL);
    //double factor = 1./m_scales[i]/m_scales[i];
    //*tmpHessian *= factor;

    m_preComputedPosPlusNewton[i] = m_vectorSpace->newVector();
    *(m_preComputedPosPlusNewton[i]) = *m_preComputingPositions[i] - tmpHessian->invertMultiply(*tmpGrad);

    m_rvs[i] = new uqGaussianVectorRVClass<V,M>(m_prefix.c_str(),
                                                *m_vectorSpace,
                                                *m_preComputedPosPlusNewton[i],
                                                *tmpHessian);
  }
  delete tmpGrad;
  delete tmpHessian;

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
