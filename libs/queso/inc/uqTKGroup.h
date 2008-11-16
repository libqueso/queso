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

  virtual       void                    print(std::ostream& os) const = 0;

protected:
  const   uqEmptyEnvironmentClass*                    m_emptyEnv;
  const   uqBaseEnvironmentClass&                     m_env;
          std::string                                 m_prefix;
  const   uqVectorSpaceClass<V,M>*                    m_vectorSpace;
          std::vector<double>                         m_scales;
  const   V*                                          m_mInfVec;
  const   V*                                          m_pInfVec;
          std::vector<uqGaussianVectorRVClass<V,M>* > m_rvs; // Gaussian, not Base... And nothing const...
};

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass()
  :
  m_emptyEnv   (new uqEmptyEnvironmentClass()),
  m_env        (*m_emptyEnv),
  m_prefix     (""),
  m_vectorSpace(NULL),
  m_scales     (0),
  m_mInfVec    (NULL),
  m_pInfVec    (NULL),
  m_rvs        (0)
{
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const std::vector<double>&     scales)
  :
  m_emptyEnv   (NULL),
  m_env        (vectorSpace.env()),
  m_prefix     (prefix),
  m_vectorSpace(&vectorSpace),
  m_scales     (scales.size(),1.),
  m_mInfVec    (vectorSpace.newVector(-INFINITY)),
  m_pInfVec    (vectorSpace.newVector( INFINITY)),
  m_rvs        (0)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    m_scales[i] = scales[i];
  }
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::~uqBaseTKGroupClass()
{
  if (m_pInfVec)  delete m_pInfVec;
  if (m_mInfVec)  delete m_mInfVec;
  if (m_emptyEnv) delete m_emptyEnv;
}

template<class V, class M>
const uqBaseEnvironmentClass&
uqBaseTKGroupClass<V,M>::env() const
{
  return m_env;
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

       void                          updateCovMatrix(const M& covMatrix);
 const uqGaussianVectorRVClass<V,M>& rv             (const V&                     position ) const;
 const uqGaussianVectorRVClass<V,M>& rv             (const std::vector<const V*>& positions) const;
       void                          print          (std::ostream& os) const;

private:
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_vectorSpace;
  using uqBaseTKGroupClass<V,M>::m_scales;
  using uqBaseTKGroupClass<V,M>::m_mInfVec;
  using uqBaseTKGroupClass<V,M>::m_pInfVec;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  M m_covMatrix;
};

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::uqScaledCovMatrixTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
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
                                                *m_mInfVec,
                                                *m_pInfVec,
                                                m_vectorSpace->zeroVector(),
                                                factor*covMatrix);
  }
}

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::~uqScaledCovMatrixTKGroupClass()
{
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
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKGroupClass<V,M>::rv(const V& position) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() == 0,
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() = 0");

  m_rvs[0]->updateExpectedValues(position);

  return (*m_rvs[0]);
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKGroupClass<V,M>::rv(const std::vector<const V*>& positions) const
{
  UQ_FATAL_TEST_MACRO(positions.size() == 0,
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "positions.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs.size() < positions.size(),
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() < positions.size()");

  unsigned int i = positions.size()-1;
  m_rvs[i]->updateExpectedValues(*positions[0]);

  return (*m_rvs[i]);
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::print(std::ostream& os) const
{
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

       void                          addPreComputingPosition   (unsigned int positionId, const V& position);
       void                          clearPreComputingPositions();
 const uqGaussianVectorRVClass<V,M>& rv                        (unsigned int                     positionId ) const;
 const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& positionIds) const;
       void                          print                     (std::ostream& os) const;

private:
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_vectorSpace;
  using uqBaseTKGroupClass<V,M>::m_scales;
  using uqBaseTKGroupClass<V,M>::m_mInfVec;
  using uqBaseTKGroupClass<V,M>::m_pInfVec;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  const uqBaseVectorPdfClass<V,M>& m_targetPdf;
  std::vector<const V*>            m_preComputingPositions;
  std::vector<const V*>            m_preComputedGrads;
//std::vector<const M*>            m_preComputedHessians; // Hessians are stored inside the rvs
};

template<class V, class M>
uqHessianCovMatricesTKGroupClass<V,M>::uqHessianCovMatricesTKGroupClass(
  const char*                      prefix,
  const uqVectorSpaceClass<V,M>&   vectorSpace,
  const std::vector<double>&       scales,
  const uqBaseVectorPdfClass<V,M>& targetPdf)
  :
  uqBaseTKGroupClass<V,M>(prefix,vectorSpace,scales),
  m_targetPdf            (targetPdf),
  m_preComputingPositions(scales.size(),NULL),
  m_preComputedGrads     (scales.size(),NULL)
//m_preComputedHessians  (scales.size(),NULL)
{
}

template<class V, class M>
uqHessianCovMatricesTKGroupClass<V,M>::~uqHessianCovMatricesTKGroupClass()
{
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::addPreComputingPosition(unsigned int positionId, const V& position)
{
  return;
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::clearPreComputingPositions()
{
  return;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqHessianCovMatricesTKGroupClass<V,M>::rv(unsigned int positionId) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() <= positionId,
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() <= positionId");

  return *m_rvs[positionId];
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqHessianCovMatricesTKGroupClass<V,M>::rv(const std::vector<unsigned int>& positionIds) const
{
  UQ_FATAL_TEST_MACRO(positionIds.size() == 0,
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "positionIds.size() = 0");

  UQ_FATAL_TEST_MACRO(m_rvs.size() < positionIds.size(),
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() < positionIds.size()");

  return *m_rvs[positionIds.size()-1];
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_TRANSITION_KERNEL_GROUP_H__
