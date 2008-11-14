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
                     const uqVectorSpaceClass<V,M>& vectorSpace);
  virtual ~uqBaseTKGroupClass();

          const uqBaseEnvironmentClass& env() const;

  virtual       void                    print(std::ostream& os) const = 0;

protected:
  const   uqEmptyEnvironmentClass* m_emptyEnv;
  const   uqBaseEnvironmentClass&  m_env;
          std::string              m_prefix;

  const   std::vector<const uqGaussianVectorRVClass<V,M>* > m_rvs; // Gaussian, not Base...
};

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass()
  :
  m_emptyEnv(new uqEmptyEnvironmentClass()),
  m_env     (*m_emptyEnv),
  m_prefix  (""),
  m_rvs     (0)
{
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::uqBaseTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace)
  :
  m_emptyEnv(NULL),
  m_env     (vectorSpace.env()),
  m_prefix  (prefix),
  m_rvs     (0)
{
}

template<class V, class M>
uqBaseTKGroupClass<V,M>::~uqBaseTKGroupClass()
{
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
                                const M&                       covMatrix,
                                const std::vector<double>&     scales);
 ~uqScaledCovMatrixTKGroupClass();

       void                          updateCovMatrix(const M& covMatrix);
 const uqGaussianVectorRVClass<V,M>& rv             (const std::vector<const V*>& positions) const;
       void                          print          (std::ostream& os) const;

private:
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  M                   m_covMatrix;
  std::vector<double> m_scales;
};

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::uqScaledCovMatrixTKGroupClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const M&                       covMatrix,
  const std::vector<double>&     scales)
  :
  uqBaseTKGroupClass<V,M>(prefix,vectorSpace),
  m_covMatrix       (covMatrix),
  m_scales          (0)
{
}

template<class V, class M>
uqScaledCovMatrixTKGroupClass<V,M>::~uqScaledCovMatrixTKGroupClass()
{
}

template<class V, class M>
void
uqScaledCovMatrixTKGroupClass<V,M>::updateCovMatrix(const M& covMatrix)
{
  return;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKGroupClass<V,M>::rv(const std::vector<const V*>& positions) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < positions.size(),
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() < positions.size()");

  return (*m_rvs[positions.size()-1]);
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
                                   const uqBaseVectorPdfClass<V,M>& targetPdf,
                                   unsigned int                     maxNumberOfStages);
 ~uqHessianCovMatricesTKGroupClass();

       void                          addPreComputingPosition   (unsigned int positionId, const V& position);
       void                          clearPreComputingPositions();
 const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& positionIds) const;
       void                          print                     (std::ostream& os) const;

private:
  using uqBaseTKGroupClass<V,M>::m_env;
  using uqBaseTKGroupClass<V,M>::m_prefix;
  using uqBaseTKGroupClass<V,M>::m_rvs;

  const uqBaseVectorPdfClass<V,M>& m_targetPdf;
  unsigned int                     m_maxNumberOfStages;
  std::vector<const V*>            m_preComputingPositions;
  std::vector<const V*>            m_preComputedGrads;
  //std::vector<const M*>            m_preComputedHessians; // Hessians are stored inside the rvs
};

template<class V, class M>
uqHessianCovMatricesTKGroupClass<V,M>::uqHessianCovMatricesTKGroupClass(
  const char*                      prefix,
  const uqVectorSpaceClass<V,M>&   vectorSpace,
  const uqBaseVectorPdfClass<V,M>& targetPdf,
  unsigned int                     maxNumberOfStages)
  :
  uqBaseTKGroupClass<V,M>     (prefix,vectorSpace),
  m_targetPdf            (targetPdf),
  m_maxNumberOfStages    (maxNumberOfStages),
  m_preComputingPositions(m_maxNumberOfStages,NULL),
  m_preComputedGrads     (m_maxNumberOfStages,NULL)
//m_preComputedHessians  (m_maxNumberOfStages,NULL)
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
uqHessianCovMatricesTKGroupClass<V,M>::rv(const std::vector<unsigned int>& positionIds) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < positionIds.size(),
                      m_env.rank(),
                      "uqScaledCovTKGroupClass<V,M>::pdf()",
                      "m_rvs.size() < positionIds.size()");

  return m_rvs[positionIds.size()-1];
}

template<class V, class M>
void
uqHessianCovMatricesTKGroupClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_TRANSITION_KERNEL_GROUP_H__
