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

#ifndef __UQ_TRANSITION_KERNEL_H__
#define __UQ_TRANSITION_KERNEL_H__

#include <uqVectorRV.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseTKClass {
public:
  uqBaseTKClass();
  uqBaseTKClass(const char*                    prefix,
                const uqVectorSpaceClass<V,M>& vectorSpace);
  virtual ~uqBaseTKClass();

          const uqBaseEnvironmentClass& env() const;

  virtual       void                    print(std::ostream& os) const = 0;

protected:
  const   uqEmptyEnvironmentClass* m_emptyEnv;
  const   uqBaseEnvironmentClass&  m_env;
          std::string              m_prefix;

  const   std::vector<const uqGaussianVectorRVClass<V,M>* > m_rvs; // Gaussian, not Base...
};

template<class V, class M>
uqBaseTKClass<V,M>::uqBaseTKClass()
  :
  m_emptyEnv(new uqEmptyEnvironmentClass()),
  m_env     (*m_emptyEnv),
  m_prefix  (""),
  m_rvs     (0)
{
}

template<class V, class M>
uqBaseTKClass<V,M>::uqBaseTKClass(
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
uqBaseTKClass<V,M>::~uqBaseTKClass()
{
  if (m_emptyEnv) delete m_emptyEnv;
}

template<class V, class M>
const uqBaseEnvironmentClass&
uqBaseTKClass<V,M>::env() const
{
  return m_env;
}

//*****************************************************
// TK with scaled cov matrix
//*****************************************************
template<class V, class M>
class uqScaledCovMatrixTKClass : public uqBaseTKClass<V,M> {
public:
  uqScaledCovMatrixTKClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& vectorSpace,
                           const M&                       covMatrix,
                           const std::vector<double>&     scales);
 ~uqScaledCovMatrixTKClass();

       void                          updateCovMatrix(const M& covMatrix);
 const uqGaussianVectorRVClass<V,M>& rv             (const std::vector<const V*>& positions) const;
       void                          print          (std::ostream& os) const;

private:
  using uqBaseTKClass<V,M>::m_env;
  using uqBaseTKClass<V,M>::m_prefix;
  using uqBaseTKClass<V,M>::m_rvs;

  M                   m_covMatrix;
  std::vector<double> m_scales;
};

template<class V, class M>
uqScaledCovMatrixTKClass<V,M>::uqScaledCovMatrixTKClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const M&                       covMatrix,
  const std::vector<double>&     scales)
  :
  uqBaseTKClass<V,M>(prefix,vectorSpace),
  m_covMatrix       (covMatrix),
  m_scales          (0)
{
}

template<class V, class M>
uqScaledCovMatrixTKClass<V,M>::~uqScaledCovMatrixTKClass()
{
}

template<class V, class M>
void
uqScaledCovMatrixTKClass<V,M>::updateCovMatrix(const M& covMatrix)
{
  return;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqScaledCovMatrixTKClass<V,M>::rv(const std::vector<const V*>& positions) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < positions.size(),
                      m_env.rank(),
                      "uqScaledCovTKClass<V,M>::pdf()",
                      "m_rvs.size() < positions.size()");

  return m_rvs[positions.size()-1];
}

template<class V, class M>
void
uqScaledCovMatrixTKClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// TK with Hessians
//*****************************************************
template<class V, class M>
class uqHessianCovMatricesTKClass : public uqBaseTKClass<V,M> {
public:
  uqHessianCovMatricesTKClass(const char*                      prefix,
                              const uqVectorSpaceClass<V,M>&   vectorSpace,
                              const uqBaseVectorPdfClass<V,M>& targetPdf,
                              unsigned int                     maxNumberOfStages);
 ~uqHessianCovMatricesTKClass();

       void                          addPreComputingPosition   (unsigned int positionId, const V& position);
       void                          clearPreComputingPositions();
 const uqGaussianVectorRVClass<V,M>& rv                        (const std::vector<unsigned int>& positionIds) const;
       void                          print                     (std::ostream& os) const;

private:
  using uqBaseTKClass<V,M>::m_env;
  using uqBaseTKClass<V,M>::m_prefix;
  using uqBaseTKClass<V,M>::m_rvs;

  const uqBaseVectorPdfClass<V,M>& m_targetPdf;
  unsigned int                     m_maxNumberOfStages;
  std::vector<const V*>            m_preComputingPositions;
  std::vector<const V*>            m_preComputedGrads;
  //std::vector<const M*>            m_preComputedHessians; // Hessians are saved inside the rvs
};

template<class V, class M>
uqHessianCovMatricesTKClass<V,M>::uqHessianCovMatricesTKClass(
  const char*                      prefix,
  const uqVectorSpaceClass<V,M>&   vectorSpace,
  const uqBaseVectorPdfClass<V,M>& targetPdf,
  unsigned int                     maxNumberOfStages)
  :
  uqBaseTKClass<V,M>     (prefix,vectorSpace),
  m_targetPdf            (targetPdf),
  m_maxNumberOfStages    (maxNumberOfStages),
  m_preComputingPositions(m_maxNumberOfStages,NULL),
  m_preComputedGrads     (m_maxNumberOfStages,NULL)
//m_preComputedHessians  (m_maxNumberOfStages,NULL)
{
}

template<class V, class M>
uqHessianCovMatricesTKClass<V,M>::~uqHessianCovMatricesTKClass()
{
}

template<class V, class M>
void
uqHessianCovMatricesTKClass<V,M>::addPreComputingPosition(unsigned int positionId, const V& position)
{
  return;
}

template<class V, class M>
void
uqHessianCovMatricesTKClass<V,M>::clearPreComputingPositions()
{
  return;
}

template<class V, class M>
const uqGaussianVectorRVClass<V,M>&
uqHessianCovMatricesTKClass<V,M>::rv(const std::vector<unsigned int>& positionIds) const
{
  UQ_FATAL_TEST_MACRO(m_rvs.size() < positionIds.size(),
                      m_env.rank(),
                      "uqScaledCovTKClass<V,M>::pdf()",
                      "m_rvs.size() < positionIds.size()");

  return m_rvs[positionIds.size()-1];
}

template<class V, class M>
void
uqHessianCovMatricesTKClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_TRANSITION_KERNEL_H__
