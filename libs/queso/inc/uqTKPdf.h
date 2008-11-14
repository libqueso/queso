/* uq/libs/queso/inc/uqTKPdf.h
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

#ifndef __UQ_TK_PDF_H__
#define __UQ_TK_PDF_H__

#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseTKPdfClass {
public:
           uqBaseTKPdfClass(const char*                    prefix,
                            const uqVectorSpaceClass<V,M>& domainSpace,
                            const V&                       expVector,
                            const M&                       covMatrix);
  virtual ~uqBaseTKPdfClass();

  const   uqVectorSpaceClass<V,M>& domainSpace    ()                     const;
  virtual double                   actualDensity  (const V& paramValues) const = 0;
  virtual double                   minus2LnDensity(const V& paramValues) const = 0;

protected:
  const   uqBaseEnvironmentClass&  m_env;
          std::string              m_prefix;
  const   uqVectorSpaceClass<V,M>& m_domainSpace;
};

template<class V, class M>
uqBaseTKPdfClass<V,M>::uqBaseTKPdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace)
  :
  m_env        (domainSpace.env()          ),
  m_prefix     ((std::string)(prefix)+"tk_"),
  m_domainSpace(domainSpace                )
{
}

template<class V, class M>
uqBaseTKPdfClass<V,M>::~uqBaseTKPdfClass()
{
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianTKPdfClass : public uqBaseTKPdfClass<V,M> {
public:
  uqGaussianTKPdfClass(const char*                    prefix,
                       const uqVectorSpaceClass<V,M>& domainSpace,
                       const V&                       expVector,
                       const M&                       covMatrix);
 ~uqGaussianTKPdfClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

  void setExpectedVector  (const V& expVector);
  void setCovarianceMatrix(const M& covMatrix);

private:
  using uqBaseTKPdfClass<V,M>::m_env;
  using uqBaseTKPdfClass<V,M>::m_prefix;
  using uqBaseTKPdfClass<V,M>::m_domainSpace;

  V* m_expVector;
  M* m_covMatrix;
};

template<class V, class M>
uqGaussianTKPdfClass<V,M>::uqGaussianTKPdfClass()
  :
  uqBaseTKPdfClass<V,M>(prefix,domainSpace  ),
  m_expVector           (new V(expectedVetor)),
  m_covMatrix          (new M(covMatrix)    )
{
}

template<class V, class M>
double
uqBaseTKPdfClass<V,M>::actualDensity(const V& paramValues) const
{
  return 0.;
}

template<class V, class M>
double
uqBaseTKPdfClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  return 0.;
}

template<class V, class M>
void
uqBaseTKPdfClass<V,M>::setExpectedVector(const V& expVector)
{
  if (m_expVector) delete m_expVector;
  m_expVector = new V(expVector);

  return;
}

template<class V, class M>
void
uqBaseTKPdfClass<V,M>::setCovarianceMatrix(const M& covMatrix)
{
  if (m_covMatrix) delete m_covMatrix;
  m_covMatrix = new M(covMatrix);

  return;
}
#endif // __UQ_TK_PDF_H__
