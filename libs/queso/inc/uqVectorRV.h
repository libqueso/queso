/* uq/libs/queso/inc/uqVectorRV.h
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

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqVectorSpace.h>
#include <uqVectorPdf.h>
#include <uqVectorRealizer.h>
#include <uqVectorCdf.h>
#include <uqVectorMdf.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorRVClass {
public:
  uqBaseVectorRVClass(const char*                    prefix,
                      const uqVectorSpaceClass<V,M>& imageSpace);
  virtual ~uqBaseVectorRVClass();

  const   uqBaseEnvironmentClass&         env        () const;
  const   uqVectorSpaceClass       <V,M>& imageSpace () const;
  const   uqBaseVectorPdfClass     <V,M>& pdf        () const;
  const   uqBaseVectorRealizerClass<V,M>& realizer   () const;
  const   uqBaseVectorCdfClass     <V,M>& cdf        () const;
  const   uqBaseVectorMdfClass     <V,M>& mdf        () const;

  virtual void                            setPdf     (const uqBaseVectorPdfClass     <V,M>& pdf     ) = 0;
  virtual void                            setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer) = 0;
  virtual void                            setCdf     (const uqBaseVectorCdfClass     <V,M>& cdf     ) = 0;
  virtual void                            setMdf     (const uqBaseVectorMdfClass     <V,M>& mdf     ) = 0;

  virtual void                            print      (std::ostream& os) const;

protected:
  const   uqBaseEnvironmentClass&         m_env;
          std::string                     m_prefix;
  const   uqVectorSpaceClass       <V,M>& m_imageSpace;
          uqBaseVectorPdfClass     <V,M>* m_pdf;
	  uqBaseVectorRealizerClass<V,M>* m_realizer;
  const   uqBaseVectorCdfClass     <V,M>* m_cdf;
  const   uqBaseVectorMdfClass     <V,M>* m_mdf;
};

template<class V, class M>
uqBaseVectorRVClass<V,M>::uqBaseVectorRVClass(
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace)
  :
  m_env       (imageSpace.env()),
  m_prefix    ((std::string)(prefix)+"rv_"),
  m_imageSpace(imageSpace),
  m_pdf       (NULL),
  m_realizer  (NULL),
  m_cdf       (NULL),
  m_mdf       (NULL)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRVClass<V,M>::~uqBaseVectorRVClass()
{
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorRVClass<V,M>::imageSpace() const
{
  return m_imageSpace;
}

template<class V, class M>
const uqBaseVectorPdfClass<V,M>&
uqBaseVectorRVClass<V,M>::pdf() const
{
  UQ_FATAL_TEST_MACRO(m_pdf == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::pdf()",
                      "m_pdf is NULL");

  return *m_pdf;
}

template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqBaseVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::realizer(), prefix=")+m_prefix,
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
const uqBaseVectorCdfClass<V,M>&
uqBaseVectorRVClass<V,M>::cdf() const
{
  UQ_FATAL_TEST_MACRO(m_cdf == NULL,
                      m_env.rank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::cdf(), prefix=")+m_prefix,
                      "m_cdf is NULL");

  return *m_cdf;
}

template<class V, class M>
const uqBaseVectorMdfClass<V,M>&
uqBaseVectorRVClass<V,M>::mdf() const
{
  UQ_FATAL_TEST_MACRO(m_mdf == NULL,
                      m_env.rank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::mdf(), prefix=")+m_prefix,
                      "m_mdf is NULL");

  return *m_mdf;
}

template <class V, class M>
const uqBaseEnvironmentClass&
uqBaseVectorRVClass<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::print(std::ostream& os) const
{
#if 0
  os << "\nComponents are:"
     << std::endl;
  for (unsigned int i = 0; i < m_components.size(); ++i) {
    os << i << " ";
    if (m_components[i]) {
      os << *(m_components[i]);
    }
    else {
      os << "NULL";
    }
    os << std::endl;
  }
#endif
  return;
}

template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqBaseVectorRVClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V, class M>
class uqGenericVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGenericVectorRVClass(const char*                           prefix,
                         const uqVectorSpaceClass       <V,M>& imageSpace);
  uqGenericVectorRVClass(const char*                           prefix,
                         const uqVectorSpaceClass       <V,M>& imageSpace,
                         const uqBaseVectorPdfClass     <V,M>& pdf,
                         const uqBaseVectorRealizerClass<V,M>& realizer,
                         const uqBaseVectorCdfClass     <V,M>& cdf,
                         const uqBaseVectorMdfClass     <V,M>& mdf);
  virtual ~uqGenericVectorRVClass();

          void setPdf     (const uqBaseVectorPdfClass     <V,M>& pdf     );
          void setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer);
          void setCdf     (const uqBaseVectorCdfClass     <V,M>& cdf     );
          void setMdf     (const uqBaseVectorMdfClass     <V,M>& mdf     );

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_cdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                     prefix,
  const uqVectorSpaceClass <V,M>& imageSpace)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGenericVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGenericVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                           prefix,
  const uqVectorSpaceClass       <V,M>& imageSpace,
  const uqBaseVectorPdfClass     <V,M>& pdf,
  const uqBaseVectorRealizerClass<V,M>& realizer,
  const uqBaseVectorCdfClass     <V,M>& cdf,
  const uqBaseVectorMdfClass     <V,M>& mdf)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGenericVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_pdf      = &pdf;
  m_realizer = &realizer;
  m_cdf      = &cdf;
  m_mdf      = &mdf;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGenericVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRVClass<V,M>::~uqGenericVectorRVClass()
{
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setPdf(const uqBaseVectorPdfClass<V,M>& pdf)
{
  m_pdf = &pdf;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setCdf(const uqBaseVectorCdfClass<V,M>& cdf)
{
  m_cdf = &cdf;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setMdf(const uqBaseVectorMdfClass<V,M>& mdf)
{
  m_mdf = &mdf;
  return;
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGaussianVectorRVClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const V&                       imageMinValues,
                          const V&                       imageMaxValues,
                          const V&                       imageExpectedValues,
                          const V&                       imageVarianceValues);
  uqGaussianVectorRVClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const V&                       imageMinValues,
                          const V&                       imageMaxValues,
                          const V&                       imageExpectedValues,
                          const M&                       covMatrix);
  virtual ~uqGaussianVectorRVClass();

          void setPdf     (const uqBaseVectorPdfClass     <V,M>& pdf     );
          void setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer);
          void setCdf     (const uqBaseVectorCdfClass     <V,M>& cdf     );
          void setMdf     (const uqBaseVectorMdfClass     <V,M>& mdf     );
	  void updateExpectedValues(const V& newExpectedValues );

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_cdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectedValues,
  const V&                       imageVarianceValues)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_pdf = new uqGaussianVectorPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSpace,
                                            imageMinValues,
                                            imageMaxValues,
                                            imageExpectedValues,
                                            imageVarianceValues);

  M lowerCholCovMatrix(imageVarianceValues);
  int ierr = lowerCholCovMatrix.chol();
  UQ_FATAL_TEST_MACRO( (ierr!=0), m_env.rank(),
		       "uqGaussianVectorRVClass<V,M>::constructor() [2]",
		       "Cholesky decomposition of covariance matrix failed.");
  lowerCholCovMatrix.zeroUpper(false);

  m_cdf         = NULL; // FIX ME: complete code
  m_mdf         = NULL; // FIX ME: complete code

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_pdf = new uqGaussianVectorPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSpace,
                                            imageMinValues,
                                            imageMaxValues,
                                            imageExpectedValues,
                                            covMatrix);

  M lowerCholCovMatrix(covMatrix);
  int ierr = lowerCholCovMatrix.chol();
  UQ_FATAL_TEST_MACRO( (ierr!=0), m_env.rank(),
		       "uqGaussianVectorRVClass<V,M>::constructor() [2]",
		       "Cholesky decomposition of covariance matrix failed.");
  lowerCholCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
						      m_imageSpace,
						      imageExpectedValues,
						      lowerCholCovMatrix);

  m_cdf         = NULL; // FIX ME: complete code
  m_mdf         = NULL; // FIX ME: complete code

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::~uqGaussianVectorRVClass()
{
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setPdf(const uqBaseVectorPdfClass<V,M>& pdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setPdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setRealizer()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setCdf(const uqBaseVectorCdfClass<V,M>& cdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setCdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setMdf(const uqBaseVectorMdfClass<V,M>& mdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setMdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateExpectedValues(const V& newExpectedValues)
{
  // we are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianVectorPdfClass     <V,M>* >(m_pdf     ) )->updateExpectedValues(newExpectedValues);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateExpectedValues(newExpectedValues);
  return;
}

//*****************************************************
// Uniform class
//*****************************************************
template<class V, class M>
class uqUniformVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqUniformVectorRVClass(const char*                     prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const V&                       imageMinValues,
                          const V&                       imageMaxValues);
  virtual ~uqUniformVectorRVClass();

          void setPdf     (const uqBaseVectorPdfClass     <V,M>& pdf     );
          void setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer);
          void setCdf     (const uqBaseVectorCdfClass     <V,M>& cdf     );
          void setMdf     (const uqBaseVectorMdfClass     <V,M>& mdf     );

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_cdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqUniformVectorRVClass<V,M>::uqUniformVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const V&                       imageMinValues,
  const V&                       imageMaxValues)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqUniformVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_pdf = new uqUniformVectorPdfClass<V,M>(m_prefix.c_str(),
                                           m_imageSpace,
                                           imageMinValues,
                                           imageMaxValues);
  m_realizer    = NULL; // FIX ME: complete code
  m_cdf         = NULL; // FIX ME: complete code
  m_mdf         = NULL; // FIX ME: complete code

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqUniformVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqUniformVectorRVClass<V,M>::~uqUniformVectorRVClass()
{
}

template<class V, class M>
void
uqUniformVectorRVClass<V,M>::setPdf(const uqBaseVectorPdfClass<V,M>& pdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqUniformVectorRVClass<V,M>::setPdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqUniformVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqUniformVectorRVClass<V,M>::setRealizer()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqUniformVectorRVClass<V,M>::setCdf(const uqBaseVectorCdfClass<V,M>& cdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqUniformVectorRVClass<V,M>::setCdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqUniformVectorRVClass<V,M>::setMdf(const uqBaseVectorMdfClass<V,M>& mdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqUniformVectorRVClass<V,M>::setMdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

#endif // __UQ_VECTOR_RV_H__
