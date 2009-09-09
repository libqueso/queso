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

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqVectorSpace.h>
#include <uqJointPdf.h>
#include <uqVectorRealizer.h>
#include <uqVectorCdf.h>
#include <uqVectorMdf.h>
#include <uqSequenceOfVectors.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorRVClass {
public:
  uqBaseVectorRVClass(const char*                  prefix,
                      const uqVectorSetClass<V,M>& imageSet);
  virtual ~uqBaseVectorRVClass();

  const   uqBaseEnvironmentClass&         env       () const;
  const   uqVectorSetClass         <V,M>& imageSet  () const;
  const   uqBaseJointPdfClass      <V,M>& pdf       () const;
  const   uqBaseVectorRealizerClass<V,M>& realizer  () const;
  const   uqBaseVectorCdfClass     <V,M>& subCdf    () const;
  const   uqBaseVectorCdfClass     <V,M>& unifiedCdf() const;
  const   uqBaseVectorMdfClass     <V,M>& mdf       () const;

  virtual void                            print     (std::ostream& os) const = 0;

protected:
  const   uqBaseEnvironmentClass&         m_env;
          std::string                     m_prefix;
  const   uqVectorSetClass         <V,M>& m_imageSet;
          uqBaseJointPdfClass      <V,M>* m_pdf;
	  uqBaseVectorRealizerClass<V,M>* m_realizer;
  const   uqBaseVectorCdfClass     <V,M>* m_subCdf;
  const   uqBaseVectorCdfClass     <V,M>* m_unifiedCdf;
  const   uqBaseVectorMdfClass     <V,M>* m_mdf;
};

template<class V, class M>
uqBaseVectorRVClass<V,M>::uqBaseVectorRVClass(
  const char*                            prefix,
  const uqVectorSetClass          <V,M>& imageSet)
  :
  m_env       (imageSet.env()),
  m_prefix    ((std::string)(prefix)+"rv_"),
  m_imageSet  (imageSet),
  m_pdf       (NULL),
  m_realizer  (NULL),
  m_subCdf    (NULL),
  m_unifiedCdf(NULL),
  m_mdf       (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRVClass<V,M>::~uqBaseVectorRVClass()
{
  //if (m_mdf       ) delete m_mdf;
  //if (m_subCdf    ) delete m_subCdf;
  //if (m_unifiedCdf) delete m_unifiedCdf;
  //if (m_realizer  ) delete m_realizer;
  //if (m_pdf       ) delete m_pdf;
}

template<class V, class M>
const uqVectorSetClass<V,M>&
uqBaseVectorRVClass<V,M>::imageSet() const
{
  return m_imageSet;
}

template<class V, class M>
const uqBaseJointPdfClass<V,M>&
uqBaseVectorRVClass<V,M>::pdf() const
{
  UQ_FATAL_TEST_MACRO(m_pdf == NULL,
                      m_env.fullRank(),
                      "uqBaseVectorRVClass<V,M>::pdf()",
                      "m_pdf is NULL");

  return *m_pdf;
}

template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqBaseVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.fullRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::realizer(), prefix=")+m_prefix,
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
const uqBaseVectorCdfClass<V,M>&
uqBaseVectorRVClass<V,M>::subCdf() const
{
  UQ_FATAL_TEST_MACRO(m_subCdf == NULL,
                      m_env.fullRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::subCdf(), prefix=")+m_prefix,
                      "m_subCdf is NULL");

  return *m_subCdf;
}

template<class V, class M>
const uqBaseVectorCdfClass<V,M>&
uqBaseVectorRVClass<V,M>::unifiedCdf() const
{
  UQ_FATAL_TEST_MACRO(m_unifiedCdf == NULL,
                      m_env.fullRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::unifiedCdf(), prefix=")+m_prefix,
                      "m_unifiedCdf is NULL");

  return *m_unifiedCdf;
}

template<class V, class M>
const uqBaseVectorMdfClass<V,M>&
uqBaseVectorRVClass<V,M>::mdf() const
{
  UQ_FATAL_TEST_MACRO(m_mdf == NULL,
                      m_env.fullRank(),
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
                         const uqVectorSetClass         <V,M>& imageSet);
  uqGenericVectorRVClass(const char*                           prefix,
                         const uqVectorSetClass         <V,M>& imageSet,
                         const uqBaseJointPdfClass      <V,M>& pdf,
                         const uqBaseVectorRealizerClass<V,M>& realizer,
                         const uqBaseVectorCdfClass     <V,M>& subCdf,
                         const uqBaseVectorCdfClass     <V,M>& unifiedCdf,
                         const uqBaseVectorMdfClass     <V,M>& mdf);
  virtual ~uqGenericVectorRVClass();

  void setPdf       (uqBaseJointPdfClass      <V,M>& pdf       );
  void setRealizer  (uqBaseVectorRealizerClass<V,M>& realizer  );
  void setSubCdf    (uqBaseVectorCdfClass     <V,M>& subCdf    );
  void setUnifiedCdf(uqBaseVectorCdfClass     <V,M>& unifiedCdf);
  void setMdf       (uqBaseVectorMdfClass     <V,M>& mdf       );

  void print(std::ostream& os) const;

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSet;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_subCdf;
  using uqBaseVectorRVClass<V,M>::m_unifiedCdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                     prefix,
  const uqVectorSetClass <V,M>& imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGenericVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGenericVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                           prefix,
  const uqVectorSetClass         <V,M>& imageSet,
  const uqBaseJointPdfClass      <V,M>& pdf,
  const uqBaseVectorRealizerClass<V,M>& realizer,
  const uqBaseVectorCdfClass     <V,M>& subCdf,
  const uqBaseVectorCdfClass     <V,M>& unifiedCdf,
  const uqBaseVectorMdfClass     <V,M>& mdf)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGenericVectorRVClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = &pdf;
  m_realizer   = &realizer;
  m_subCdf     = &subCdf;
  m_unifiedCdf = &unifiedCdf;
  m_mdf        = &mdf;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGenericVectorRVClass<V,M>::constructor() [2]"
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
uqGenericVectorRVClass<V,M>::setPdf(uqBaseJointPdfClass<V,M>& pdf)
{
  m_pdf = &pdf;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setRealizer(uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setSubCdf(uqBaseVectorCdfClass<V,M>& subCdf)
{
  m_subCdf = &subCdf;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setUnifiedCdf(uqBaseVectorCdfClass<V,M>& unifiedCdf)
{
  m_unifiedCdf = &unifiedCdf;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setMdf(uqBaseVectorMdfClass<V,M>& mdf)
{
  m_mdf = &mdf;
  return;
}

template <class V, class M>
void
uqGenericVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqGenericVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGaussianVectorRVClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& imageSet,
                          const V&                     imageExpVector,
                          const V&                     imageVarVector);
  uqGaussianVectorRVClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& imageSet,
                          const V&                     imageExpVector,
                          const M&                     covMatrix);
  virtual ~uqGaussianVectorRVClass();

  void updateExpVector(const V& newExpVector);
  void updateCovMatrix(const M& newCovMatrix);
  
  void print(std::ostream& os) const;

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSet;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_subCdf;
  using uqBaseVectorRVClass<V,M>::m_unifiedCdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     imageExpVector,
  const V&                     imageVarVector)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf = new uqGaussianJointPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            imageExpVector,
                                            imageVarVector);

  M lowerCholCovMatrix(imageVarVector);
  int iRC = lowerCholCovMatrix.chol();
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
                      "uqGaussianVectorRVClass<V,M>::constructor() [1]",
                      "Cholesky decomposition of covariance matrix failed.");
  lowerCholCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
						      m_imageSet,
						      imageExpVector,
						      lowerCholCovMatrix);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     imageExpVector,
  const M&                     covMatrix)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorRVClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf = new uqGaussianJointPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            imageExpVector,
                                            covMatrix);

  M lowerCholCovMatrix(covMatrix);
  int iRC = lowerCholCovMatrix.chol();
  if (iRC) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [2]: covMatrix contents are\n";
      *m_env.subDisplayFile() << covMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
  }
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
		      "uqGaussianVectorRVClass<V,M>::constructor() [2]",
		      "Cholesky decomposition of covariance matrix failed.");
  lowerCholCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
						      m_imageSet,
						      imageExpVector,
						      lowerCholCovMatrix);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [2]"
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
uqGaussianVectorRVClass<V,M>::updateExpVector(const V& newExpVector)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass     <V,M>* >(m_pdf     ) )->updateExpVector(newExpVector);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateExpVector(newExpVector);
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateCovMatrix(const M& newCovMatrix)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass     <V,M>* >(m_pdf     ) )->updateCovMatrix(newCovMatrix);

  M newLowerCholCovMatrix(newCovMatrix);
  int iRC = newLowerCholCovMatrix.chol();
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
                      "uqGaussianVectorRVClass<V,M>::updateCovMatrix()",
                      "Cholesky decomposition of covariance matrix failed.");
  newLowerCholCovMatrix.zeroUpper(false);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLowerCholCovMatrix(newLowerCholCovMatrix);
  return;
}

template <class V, class M>
void
uqGaussianVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqGaussianVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Uniform class
//*****************************************************
template<class V, class M>
class uqUniformVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqUniformVectorRVClass(const char*                  prefix,
                         const uqVectorSetClass<V,M>& imageSet);
  virtual ~uqUniformVectorRVClass();

  void print(std::ostream& os) const;

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSet;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_subCdf;
  using uqBaseVectorRVClass<V,M>::m_unifiedCdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;
};

template<class V, class M>
uqUniformVectorRVClass<V,M>::uqUniformVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqUniformVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = new uqUniformJointPdfClass<V,M>(m_prefix.c_str(),
                                                 m_imageSet);
  m_realizer   = new uqUniformVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                       m_imageSet);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqUniformVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqUniformVectorRVClass<V,M>::~uqUniformVectorRVClass()
{
}

template <class V, class M>
void
uqUniformVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqUniformVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

template <class P_V, class P_M, class Q_V, class Q_M>
void
uqComputeCovCorrMatricesBetweenVectorRvs(
  const uqBaseVectorRVClass<P_V,P_M>& paramRv,
  const uqBaseVectorRVClass<Q_V,Q_M>& qoiRv,
        unsigned int                  localNumSamples,
        P_M&                          pqCovMatrix,
        P_M&                          pqCorrMatrix)
{
  // Check input data consistency
  const uqBaseEnvironmentClass& env = paramRv.env();

  bool useOnlyInter0Comm = (paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1) &&
                           (qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage()   == 1);

  UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorRvs()",
                      "parallel vectors not supported yet");

  unsigned int numRows = paramRv.imageSet().vectorSpace().dim();
  unsigned int numCols = qoiRv.imageSet().vectorSpace().dim();

  UQ_FATAL_TEST_MACRO((numRows != pqCovMatrix.numRows()) || (numCols != pqCovMatrix.numCols()),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorRvs()",
                      "inconsistent dimensions for covariance matrix");

  UQ_FATAL_TEST_MACRO((numRows != pqCorrMatrix.numRows()) || (numCols != pqCorrMatrix.numCols()),
                      env.fullRank(),
                      "uqComputeCorrelationBetweenVectorRvs()",
                      "inconsistent dimensions for correlation matrix");

  UQ_FATAL_TEST_MACRO((localNumSamples > paramRv.realizer().period()) || (localNumSamples > qoiRv.realizer().period()),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorRvs()",
                      "localNumSamples is too large");

  // For both P and Q vector sequences: fill them
  P_V tmpP(paramRv.imageSet().vectorSpace().zeroVector());
  Q_V tmpQ(qoiRv.imageSet().vectorSpace().zeroVector());

  uqSequenceOfVectorsClass<P_V,P_M> localWorkingPSeq(paramRv.imageSet().vectorSpace(),
                                                     localNumSamples,
                                                     "covTmpP");
  uqSequenceOfVectorsClass<Q_V,Q_M> localWorkingQSeq(qoiRv.imageSet().vectorSpace(),
                                                     localNumSamples,
                                                     "covTmpQ");
  for (unsigned int k = 0; k < localNumSamples; ++k) {
    paramRv.realizer().realization(tmpP);
    localWorkingPSeq.setPositionValues(k,tmpP);

    qoiRv.realizer().realization(tmpQ);
    localWorkingQSeq.setPositionValues(k,tmpQ);
  }

  uqComputeCovCorrMatricesBetweenVectorSequences(localWorkingPSeq,
                                                 localWorkingQSeq,
                                                 localNumSamples,
                                                 pqCovMatrix,
                                                 pqCorrMatrix);

  return;
}
#endif // __UQ_VECTOR_RV_H__
