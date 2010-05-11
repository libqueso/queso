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
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet)
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
                          const V&                     lawExpVector,
                          const V&                     lawVarVector);
  uqGaussianVectorRVClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& imageSet,
                          const V&                     lawExpVector,
                          const M&                     lawCovMatrix);
  virtual ~uqGaussianVectorRVClass();

  void updateLawExpVector(const V& newLawExpVector);
  void updateLawCovMatrix(const M& newLawCovMatrix);
  
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
  const V&                     lawExpVector,
  const V&                     lawVarVector)
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
                                            lawExpVector,
                                            lawVarVector);

  M lowerCholLawCovMatrix(lawVarVector);
  int iRC = lowerCholLawCovMatrix.chol();
  if (iRC) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [1]: lawVarVector contents are\n";
      *m_env.subDisplayFile() << lawVarVector; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (lawVarVector);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = matU.svd(matVt,vecS);
    vecS.cwSqrt();
    lowerCholLawCovMatrix = matU * leftDiagScaling(vecS,matVt);
#if 1 // For debug only
    M matOrig (lawVarVector);
    M matCheck(lowerCholLawCovMatrix * lowerCholLawCovMatrix.transpose());
    M matDiff (matOrig - matCheck);
    double frobOrig  = matOrig.normFrob();
    double frobCheck = matCheck.normFrob();
    double frobDiff  = matDiff.normFrob();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [1]"
                              << ": frobOrig = "  << frobOrig
                              << ", frobCheck = " << frobCheck
                              << ", frobDiff = "  << frobDiff
                              << ", diff/orig = " << frobDiff/frobOrig
                              << std::endl;
    }
#endif
  }
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
                      "uqGaussianVectorRVClass<V,M>::constructor() [1]",
                      "Cholesky decomposition of covariance matrix failed.");
  lowerCholLawCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
						      m_imageSet,
						      lawExpVector,
						      lowerCholLawCovMatrix);

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
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
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
                                            lawExpVector,
                                            lawCovMatrix);

  M lowerCholLawCovMatrix(lawCovMatrix);
  int iRC = lowerCholLawCovMatrix.chol();
  if (iRC) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [2]: lawCovMatrix contents are\n";
      *m_env.subDisplayFile() << lawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (lawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = matU.svd(matVt,vecS);
    vecS.cwSqrt();
    lowerCholLawCovMatrix = matU * leftDiagScaling(vecS,matVt);
#if 1 // For debug only
    M matOrig (lawCovMatrix);
    M matCheck(lowerCholLawCovMatrix * lowerCholLawCovMatrix.transpose());
    M matDiff (matOrig - matCheck);
    double frobOrig  = matOrig.normFrob();
    double frobCheck = matCheck.normFrob();
    double frobDiff  = matDiff.normFrob();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [2]"
                              << ": frobOrig = "  << frobOrig
                              << ", frobCheck = " << frobCheck
                              << ", frobDiff = "  << frobDiff
                              << ", diff/orig = " << frobDiff/frobOrig
                              << std::endl;
    }
#endif
  }
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
		      "uqGaussianVectorRVClass<V,M>::constructor() [2]",
		      "Cholesky decomposition of covariance matrix failed.");
  lowerCholLawCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
						      m_imageSet,
						      lawExpVector,
						      lowerCholLawCovMatrix);

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
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass      <V,M>* >(m_pdf     ) )->updateLawExpVector(newLawExpVector);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLawExpVector(newLawExpVector);
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass<V,M>* >(m_pdf) )->updateLawCovMatrix(newLawCovMatrix);

  M newLowerCholLawCovMatrix(newLawCovMatrix);
  int iRC = newLowerCholLawCovMatrix.chol();
  if (iRC) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::updateLawCovMatrix(): newLawCovMatrix contents are\n";
      *m_env.subDisplayFile() << newLawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (newLawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = matU.svd(matVt,vecS);
    vecS.cwSqrt();
    newLowerCholLawCovMatrix = matU * leftDiagScaling(vecS,matVt);
#if 1 // For debug only
    M matOrig (newLawCovMatrix);
    M matCheck(newLowerCholLawCovMatrix * newLowerCholLawCovMatrix.transpose());
    M matDiff (matOrig - matCheck);
    double frobOrig  = matOrig.normFrob();
    double frobCheck = matCheck.normFrob();
    double frobDiff  = matDiff.normFrob();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::updateLawCovMatrix()"
                              << ": frobOrig = "  << frobOrig
                              << ", frobCheck = " << frobCheck
                              << ", frobDiff = "  << frobDiff
                              << ", diff/orig = " << frobDiff/frobOrig
                              << std::endl;
    }
#endif
  }
  UQ_FATAL_TEST_MACRO(iRC,
                      m_env.fullRank(),
                      "uqGaussianVectorRVClass<V,M>::updateLawCovMatrix()",
                      "Cholesky decomposition of covariance matrix failed.");
  newLowerCholLawCovMatrix.zeroUpper(false);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLowerCholLawCovMatrix(newLowerCholLawCovMatrix);
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
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}

template <class V, class M>
void
uqUniformVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqUniformVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// InverseGamma class
//*****************************************************
template<class V, class M>
class uqInverseGammaVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqInverseGammaVectorRVClass(const char*                  prefix,
                              const uqVectorSetClass<V,M>& imageSet,
                              const V&                     alpha,
                              const V&                     beta);
  virtual ~uqInverseGammaVectorRVClass();

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
uqInverseGammaVectorRVClass<V,M>::uqInverseGammaVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqInverseGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = new uqInverseGammaJointPdfClass<V,M>(m_prefix.c_str(),
                                                      m_imageSet,
                                                      alpha,
                                                      beta);
  m_realizer   = new uqInverseGammaVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                            m_imageSet,
                                                            alpha,
                                                            beta);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqInverseGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqInverseGammaVectorRVClass<V,M>::~uqInverseGammaVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}

template <class V, class M>
void
uqInverseGammaVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqInverseGammaVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Concatenated class
//*****************************************************
template<class V, class M>
class uqConcatenatedVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqConcatenatedVectorRVClass(const char*                     prefix,
                              const uqBaseVectorRVClass<V,M>& rv1,
                              const uqBaseVectorRVClass<V,M>& rv2,
                              const uqVectorSetClass<V,M>&    imageSet);
  virtual ~uqConcatenatedVectorRVClass();

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

  const uqBaseVectorRVClass<V,M>& m_rv1;
  const uqBaseVectorRVClass<V,M>& m_rv2;
};

template<class V, class M>
uqConcatenatedVectorRVClass<V,M>::uqConcatenatedVectorRVClass(
  const char*                     prefix,
  const uqBaseVectorRVClass<V,M>& rv1,
  const uqBaseVectorRVClass<V,M>& rv2,
  const uqVectorSetClass<V,M>&    imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"concat").c_str(),imageSet),
  m_rv1(rv1),
  m_rv2(rv2)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = new uqConcatenatedJointPdfClass<V,M>(m_prefix.c_str(),
                                                      m_rv1.pdf(),
                                                      m_rv2.pdf(),
                                                      m_imageSet);

  m_realizer   = new uqConcatenatedVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                            m_rv1.realizer(),
                                                            m_rv2.realizer(),
                                                            m_imageSet);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqConcatenatedVectorRVClass<V,M>::~uqConcatenatedVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}

template <class V, class M>
void
uqConcatenatedVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqConcatenatedVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
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

  bool useOnlyInter0Comm = (paramRv.imageSet().vectorSpace().numOfProcsForStorage() == 1) &&
                           (qoiRv.imageSet().vectorSpace().numOfProcsForStorage()   == 1);

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
