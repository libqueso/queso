//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqVectorSpace.h>
#include <uqJointPdf.h>
#include <uqVectorRealizer.h>
#include <uqVectorCdf.h>
#include <uqVectorMdf.h>
#include <uqSequenceOfVectors.h>
#include <uqInfoTheory.h>
#include <gsl/gsl_sf_psi.h> // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

//*****************************************************
// Base class [RV-00]
//*****************************************************

/*! \file uqVectorRVClass.h
 * \brief A templated class for handling vector random variables (RV).
 * 
 * \class uqBaseVectorRVClass
 * \brief A templated base class for handling vector RV.
 *
 * This class allows two basic but quite crucial functionalities: to compute the value of the
 * PDF of a random variable (RV) at a point and to generate realizations (samples) from such PDF. */

template<class V, class M>
class uqBaseVectorRVClass {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new instance of uqBaseVectorRVClass, given a prefix and the image set of the
   * vector RV.
   */
  uqBaseVectorRVClass(const char*                  prefix,
                      const uqVectorSetClass<V,M>& imageSet);
  
  //! Virtual destructor.
  virtual ~uqBaseVectorRVClass();
  //@}
  
  //! @name Random variable-handling methods
  //@{
  //! QUESO environment; access to private attribute m_env.
  const   uqBaseEnvironmentClass&         env       () const;
  
  //! Image set of the vector RV; access to private attribute m_imageSet.
  const   uqVectorSetClass         <V,M>& imageSet  () const;
  
  //! Posterior Density Function of the vector RV; access to private attribute m_pdf.
  const   uqBaseJointPdfClass      <V,M>& pdf       () const;
  
  //! Finds a realization (sample) of the PDF of this vector RV; access to private attribute m_realizer.
  const   uqBaseVectorRealizerClass<V,M>& realizer  () const;
  
  //! Finds the Cumulative Distribution Function of this vector RV, considering only the sub-sequence of data; access to private attribute m_subCdf.
  const   uqBaseVectorCdfClass     <V,M>& subCdf    () const;
  
  //! Finds the Cumulative Distribution Function of this vector RV, considering the unified sequence of data; access to private attribute m_unifiedCdf.
  const   uqBaseVectorCdfClass     <V,M>& unifiedCdf() const;
  
  //! Finds the Mass Density Function of this vector RV; access to private attribute m_mdf.
  const   uqBaseVectorMdfClass     <V,M>& mdf       () const;
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  virtual void                            print     (std::ostream& os) const = 0;
  //@}
  
#ifdef QUESO_HAS_ANN
  virtual double                          estimateENT_ANN() const;
  /*
  virtual double                          estimateENT_ANN( unsigned int k, double eps ) const;
  virtual double                          estimateENTSubset_ANN( const unsigned int dimSel[] ) const;
  virtual double                          estimateENTSubset_ANN( const unsigned int dimSel[], unsigned int k, double eps ) const;
  */
#endif // QUESO_HAS_ANN
  //@}
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
// Default constructor -----------------------------
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
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqBaseVectorRVClass<V,M>::~uqBaseVectorRVClass()
{
  //if (m_mdf       ) delete m_mdf;
  //if (m_subCdf    ) delete m_subCdf;
  //if (m_unifiedCdf) delete m_unifiedCdf;
  //if (m_realizer  ) delete m_realizer;
  //if (m_pdf       ) delete m_pdf;
}
// RV handling-methods ------------------------------

template <class V, class M>
const uqBaseEnvironmentClass&
uqBaseVectorRVClass<V,M>::env() const
{
  return m_env;
}
//---------------------------------------------------
template<class V, class M>
const uqVectorSetClass<V,M>&
uqBaseVectorRVClass<V,M>::imageSet() const
{
  return m_imageSet;
}
//---------------------------------------------------
template<class V, class M>
const uqBaseJointPdfClass<V,M>&
uqBaseVectorRVClass<V,M>::pdf() const
{
  UQ_FATAL_TEST_MACRO(m_pdf == NULL,
                      m_env.worldRank(),
                      "uqBaseVectorRVClass<V,M>::pdf()",
                      "m_pdf is NULL");

  return *m_pdf;
}
//---------------------------------------------------
template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqBaseVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.worldRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::realizer(), prefix=")+m_prefix,
                      "m_realizer is NULL");

  return *m_realizer;
}
//---------------------------------------------------
template<class V, class M>
const uqBaseVectorCdfClass<V,M>&
uqBaseVectorRVClass<V,M>::subCdf() const
{
  UQ_FATAL_TEST_MACRO(m_subCdf == NULL,
                      m_env.worldRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::subCdf(), prefix=")+m_prefix,
                      "m_subCdf is NULL");

  return *m_subCdf;
}
//---------------------------------------------------
template<class V, class M>
const uqBaseVectorCdfClass<V,M>&
uqBaseVectorRVClass<V,M>::unifiedCdf() const
{
  UQ_FATAL_TEST_MACRO(m_unifiedCdf == NULL,
                      m_env.worldRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::unifiedCdf(), prefix=")+m_prefix,
                      "m_unifiedCdf is NULL");

  return *m_unifiedCdf;
}
//---------------------------------------------------
template<class V, class M>
const uqBaseVectorMdfClass<V,M>&
uqBaseVectorRVClass<V,M>::mdf() const
{
  UQ_FATAL_TEST_MACRO(m_mdf == NULL,
                      m_env.worldRank(),
                      (std::string)("uqBaseVectorRVClass<V,M>::mdf(), prefix=")+m_prefix,
                      "m_mdf is NULL");

  return *m_mdf;
}

//---------------------------------------------------
// Operator declared outside class definition -------
//---------------------------------------------------
template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqBaseVectorRVClass<V,M>& obj)
{
  obj.print(os);

  return os;
}
//---------------------------------------------------
#ifdef QUESO_HAS_ANN
template <class V, class M>
double 
uqBaseVectorRVClass<V,M>::estimateENT_ANN() const
{
  ANNpointArray data;
  double* dists;
  double ENT_est;

  // FIXME: these default values should be stored in the
  // QUESO input file ( create a uqInfoTheoryOptionsClass )
  unsigned int k = UQ_INFTH_ANN_KNN;
  double eps = UQ_INFTH_ANN_EPS;

  // here it is assumed that the entropy for the
  // entire joint RV will be computed 
  unsigned int dim = this->imageSet().vectorSpace().dimGlobal();

  // FIXME: get the number already stored, otherwise
  // use the default value
  unsigned int N = this->realizer().subPeriod();
  if( N == 0 ) {
    N = UQ_INFTH_ANN_NO_SMP;
  }

  // allocate memory
  data = annAllocPts(N,dim);
  dists = new double[N];
  
  // copy samples in the ANN data structure
  V smpRV( this->imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < N; i++ ) {
    // get a sample from the distribution
    this->realizer().realization( smpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dim; j++ ) {
      data[ i ][ j ] = smpRV[ j ];
    }
  }
  
  // get distance to knn for each point
  // (k+1) because the 1st nn is itself
  distANN_XY( data, data, dists, dim, dim, N, N, k+1, eps );

  // compute the entropy estimate using the L-infinity (Max) norm
  // so no need for the adjustment of the mass of the hyperball
  // this has to be enforced before compiling the ANN lib
  double sum_log_dist = 0.0;
  for( unsigned int i = 0; i < N; i++ ) {
    if( dists[ i ] > 0 ) {
      sum_log_dist += log( 2.0*dists[ i ] );
    }
  }
  ENT_est = - gsl_sf_psi_int( k ) + gsl_sf_psi_int( N ) + (double)dim / (double)N * sum_log_dist; // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

  // deallocate memory
  delete [] dists;
  annDeallocPts( data );  

  return ENT_est;
}
#endif // QUESO_HAS_ANN

//*****************************************************
// Generic class [RV-01]
//*****************************************************
 /*! \class uqGenericVectorRVClass
 * \brief A templated class for handling generic vector RVs.
 *
 * This class allows the user to compute the value of the PDF of a generic random variable (RV)
 * and to generate realizations (samples) from such PDF.  This is the class used by QUESO to 
 * store the solution of an statistical inverse problem. */
 
template<class V, class M>
class uqGenericVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new instance, given a prefix and the image set of the vector RV.*/
  uqGenericVectorRVClass(const char*                           prefix,
                         const uqVectorSetClass         <V,M>& imageSet);
  
  //! Constructor
  /*! Constructs a new instance, given all the attributes that characterize the vector RV: prefix, image set, pdf, etc.*/
  uqGenericVectorRVClass(const char*                           prefix,
                         const uqVectorSetClass         <V,M>& imageSet,
                         const uqBaseJointPdfClass      <V,M>& pdf,
                         const uqBaseVectorRealizerClass<V,M>& realizer,
                         const uqBaseVectorCdfClass     <V,M>& subCdf,
                         const uqBaseVectorCdfClass     <V,M>& unifiedCdf,
                         const uqBaseVectorMdfClass     <V,M>& mdf);
  //! Virtual destructor
  virtual ~uqGenericVectorRVClass();
  //@}
  
    //! @name Random variable-handling methods
  //@{
  //! Sets the PDF of \c this vector RV  to \c pdf.  
  void setPdf       (uqBaseJointPdfClass      <V,M>& pdf       );
  
  //! Sets the realizer of \c this vector RV  to \c realizer.  
  void setRealizer  (uqBaseVectorRealizerClass<V,M>& realizer  );
  
  //! Sets the CDF of the sub-sequence of \c this vector RV  to \c subCdf.
  void setSubCdf    (uqBaseVectorCdfClass     <V,M>& subCdf    );
  
  //! Sets the CDF of the unified sequence of \c this vector RV  to \c unifiedCdf.
  void setUnifiedCdf(uqBaseVectorCdfClass     <V,M>& unifiedCdf);
  
  //! Sets the MDF of  \c this vector RV  to \c Mdf.
  void setMdf       (uqBaseVectorMdfClass     <V,M>& mdf       );
  //@}  
   
    //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

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
// Default constructor -----------------------------
template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                     prefix,
  const uqVectorSetClass <V,M>& imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
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
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericVectorRVClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = &pdf;
  m_realizer   = &realizer;
  m_subCdf     = &subCdf;
  m_unifiedCdf = &unifiedCdf;
  m_mdf        = &mdf;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericVectorRVClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqGenericVectorRVClass<V,M>::~uqGenericVectorRVClass()
{
}
// Random variable-handling methods ----------------
template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setPdf(uqBaseJointPdfClass<V,M>& pdf)
{
  m_pdf = &pdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setRealizer(uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setSubCdf(uqBaseVectorCdfClass<V,M>& subCdf)
{
  m_subCdf = &subCdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setUnifiedCdf(uqBaseVectorCdfClass<V,M>& unifiedCdf)
{
  m_unifiedCdf = &unifiedCdf;
  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setMdf(uqBaseVectorMdfClass<V,M>& mdf)
{
  m_mdf = &mdf;
  return;
}
//--------------------------------------------------
template <class V, class M>
void
uqGenericVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqGenericVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Gaussian class [RV-03]
//*****************************************************

/*!
 * \class uqGaussianVectorRVClass
 * \brief A class representing a Gaussian vector RV.
 * 
 * This class allows the user to compute the value of a Gaussian PDF and to generate realizations
 * (samples) from it.
 * 
 * In probability theory, the normal (or Gaussian) distribution is a continuous probability 
 * distribution, defined by the formula:
 * \f[    f(x| \mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} e^{ -\frac{(x-\mu)^2}{2\sigma^2} }. \f]
 * 
 * The parameter \f$ \mu \f$  in this formula is the mean or expectation of the distribution (and also 
 * its median and mode). The parameter \f$ \sigma \f$  is its standard deviation; its variance is therefore
 * \f$ \sigma^2 \f$ . */

template<class V, class M>
class uqGaussianVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor  
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and diagonal covariance matrix
   * \c lawVarVector whose variates live in \c imageSet.*/
  uqGaussianVectorRVClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& imageSet,
                          const V&                     lawExpVector,
                          const V&                     lawVarVector);
  
  //! Constructor  
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and covariance matrix
   * \c lawCovMatrix whose variates live in \c imageSet.*/
  uqGaussianVectorRVClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& imageSet,
                          const V&                     lawExpVector,
                          const M&                     lawCovMatrix);
  
  //! Virtual destructor
  virtual ~uqGaussianVectorRVClass();
  //@}

  //! @name Statistical methods
  //@{
  //! Updates the vector that contains the mean values.
  void updateLawExpVector(const V& newLawExpVector);
  
  //! Updates the covariance matrix.
  /*! This method tries to use Cholesky decomposition; and if it fails, the method then 
   *  calls a SVD decomposition.*/
  void updateLawCovMatrix(const M& newLawCovMatrix);
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
 //@}
  
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
// Constructor---------------------------------------
template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((lawVarVector.getMinValue() <= 0.0),
                      m_env.worldRank(),
                      "uqGaussianVectorRVClass<V,M>::constructor() [1]",
                      "Covariance matrix is not symmetric positive definite.");

  m_pdf = new uqGaussianJointPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            lawExpVector,
                                            lawVarVector);

  V cholDiag(lawVarVector);
  cholDiag.cwSqrt();
  M lowerCholLawCovMatrix(cholDiag);
  lowerCholLawCovMatrix.zeroUpper(false);

  m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                      m_imageSet,
                                                      lawExpVector,
                                                      lowerCholLawCovMatrix);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor---------------------------------------
template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
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
  lowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In uqGaussianVectorRVClass<V,M>::constructor() [2]: chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::constructor() [2]: chol failed; will use svd; lawCovMatrix contents are\n";
      *m_env.subDisplayFile() << lawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (lawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = lawCovMatrix.svd(matU,vecS,matVt);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.worldRank(),
                        "uqGaussianVectorRVClass<V,M>::constructor() [2]",
		        "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        matU,
                                                        vecS, // already square rooted
                                                        matVt);
  }
  else {
    m_realizer = new uqGaussianVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        lowerCholLawCovMatrix);
  }

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqGaussianVectorRVClass<V,M>::~uqGaussianVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// Statistical methods-------------------------------
template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass      <V,M>* >(m_pdf     ) )->updateLawExpVector(newLawExpVector);
  ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLawExpVector(newLawExpVector);
  return;
}
//---------------------------------------------------
template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< uqGaussianJointPdfClass<V,M>* >(m_pdf) )->updateLawCovMatrix(newLawCovMatrix);

  M newLowerCholLawCovMatrix(newLawCovMatrix);
  int iRC = newLowerCholLawCovMatrix.chol();
  newLowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In uqGaussianVectorRVClass<V,M>::updateLawCovMatrix(): chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGaussianVectorRVClass<V,M>::updateLawCovMatrix(): chol failed; will use svd; newLawCovMatrix contents are\n";
      *m_env.subDisplayFile() << newLawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (newLawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = newLawCovMatrix.svd(matU,vecS,matVt);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.worldRank(),
                        "uqGaussianVectorRVClass<V,M>::updateLawCovMatrix()",
                        "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLowerCholLawCovMatrix(matU,
                                                                                                     vecS, // already square rooted
                                                                                                     matVt);
  }
  else {
    ( dynamic_cast< uqGaussianVectorRealizerClass<V,M>* >(m_realizer) )->updateLowerCholLawCovMatrix(newLowerCholLawCovMatrix);
  }
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqGaussianVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqGaussianVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}


//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
template<class V, class M>
void
uqComputeConditionalGaussianVectorRV(
  const V& muVec1,
  const V& muVec2,
  const M& sigmaMat11,
  const M& sigmaMat12,
  const M& sigmaMat21,
  const M& sigmaMat22,
  const V& sampleVec2,
        V& muVec1_cond_on_2,
        M& sigmaMat11_cond_on_2)
{
  const uqBaseEnvironmentClass& env = muVec1.env();
  unsigned int dim1 = muVec1.sizeLocal();
  unsigned int dim2 = muVec2.sizeLocal();

  UQ_FATAL_TEST_MACRO((sigmaMat11.numRowsLocal() != dim1) || (sigmaMat11.numCols() != dim1),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sigmaMat11");

  UQ_FATAL_TEST_MACRO((sigmaMat12.numRowsLocal() != dim1) || (sigmaMat12.numCols() != dim2),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sigmaMat12");

  UQ_FATAL_TEST_MACRO((sigmaMat21.numRowsLocal() != dim2) || (sigmaMat21.numCols() != dim1),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sigmaMat21");

  UQ_FATAL_TEST_MACRO((sigmaMat22.numRowsLocal() != dim2) || (sigmaMat22.numCols() != dim2),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sigmaMat22");

  // Check transpose operation
  M mat_tt(sigmaMat12);
  mat_tt.cwSet(0.);
  mat_tt.fillWithTranspose(0,0,sigmaMat21,true,true);
  double auxNorm = (mat_tt - sigmaMat12).normFrob();
  if (auxNorm >= 1.e-12) {
    if (env.subDisplayFile()) {
      *env.subDisplayFile() << "In uqComputeConditionalGaussianVectorRV()"
                            << ": WARNING, ||sigmaMat21^T - sigmaMat12||_2 = " << auxNorm
                            << std::endl;
    }
  }
  UQ_FATAL_TEST_MACRO(auxNorm >= 1.e-12,
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "sigmaMat12 and sigmaMat21 are not transpose of each other");

  UQ_FATAL_TEST_MACRO((sampleVec2.sizeLocal() != dim2),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sampleVec2");

  UQ_FATAL_TEST_MACRO((muVec1_cond_on_2.sizeLocal() != dim1),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid muVec1_cond_on_2");

  UQ_FATAL_TEST_MACRO((sigmaMat11_cond_on_2.numRowsLocal() != dim1) || (sigmaMat11_cond_on_2.numCols() != dim1),
                      env.worldRank(),
                      "uqComputeConditionalGaussianVectorRV()",
                      "invalid sigmaMat11_cond_on_2");

  muVec1_cond_on_2     = muVec1     + sigmaMat12 * sigmaMat22.invertMultiply(sampleVec2 - muVec2);
  sigmaMat11_cond_on_2 = sigmaMat11 - sigmaMat12 * sigmaMat22.invertMultiply(sigmaMat21);

  return;
}

//*****************************************************
// Uniform class [RV-04]
//*****************************************************
/*!
 * \class uqUniformVectorRVClass
 * \brief A class representing a uniform vector RV.
 * 
 * This class allows the user to compute the value of a uniform PDF and to generate realizations
 * (samples) from it. It is used, for instance, to create a uniform prior PDF. */

template<class V, class M>
class uqUniformVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a uniform vector RV, given a prefix and the image set of the vector RV.*/
  uqUniformVectorRVClass(const char*                  prefix,
                         const uqVectorSetClass<V,M>& imageSet);
  //! Virtual destructor
  virtual ~uqUniformVectorRVClass();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

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
// Default constructor-------------------------------
template<class V, class M>
uqUniformVectorRVClass<V,M>::uqUniformVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
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

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqUniformVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqUniformVectorRVClass<V,M>::~uqUniformVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods --------------------------------------
template <class V, class M>
void
uqUniformVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqUniformVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Beta class [RV-05]
//*****************************************************
/*!
 * \class uqBetaVectorRVClass
 * \brief A class representing a vector RV constructed via Beta distribution.
 * 
 * This class allows the user to compute the value of a Beta PDF and to generate realizations
 * (samples) from it.\n
 * 
 * The beta probability density function for a given value x and given pair of parameters 
 * \b a and \b b is: 
 *  \f[ y=f(x|a,b)= \frac{1}{B(a,b)} x^{a-1}(1-x)^{b-1}, \f]
 * where <b>B(·)</b> is the Beta function:
 * \f[  B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}=\frac{(a-1)!(b-1)!}{(a+b-1)!}.\f] 
 * The parameters \b a and \b b must all be positive, and the values \c x  must lie on the 
 * interval [0, 1].
 */


template<class V, class M>
class uqBetaVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct a Beta vector RV with parameters \c a>0  and \c b>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to [0,1], which
   * is a requirement imposed by the Beta distribution. If this condition is not satisfied, an error 
   * message will be displayed and the program will exit. */
  uqBetaVectorRVClass(const char*                  prefix,
                      const uqVectorSetClass<V,M>& imageSet,
                      const V&                     alpha,
                      const V&                     beta);
  
  //! Virtual destructor
  virtual ~uqBetaVectorRVClass();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

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
// Constructor---------------------------------------
template<class V, class M>
uqBetaVectorRVClass<V,M>::uqBetaVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBetaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-22 : --------------------------

  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&imageSet);

  double smallerOfMaxValues = imageBox->maxValues().getMinValue();	
  double biggerOfMaxValues = imageBox->maxValues().getMaxValue();	
  double smallerOfMinValues = imageBox->minValues().getMinValue();
  double biggerOfMinValues = imageBox->minValues().getMaxValue();	
	
 // Beta dist is defined only in [0,1]		
 if( (smallerOfMinValues < 0) || ( biggerOfMaxValues > 1 ) ) 
 {		
   std::cerr << "In uqBetaVectorRVClass<V,M>::constructor()\n" 
			 << "Beta distribution is defined only in [0, 1].\n"
			 << "The data provided is: \n"
			 << *imageBox 
   			 << "Sampling will not cover all interval.\n"   
   			 << std::endl;

 // if at least one of the min values > 1 then exit
    UQ_FATAL_TEST_MACRO(biggerOfMinValues > 1,
                      m_env.worldRank(),
                      "In uqBetaVectorRVClass<V,M>::constructor()",
                      "invalid input: Beta distribution is only defined in [0, 1], and max(m_minValues)>1");  
 // if at least one of the max values < 0 then exit                      
    UQ_FATAL_TEST_MACRO(smallerOfMaxValues < 0, //biggerOfMaxValues
                      m_env.worldRank(),
                      "In uqBetaVectorRVClass<V,M>::constructor()",
                      "invalid input: Beta distribution is only defined in [0, 1], and min(m_maxValues)<0");                
 }	
  // end kemelli 2013-April-22 --------------------------
  
  m_pdf        = new uqBetaJointPdfClass<V,M>(m_prefix.c_str(),
                                              m_imageSet,
                                              alpha,
                                              beta);
  m_realizer   = new uqBetaVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                    m_imageSet,
                                                    alpha,
                                                    beta);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBetaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqBetaVectorRVClass<V,M>::~uqBetaVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqBetaVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqBetaVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Gamma class [RV-06]
//*****************************************************
/*!
 * \class uqGammaVectorRVClass
 * \brief A class representing a vector RV constructed via Gamma distribution.
 * 
 * This class allows the user to compute the value of a Gamma PDF and to generate realizations
 * (samples) from it.\n
 * 
 * The gamma probability density function for a given value x and given pair of parameters 
 * \b a and \b b is: 
 *  \f[ y=f(x|a,b)= \frac{1}{b^{a}\Gamma(a)} x^{a-1} e^{\frac{x}{b}}, \f]
 * where \f$ \Gamma(.) \f$ is the Gamma function:
 * \f[  B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}=\frac{(a-1)!(b-1)!}{(a+b-1)!}.\f] 
 * The parameters \b a and \b b must all be positive, and the values \c x  must lie on the 
 * interval \f$ (0, \infty)\f$. */
 	
template<class V, class M>
class uqGammaVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  
      //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct a Gamma vector RV with parameters \c a>0  and \c b>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to 
   * \f$ (0, \infty)\f$, which is a requirement imposed by the Gamma distribution. If this condition 
   * is not satisfied, an error  message will be displayed and the program will exit. */
  uqGammaVectorRVClass(const char*                  prefix,
                       const uqVectorSetClass<V,M>& imageSet,
                       const V&                     a,
                       const V&                     b);
  //! Virtual destructor
  virtual ~uqGammaVectorRVClass();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}
  
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
// Constructor---------------------------------------
// TODO: Check: the constructor receives imageSet, but uses m_imageSet to assign m_pdf and m_realizer. (Kemelli 2013/4/22)
template<class V, class M>
uqGammaVectorRVClass<V,M>::uqGammaVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     a,
  const V&                     b)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
 
// begin kemelli 2013-April-22 -------------------------- 
// better to check for the parameter values in the constructor, 
// rather than in uqGammaVectorRealizerClass<V,M>::realization

  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&imageSet);
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();	
  double smallerOfMinValues = imageBox->minValues().getMinValue();
	
 // Gamma dist is defined only in (0,inf)		
 if( smallerOfMinValues < 0 ) 
 {		
   std::cerr << "In uqGammaVectorRVClass<V,M>::constructor()\n" 
			 << "Gamma distribution is only defined in (0, infinity).\n"
			 << "The data provided is: \n"
			 << *imageBox 
   			 << "Sampling will not cover all interval.\n"   
   			 << std::endl;

    UQ_FATAL_TEST_MACRO(smallerOfMaxValues < 0,
                      m_env.worldRank(),
                      "uqGammaVectorRealizerClass<V,M>::constructor()",
                      "invalid input: Gamma distribution is only defined in (0, infinity), and min(m_maxValues)<0");      
 }	
// end kemelli 2013-April-22 --------------------------

  m_pdf        = new uqGammaJointPdfClass<V,M>(m_prefix.c_str(),
                                               m_imageSet,
                                               a,
                                               b);
  m_realizer   = new uqGammaVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                     m_imageSet,
                                                     a,
                                                     b);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqGammaVectorRVClass<V,M>::~uqGammaVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqGammaVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqGammaVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// InverseGamma class [RV-07]
//*****************************************************
/*!
 * \class uqInverseGammaVectorRVClass
 * \brief A class representing a vector RV constructed via Inverse Gamma distribution.
 * 
 * This class allows the user to compute the value of a Inverse Gamma PDF and to generate realizations
 * (samples) from it.\n
 * 
 * The Inverse Gamma probability density function  is defined over the support x > 0 given a shape parameters 
 * \b a and a scale parameter \b b is: 
 *  \f[ y=f(x| a,b)= \frac{b^a}{\Gamma(a)} x^{-a - 1}\exp\left(-\frac{b}{x}\right)\f]
 * where \f$ \Gamma(.) \f$ is the Gamma function:
 * \f[  B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}=\frac{(a-1)!(b-1)!}{(a+b-1)!}.\f] 
 * The parameters \b a and \b b must all be positive, and the values \c x  must lie on the 
 * interval \f$ (0, \infty)\f$. */
template<class V, class M>
class uqInverseGammaVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct an Inverse Gamma vector RV with parameters \c alpha>0  and \c beta>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to 
   * \f$ (0, \infty)\f$, which is a requirement imposed by the Inverse Gamma distribution. If this 
   * condition is not satisfied, an error  message will be displayed and the program will exit. */
  uqInverseGammaVectorRVClass(const char*                  prefix,
                              const uqVectorSetClass<V,M>& imageSet,
                              const V&                     alpha,
                              const V&                     beta);
  //! Virtual destructor
  virtual ~uqInverseGammaVectorRVClass();
  //@}
    //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}
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
// Constructor---------------------------------------
template<class V, class M>
uqInverseGammaVectorRVClass<V,M>::uqInverseGammaVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqInverseGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-22 -------------------------- 
// InverseGamma dist is defined only in (0,inf)
  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&imageSet);
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();	
  double smallerOfMinValues = imageBox->minValues().getMinValue();
 		
 if( smallerOfMinValues < 0 ) 
 {		
   std::cerr << "In uqInverseGammaVectorRVClass<V,M>::constructor()\n" 
			 << "Inverse Gamma distribution is only defined in (0, infinity).\n"
			 << "The data provided is: \n"
			 << *imageBox 
   			 << "Sampling will not cover all interval.\n"   
   			 << std::endl;

    UQ_FATAL_TEST_MACRO(smallerOfMaxValues < 0,
                      m_env.worldRank(),
                      "uqInverseGammaVectorRealizerClass<V,M>::constructor()",
                      "invalid input: Inverse Gamma distribution is only defined in (0, infinity), and min(m_maxValues)<0");      
 }	
// end kemelli 2013-April-22 --------------------------
  
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

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqInverseGammaVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqInverseGammaVectorRVClass<V,M>::~uqInverseGammaVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqInverseGammaVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqInverseGammaVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// Wigner class [RV-09]
//*****************************************************
/*!
 * \class uqWignerVectorRVClass
 * \brief A class representing a vector RV constructed via Wigner distribution.
 * 
 * This class allows the user to compute the value of a Wigner PDF and to generate realizations
 * (samples) from it.\n
 * 
 * \todo: uqWignerVectorRealizerClass.realization() is not yet available, thus this class does 
 * nothing. */
 
template<class V, class M>
class uqWignerVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  
    
  //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  uqWignerVectorRVClass(const char*                  prefix,
                        const uqVectorSetClass<V,M>& imageSet,
                        const V&                     centerPos,
                        double                       radius);
  //! Virtual destructor
  virtual ~uqWignerVectorRVClass();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}
  
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
// Default constructor -----------------------------
template<class V, class M>
uqWignerVectorRVClass<V,M>::uqWignerVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     centerPos,
  double                       radius)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqWignerVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(radius <= 0.,
                      m_env.worldRank(),
                      "uqWignerVectorRVClass<V,M>::constructor()",
                      "invalid radius");

  m_pdf        = new uqWignerJointPdfClass<V,M>(m_prefix.c_str(),
                                                m_imageSet,
                                                centerPos,
                                                radius);
  m_realizer   = new uqWignerVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                      m_imageSet,
                                                      centerPos,
                                                      radius);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqWignerVectorRVClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
uqWignerVectorRVClass<V,M>::~uqWignerVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
//--------------------------------------------------
template <class V, class M>
void
uqWignerVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqWignerVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

//*****************************************************
// LogNormal class [RV-10]
//*****************************************************
/*!
 * \class uqLogNormalVectorRVClass
 * \brief A class representing a LogNormal vector RV.
 * 
 * This class allows the user to compute the value of a LogNormal PDF and to generate realizations
 * (samples) from it.\n
 * 
 * The probability density function of a log-normal distribution is: 
 *   \f[ f(x|\mu,\sigma) = \frac{1}{x \sigma \sqrt{2 \pi}}\, e^{-\frac{(\ln x - \mu)^2}{2\sigma^2}}, x>0 \f]
 * where  the parameters denoted \f$ \mu \f$ and \f$ \sigma \f$ are, respectively, the mean and standard 
 * deviation of the  variable's natural logarithm; and \c x>0.  */

template<class V, class M>
class uqLogNormalVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor  
  /*! Construct a LogNormal vector RV with mean \c lawExpVector and diagonal covariance matrix
   * \c lawVarVector whose variates live in \c imageSet.*/
  uqLogNormalVectorRVClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& imageSet,
                           const V&                     lawExpVector,
                           const V&                     lawVarVector);
  
  //! Virtual destructor
  virtual ~uqLogNormalVectorRVClass();
  //@}
 
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

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
// Constructor---------------------------------------
template<class V, class M>
uqLogNormalVectorRVClass<V,M>::uqLogNormalVectorRVClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& imageSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-23 -------------------------- 
// LogNormal dist is defined only in (0,inf)

  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&imageSet);
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();	
  double smallerOfMinValues = imageBox->minValues().getMinValue();
 		
 if( smallerOfMinValues < 0 ) 
 {		
   std::cerr << "In uqLogNormalVectorRVClass<V,M>::constructor()\n" 
			 << "LogNormal distribution is only defined in (0, infinity).\n"
			 << "The data provided is: \n"
			 << *imageBox 
   			 << "Sampling will not cover all interval.\n"   
   			 << std::endl;


    UQ_FATAL_TEST_MACRO(smallerOfMaxValues < 0,
                      m_env.worldRank(),
                      "uqLogNormalVectorRealizerClass<V,M>::constructor()",
                      "invalid input: LogNormal distribution is only defined in (0, infinity), and min(m_maxValues)<0");      
              
 }	
// end kemelli 2013-April-23 --------------------------

  m_pdf = new uqLogNormalJointPdfClass<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            lawExpVector,
                                            lawVarVector);

  M lowerCholLawCovMatrix(lawVarVector);
  int iRC = lowerCholLawCovMatrix.chol();
  lowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In uqLogNormalVectorRVClass<V,M>::constructor() [1]: chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqLogNormalVectorRVClass<V,M>::constructor() [1]: chol failed; will use svd; lawVarVector contents are\n";
      *m_env.subDisplayFile() << lawVarVector; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matLaw(lawVarVector);
    M matU  (lawVarVector);
    M matVt (m_imageSet.vectorSpace().zeroVector());
    V vecS  (m_imageSet.vectorSpace().zeroVector());
    iRC = matLaw.svd(matU,vecS,matVt);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.worldRank(),
                        "uqLogNormalVectorRVClass<V,M>::constructor() [1]",
                        "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    m_realizer = new uqLogNormalVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                         m_imageSet,
                                                         lawExpVector,
                                                         matU,
                                                         vecS, // already square rooted
                                                         matVt);
  }
  else {
    m_realizer = new uqLogNormalVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        lowerCholLawCovMatrix);
  }

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalVectorRVClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqLogNormalVectorRVClass<V,M>::~uqLogNormalVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqLogNormalVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqLogNormalVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}


//*****************************************************
// Concatenated class [RV-11]
//*****************************************************
/*!
 * \class uqConcatenatedVectorRVClass
 * \brief A class representing concatenated vector RVs.
 * 
 * This class allows the user to concatenate two vector RV of different types and to generate realizations
 * (samples) from this concatenated vector RV. It is used, for instance, to concatenate priors from two or 
 * more RVs, where one of them has a uniform distribution whereas the other one(s) has a Gaussian distribution. */

template<class V, class M>
class uqConcatenatedVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Concatenates two RVs: \c rv1 and \c rv2 into one vector RV, given a prefix and the image set of the vector RV.*/
  uqConcatenatedVectorRVClass(const char*                     prefix,
                              const uqBaseVectorRVClass<V,M>& rv1,
                              const uqBaseVectorRVClass<V,M>& rv2,
                              const uqVectorSetClass<V,M>&    imageSet);
  
  //! Constructor
  /*! Concatenates a sequence of RVs, given by: <c> std::vector<const uqBaseVectorRVClass<V,M>* >& rvs </c>
   * into one single vector RV, given a prefix and the image set of the resulting vector RV.*/
  uqConcatenatedVectorRVClass(const char*                                          prefix,
                              const std::vector<const uqBaseVectorRVClass<V,M>* >& rvs,
                              const uqVectorSetClass<V,M>&                         imageSet);
  
  //! Virtual destructor
  virtual ~uqConcatenatedVectorRVClass();
  //@}
  
    //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}
  
  
private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSet;
  using uqBaseVectorRVClass<V,M>::m_pdf;
  using uqBaseVectorRVClass<V,M>::m_realizer;
  using uqBaseVectorRVClass<V,M>::m_subCdf;
  using uqBaseVectorRVClass<V,M>::m_unifiedCdf;
  using uqBaseVectorRVClass<V,M>::m_mdf;

  std::vector<const uqBaseVectorRVClass      <V,M>* > m_rvs;
  std::vector<const uqBaseJointPdfClass      <V,M>* > m_pdfs;
  std::vector<const uqBaseVectorRealizerClass<V,M>* > m_realizers;
};
// Constructor -------------------------------------
template<class V, class M>
uqConcatenatedVectorRVClass<V,M>::uqConcatenatedVectorRVClass(
  const char*                     prefix,
  const uqBaseVectorRVClass<V,M>& rv1,
  const uqBaseVectorRVClass<V,M>& rv2,
  const uqVectorSetClass<V,M>&    imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"concat").c_str(),imageSet),
  m_rvs                   (2,(const uqBaseVectorRVClass      <V,M>*) NULL),
  m_pdfs                  (2,(const uqBaseJointPdfClass      <V,M>*) NULL),
  m_realizers             (2,(const uqBaseVectorRealizerClass<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedVectorRVClass<V,M>::constructor(1)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_rvs[0]       = &rv1;
  m_rvs[1]       = &rv2;
  m_pdfs[0]      = &(m_rvs[0]->pdf());
  m_pdfs[1]      = &(m_rvs[1]->pdf());
  m_realizers[0] = &(m_rvs[0]->realizer());
  m_realizers[1] = &(m_rvs[1]->realizer());

  m_pdf          = new uqConcatenatedJointPdfClass<V,M>(m_prefix.c_str(),
                                                        *(m_pdfs[0]),
                                                        *(m_pdfs[1]),
                                                        m_imageSet);

  m_realizer     = new uqConcatenatedVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                              *(m_realizers[0]),
                                                              *(m_realizers[1]),
                                                              m_imageSet);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedVectorRVClass<V,M>::constructor(1)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
uqConcatenatedVectorRVClass<V,M>::uqConcatenatedVectorRVClass(
  const char*                                          prefix,
  const std::vector<const uqBaseVectorRVClass<V,M>* >& rvs,
  const uqVectorSetClass<V,M>&                         imageSet)
  :
  uqBaseVectorRVClass<V,M>(((std::string)(prefix)+"concat").c_str(),imageSet),
  m_rvs                   (rvs.size(),(const uqBaseVectorRVClass      <V,M>*) NULL),
  m_pdfs                  (rvs.size(),(const uqBaseJointPdfClass      <V,M>*) NULL),
  m_realizers             (rvs.size(),(const uqBaseVectorRealizerClass<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedVectorRVClass<V,M>::constructor(2)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    m_rvs [i]      = rvs[i];
    m_pdfs[i]      = &(m_rvs[i]->pdf());
    m_realizers[i] = &(m_rvs[i]->realizer());
  }

  m_pdf        = new uqConcatenatedJointPdfClass<V,M>(m_prefix.c_str(),
                                                      m_pdfs,
                                                      m_imageSet);

  unsigned int minPeriod = m_realizers[0]->subPeriod();
  for (unsigned int i = 0; i < m_realizers.size(); ++i) {
    if (minPeriod > m_realizers[i]->subPeriod()) {
      minPeriod = m_realizers[i]->subPeriod();
    }
  }

  m_realizer   = new uqConcatenatedVectorRealizerClass<V,M>(m_prefix.c_str(),
                                                            m_realizers,
                                                            minPeriod,
                                                            m_imageSet);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedVectorRVClass<V,M>::constructor(2)"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqConcatenatedVectorRVClass<V,M>::~uqConcatenatedVectorRVClass()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
//--------------------------------------------------
template <class V, class M>
void
uqConcatenatedVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "uqConcatenatedVectorRVClass<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}
//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
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
                      env.worldRank(),
                      "uqComputeCovCorrMatricesBetweenVectorRvs()",
                      "parallel vectors not supported yet");

  unsigned int numRows = paramRv.imageSet().vectorSpace().dim();
  unsigned int numCols = qoiRv.imageSet().vectorSpace().dim();

  UQ_FATAL_TEST_MACRO((numRows != pqCovMatrix.numRows()) || (numCols != pqCovMatrix.numCols()),
                      env.worldRank(),
                      "uqComputeCovCorrMatricesBetweenVectorRvs()",
                      "inconsistent dimensions for covariance matrix");

  UQ_FATAL_TEST_MACRO((numRows != pqCorrMatrix.numRows()) || (numCols != pqCorrMatrix.numCols()),
                      env.worldRank(),
                      "uqComputeCorrelationBetweenVectorRvs()",
                      "inconsistent dimensions for correlation matrix");

  UQ_FATAL_TEST_MACRO((localNumSamples > paramRv.realizer().period()) || (localNumSamples > qoiRv.realizer().period()),
                      env.worldRank(),
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
