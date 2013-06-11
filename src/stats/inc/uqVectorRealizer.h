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

#ifndef __UQ_REALIZER_H__
#define __UQ_REALIZER_H__

#include <uqVectorSequence.h>
#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a probability density routine.
//*****************************************************

//*****************************************************
// Base class [R-00]
//*****************************************************
/*! \file uqVectorRealizerClass.h
 * \brief A templated class for sampling from vector RVs (holding probability density distributions).
 * 
 * \class uqBaseVectorRealizerClass
 * \brief A templated (base) class for handling sampling from vector RVs.
 *
 * A realizer is an object that, simply put, contains a realization() operation 
 * that returns a sample of a vector RV. This is the base class. QUESO also support
 * uniform, Gaussian, Beta, Gamma, Inverse Gamma and LogNormal realizers, as described
 * and implemented in the derived classes. */

template<class V, class M>
class uqBaseVectorRealizerClass {
public:
  
      //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer.*/
  uqBaseVectorRealizerClass(const char*                  prefix,
			    const uqVectorSetClass<V,M>& unifiedImageSet,
			    unsigned int                 subPeriod);

  //! Virtual destructor
  virtual ~uqBaseVectorRealizerClass();
  //@}
  
  //! @name Realization-related methods
  //@{
  //! Image set where the relizations lie.  Access to protected attribute m_unifiedImageSet.
  const   uqVectorSetClass<V,M>& unifiedImageSet()              const;
  
  //! Sub-period of the realization. Access to protected attribute m_subPeriod.
  unsigned int           subPeriod      ()              const;
  
  //! Performs a realization (sample) from a probability density function. See template specialization.
  virtual void                   realization    (V& nextValues) const = 0;
  //@}
  
protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
  const uqVectorSetClass<V,M>&  m_unifiedImageSet;
        unsigned int            m_subPeriod;
};
// Default constructor -----------------------------
template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  unsigned int                 subPeriod)
  :
  m_env            (unifiedImageSet.env()),
  m_prefix         ((std::string)(prefix)+"re_"),
  m_unifiedImageSet(unifiedImageSet),
  m_subPeriod      (subPeriod)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorRealizerClass<V,M>::constructor() [4]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorRealizerClass<V,M>::constructor() [4]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqBaseVectorRealizerClass<V,M>::~uqBaseVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
unsigned int
uqBaseVectorRealizerClass<V,M>::subPeriod() const
{
  return m_subPeriod;
}
//--------------------------------------------------
template<class V, class M>
const uqVectorSetClass<V,M>&
uqBaseVectorRealizerClass<V,M>::unifiedImageSet() const
{
  return m_unifiedImageSet;
}

//*****************************************************
// Generic class [R-01]
//*****************************************************
/*! 
 * \class uqGenericVectorRealizerClass
 * \brief A class for handling sampling from generic probability density distributions.
 *
 * A realizer is an object that, simply put, contains a realization() operation 
 * that returns a sample of a vector RV, or, particularly, a generic probability 
 * density distribution. This is the class that handles generic sampling, used, 
 * for example, to sample, posterior PDFs (the solution of a Bayesian problem).*/

template<class V, class M>
class uqGenericVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer, 
   * the sub period for the realizations and a pointer to a  generic routine. */ 
  uqGenericVectorRealizerClass(const char*                  prefix,
                               const uqVectorSetClass<V,M>& unifiedImageSet,
                               unsigned int                 subPeriod,
                               double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
                               const void* routineDataPtr);
 //! Destructor
  ~uqGenericVectorRealizerClass();
  //@}
  
  //! @name Realization-related methods
  //! Draws a realization.
  /*! This function draws a realization of \c this considering the generic routine \c m_routinePtr
   * and saves it in \c nextValues.*/
  void realization(V& nextValues) const;
  //@}

private:
  double (*m_routinePtr)(const void* routineDataPtr, V& nextParamValues);
  const void* m_routineDataPtr;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
};
// Default constructor -----------------------------
template<class V, class M>
uqGenericVectorRealizerClass<V,M>::uqGenericVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  unsigned int                 subPeriod,
  double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
  const void* routineDataPtr)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,subPeriod),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGenericVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGenericVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqGenericVectorRealizerClass<V,M>::~uqGenericVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqGenericVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  m_routinePtr(m_routineDataPtr,nextValues);
  return;
}

//*****************************************************
// Gaussian class [R-03]
//*****************************************************
/*! 
 * \class uqGaussianVectorRealizerClass
 * \brief A class for handling sampling from Gaussian probability density distributions.
 *
 * This class handles sampling from a Gaussian probability density distribution.*/

template<class V, class M>
class uqGaussianVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer, a
   * vector of mean values, \c lawExpVector, and a lower triangular matrix resulting from 
   * Cholesky decomposition of the covariance matrix, \c lowerCholLawCovMatrix.  */ 
  uqGaussianVectorRealizerClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& unifiedImageSet,
                                const V&                     lawExpVector, // vector of mean values
                                const M&                     lowerCholLawCovMatrix); // lower triangular matrix resulting from Cholesky decomposition of the covariance matrix

  //! Constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer, a
   * vector of mean values, \c lawExpVector, and a set of two matrices and one vector 
   * resulting from the Single Value Decomposition of the covariance matrix, \c matU, 
   * \c vecSsqrt and \c matVt.  */ 
  uqGaussianVectorRealizerClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& unifiedImageSet,
                                const V&                     lawExpVector, // vector of mean values
                                const M&                     matU,
                                const V&                     vecSsqrt,
                                const M&                     matVt);
  //! Destructor
  ~uqGaussianVectorRealizerClass();
  //@}

  //! @name Realization-related methods
  //@{
  //! Access to the vector of mean values and private attribute:  m_unifiedLawExpVector. 
  const V&   unifiedLawExpVector        ()              const;
  
  //! Access to the vector of variance values and private attribute:  m_unifiedLawVarVector. 
  const V&   unifiedLawVarVector        ()              const;
     
  //! Draws a realization.
  /*! This function draws a realization of a Gaussian distribution of mean \c m_unifiedLawExpVector 
   * and variance \c m_unifiedLawVarVector and saves it in \c nextValues.*/
  void realization                (V& nextValues) const;
  
  //! Updates the mean with the new value \c newLawExpVector.  
  void updateLawExpVector         (const V& newLawExpVector);
  
  //! Updates the lower triangular matrix from Cholesky decomposition of the covariance matrix to the new value \c newLowerCholLawCovMatrix.
  /*! The lower triangular matrix results resulting from a Cholesky decomposition of the 
   * covariance matrix. This routine deletes old expected values: m_lowerCholLawCovMatrix;
   *   m_matU, m_vecSsqrt, m_matVt.*/
  void updateLowerCholLawCovMatrix(const M& newLowerCholLawCovMatrix);
  
  //! Updates the SVD matrices from SVD decomposition of the covariance matrix to the new values: \c matU, \c vecSsqrt, and \c matVt.
  /*! The lower triangular matrix results resulting from a Cholesky decomposition of the 
   * covariance matrix. This routine deletes old expected values: m_lowerCholLawCovMatrix;
   *   m_matU, m_vecSsqrt, m_matVt. */
    void updateLowerCholLawCovMatrix(const M& matU,
				   const V& vecSsqrt,
				   const M& matVt);
  //@}
  
private:
  V* m_unifiedLawExpVector;
  V* m_unifiedLawVarVector;
  M* m_lowerCholLawCovMatrix;
  M* m_matU;
  V* m_vecSsqrt;
  M* m_matVt;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
};
// Constructor -------------------------------------
template<class V, class M>
uqGaussianVectorRealizerClass<V,M>::uqGaussianVectorRealizerClass(const char* prefix,
								  const uqVectorSetClass<V,M>& unifiedImageSet,
								  const V& lawExpVector,
								  const M& lowerCholLawCovMatrix)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(new M(lowerCholLawCovMatrix)),
  m_matU                 (NULL),
  m_vecSsqrt             (NULL),
  m_matVt                (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorRealizerClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  *m_unifiedLawExpVector = lawExpVector; // ????

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRealizerClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
uqGaussianVectorRealizerClass<V,M>::uqGaussianVectorRealizerClass(const char* prefix,
								  const uqVectorSetClass<V,M>& unifiedImageSet,
								  const V& lawExpVector,
								  const M& matU,
								  const V& vecSsqrt,
								  const M& matVt)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(NULL),
  m_matU                 (new M(matU)),
  m_vecSsqrt             (new V(vecSsqrt)),
  m_matVt                (new M(matVt))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorRealizerClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  *m_unifiedLawExpVector = lawExpVector; // ????

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorRealizerClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqGaussianVectorRealizerClass<V,M>::~uqGaussianVectorRealizerClass()
{
  delete m_matVt;
  delete m_vecSsqrt;
  delete m_matU;
  delete m_lowerCholLawCovMatrix;
  delete m_unifiedLawVarVector;
  delete m_unifiedLawExpVector;
}
// Realization-related methods----------------------
template <class V, class M>
const V&
uqGaussianVectorRealizerClass<V,M>::unifiedLawExpVector() const
{
  return *m_unifiedLawExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
uqGaussianVectorRealizerClass<V,M>::unifiedLawVarVector() const
{
  return *m_unifiedLawVarVector;
}
//--------------------------------------------------
template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  V iidGaussianVector(m_unifiedImageSet.vectorSpace().zeroVector());

  bool outOfSupport = true;
  do {
    iidGaussianVector.cwSetGaussian(0.0, 1.0);

    if (m_lowerCholLawCovMatrix) {
      nextValues = (*m_unifiedLawExpVector) + (*m_lowerCholLawCovMatrix)*iidGaussianVector;
    }
    else if (m_matU && m_vecSsqrt && m_matVt) {
      nextValues = (*m_unifiedLawExpVector) + (*m_matU)*( (*m_vecSsqrt) * ((*m_matVt)*iidGaussianVector) );
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqGaussianVectorRealizerClass<V,M>::realization()",
                          "inconsistent internal state");
    }

    outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));
  } while (outOfSupport); // prudenci 2011-Oct-04

  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_unifiedLawExpVector;

  m_unifiedLawExpVector = new V(newLawExpVector);
 
  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::updateLowerCholLawCovMatrix(const M& newLowerCholLawCovMatrix)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lowerCholLawCovMatrix;
  delete m_matU;
  delete m_vecSsqrt;
  delete m_matVt;

  m_lowerCholLawCovMatrix = new M(newLowerCholLawCovMatrix);
  m_matU                  = NULL;
  m_vecSsqrt              = NULL;
  m_matVt                 = NULL;

  return;
}
//--------------------------------------------------
template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::updateLowerCholLawCovMatrix(
  const M& matU,
  const V& vecSsqrt,
  const M& matVt)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lowerCholLawCovMatrix;
  delete m_matU;
  delete m_vecSsqrt;
  delete m_matVt;

  m_lowerCholLawCovMatrix = NULL;
  m_matU                  = new M(matU);
  m_vecSsqrt              = new V(vecSsqrt);
  m_matVt                 = new M(matVt);

  return;
}

//*****************************************************
// Sequential class [R-nn]
//*****************************************************
/*! 
 * \class uqSequentialVectorRealizerClass
 * \brief A class for handling sequential draws (sampling) from probability density distributions.
 *
 * This class handles sequential sampling (it returns the next value of the chain) from a 
 * probability density distribution.*/

template<class V, class M>
class uqSequentialVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
  //!@name Constructor/Destructor methods
  //! Default constructor.
  uqSequentialVectorRealizerClass(const char*                           prefix,
                                  const uqBaseVectorSequenceClass<V,M>& chain);
  
  //! Destructor.
  ~uqSequentialVectorRealizerClass();
  //@}
 
  //!@name Sampling-related methdos
  //! Returns the unified mean vector; access to private attribute m_unifiedSampleExpVector.
  const V&   unifiedSampleExpVector()              const;
  
  //! Returns the unified variance vector; access to private attribute m_unifiedSampleVarVector.
  const V&   unifiedSampleVarVector()              const;
  
  //! Draws the next value from this chain (\c m_chain) and saves it in \c nextValues
  void realization           (V& nextValues) const;
  //@}

private:
  const   uqBaseVectorSequenceClass<V,M>& m_chain;
  mutable unsigned int                    m_currentChainPos;
          V*                              m_unifiedSampleExpVector;
          V*                              m_unifiedSampleVarVector;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
};
// Constructor -------------------------------------
template<class V, class M>
uqSequentialVectorRealizerClass<V,M>::uqSequentialVectorRealizerClass(
  const char*                           prefix,
  const uqBaseVectorSequenceClass<V,M>& chain)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"seq").c_str(),chain.unifiedBoxPlain(),chain.subSequenceSize()),
  m_chain                 (chain),
  m_currentChainPos       (0),
  m_unifiedSampleExpVector(new V(chain.unifiedMeanPlain()          )), // IMPORTANT
  m_unifiedSampleVarVector(new V(chain.unifiedSampleVariancePlain()))  // IMPORTANT
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqSequentialVectorRealizerClass<V,M>::constructor()"
                            << ": m_chain.subSequenceSize() = " << m_chain.subSequenceSize()
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqSequentialVectorRealizerClass<V,M>::~uqSequentialVectorRealizerClass()
{
  delete m_unifiedSampleVarVector;
  delete m_unifiedSampleExpVector;
}
// Realization-related methods----------------------
template<class V, class M>
void
uqSequentialVectorRealizerClass<V,M>::realization(V& nextParamValues) const
{
  m_chain.getPositionValues(m_currentChainPos++,nextParamValues);
  if (m_currentChainPos >= m_subPeriod) m_currentChainPos = 0;

  return;
}
//--------------------------------------------------
template <class V, class M>
const V&
uqSequentialVectorRealizerClass<V,M>::unifiedSampleExpVector() const
{
  return *m_unifiedSampleExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
uqSequentialVectorRealizerClass<V,M>::unifiedSampleVarVector() const
{
  return *m_unifiedSampleVarVector;
}

//*****************************************************
// Uniform class [R-04]
//*****************************************************
/*! 
 * \class uqUniformVectorRealizerClass
 * \brief A class for handling sampling from a Uniform probability density distribution.
 *
 * This class handles sampling from a uniform probability density distribution.*/

template<class V, class M>
class uqUniformVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer.  */
  uqUniformVectorRealizerClass(const char*                  prefix,
                               const uqVectorSetClass<V,M>& unifiedImageSet);
  //! Destructor
  ~uqUniformVectorRealizerClass();
  //@}
  
  //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a uniform distribution and saves it in \c nextValues. It 
   * internally finds the minimum and the maximum values of the distribution.
   */  
  void realization(V& nextValues) const;
  //@}
  
private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
};
// Constructor -------------------------------------
template<class V, class M>
uqUniformVectorRealizerClass<V,M>::uqUniformVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqUniformVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqUniformVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqUniformVectorRealizerClass<V,M>::~uqUniformVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqUniformVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&m_unifiedImageSet);

  UQ_FATAL_TEST_MACRO(imageBox == NULL,
                      m_env.worldRank(),
                      "uqUniformVectorRealizerClass<V,M>::realization()",
                      "only box images are supported right now");
  
  nextValues.cwSetUniform(imageBox->minValues(),imageBox->maxValues());
  return;
}

//*****************************************************
// Beta class [R-05]
//*****************************************************
/*! 
 * \class uqBetaVectorRealizerClass
 * \brief A class for handling sampling from a Beta probability density distribution.
 *
 * This class handles sampling from a Beta probability density distribution, of 
 * parameters \c alpha and \c beta.*/

template<class V, class M>
class uqBetaVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
   //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, and the
   * Beta distribution parameters \c alpha and \c beta, which are assigned to private attributes
   * m_alpha and m_beta.  */
  uqBetaVectorRealizerClass(const char*                  prefix,
                            const uqVectorSetClass<V,M>& unifiedImageSet,
                            const V&                     alpha,
                            const V&                     beta);
 
  //! Destructor
  ~uqBetaVectorRealizerClass();
  //@}

    //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a Beta distribution and saves it in \c nextValues. 
   * It internally checks whether the image set, where the realization should be drawn, belongs 
   * to the interval (0, 1] - which is the range where Beta distribution is defined over. */ 
  void realization(V& nextValues) const;
  //@}
private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;

  V m_alpha;
  V m_beta;
};
// Constructor -------------------------------------
template<class V, class M>
uqBetaVectorRealizerClass<V,M>::uqBetaVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBetaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBetaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqBetaVectorRealizerClass<V,M>::~uqBetaVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqBetaVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  nextValues.cwSetBeta(m_alpha,m_beta);
  return;
}

//*****************************************************
// Gamma class [R-06]
//*****************************************************
/*! 
 * \class uqGammaVectorRealizerClass
 * \brief A class for handling sampling from a Gamma probability density distribution.
 *
 * This class handles sampling from a Gamma probability density distribution, of 
 * parameters \c a and \c b.*/
template<class V, class M>
class uqGammaVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, and the
   * Gamma distribution parameters \c a and \c b, which are assigned to private attributes
   * m_alpha and m_beta.  */
  uqGammaVectorRealizerClass(const char*                  prefix,
                             const uqVectorSetClass<V,M>& unifiedImageSet,
                             const V&                     a,
                             const V&                     b);
  
  //! Destructor
 ~uqGammaVectorRealizerClass();
  //@}
 
   //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a Gamma distribution and saves it in \c nextValues. 
   * It internally checks whether the image set, where the realization should be drawn, belongs 
   * to the interval (0, infinity) - which is the range where Gamma distribution is defined over. */ 
  void realization(V& nextValues) const;

private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;

  V m_a;
  V m_b;
};
// Constructor -------------------------------------
template<class V, class M>
uqGammaVectorRealizerClass<V,M>::uqGammaVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  const V&                     a,
  const V&                     b)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_a(a),
  m_b(b)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGammaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGammaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqGammaVectorRealizerClass<V,M>::~uqGammaVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqGammaVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  // nextValues.cwSetGamma(m_a,m_b);

// begin kemelli 2013-April-22 :
  const uqBoxSubsetClass<V,M>* imageBox = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&this->m_unifiedImageSet);
//  double biggerOfMaxValues = imageBox->maxValues().getMaxValue();
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();	
  double smallerOfMinValues = imageBox->minValues().getMinValue();
	
 // Gamma dist belongs to (0,inf)		
 if( smallerOfMinValues < 0 ) //(biggerOfMinValues < 0) || 
 {		
   std::cerr << "In uqGammaVectorRealizerClass<V,M>::realization()\n" 
			 << "Gamma distribution is only defined in (0, infinity).\n"
			 << "The data provided is: \n"
			 << *imageBox 
   			 << "Sampling will not cover all inteval.\n"   
   			 << std::endl;


    UQ_FATAL_TEST_MACRO(smallerOfMaxValues < 0,
                      m_env.worldRank(),
                      "uqGammaVectorRealizerClass<V,M>::realization()",
                      "invalid input: Gamma distribution is only defined in (0, infinity), and min(m_maxValues)<0. ");      
                      
 //  UQ_FATAL_TEST_MACRO(biggerOfMaxValues < 0,
 //                     m_env.worldRank(),
 //                     "uqGammaVectorRealizerClass<V,M>::realization()",
 //                     "invalid input: Gamma distribution is only defined in (0, infinity). ");            
              
 }	

  // end kemelli 2013-April-22 
  
  // begin kemelli 2013-April-18 : 
  bool outOfSupport = true;
  do {

	nextValues.cwSetGamma(m_a,m_b);
	outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));

  } while (outOfSupport); 

  // end kemelli 2013-April-18 :
  
  return;
}

//*****************************************************
// InverseGamma class [R-07]
//*****************************************************
/*! 
 * \class uqInverseGammaVectorRealizerClass
 * \brief A class for handling sampling from an Inverse Gamma probability density distribution.
 *
 * This class handles sampling from an Inverse Gamma probability density distribution, of 
 * parameters \c alpha and \c beta.*/

template<class V, class M>
class uqInverseGammaVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, and the
   * Beta distribution parameters \c a and \c b, which are assigned to private attributes
   * m_alpha and m_beta.  */
  uqInverseGammaVectorRealizerClass(const char*                  prefix,
                                    const uqVectorSetClass<V,M>& unifiedImageSet,
                                    const V&                     alpha,
                                    const V&                     beta);
  //! Destructor
  ~uqInverseGammaVectorRealizerClass();
  //@}

     //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of an Inverse Gamma distribution and saves it in \c nextValues. 
   * It internally checks whether the image set, where the realization should be drawn, belongs 
   * to the interval (0, infinity) - which is the range where Gamma distribution is defined over. */ 
  void realization(V& nextValues) const;
  //@}

private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;

  V m_alpha;
  V m_beta;
};
// Constructor -------------------------------------
template<class V, class M>
uqInverseGammaVectorRealizerClass<V,M>::uqInverseGammaVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqInverseGammaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqInverseGammaVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqInverseGammaVectorRealizerClass<V,M>::~uqInverseGammaVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqInverseGammaVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  nextValues.cwSetInverseGamma(m_alpha,m_beta);
  return;
}

//*****************************************************
// Wigner class [R-09]
//*****************************************************
/*!
 * \class uqWignerVectorRealizerClass
 * \brief A class for handling sampling from a Wigner probability density distribution.
 *
 * This class \b will handle sampling from an Wigner probability density distribution, with a
 * given center position and a radius.
 * 
 * \todo: The metodo uqWignerVectorRealizerClass:realization() is not yet available, 
 * thus this class does  nothing. */

template<class V, class M>
class uqWignerVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, the
   * center position \c centerPos, and a radius \c radius.*/  
  uqWignerVectorRealizerClass(const char*                  prefix,
                              const uqVectorSetClass<V,M>& unifiedImageSet,
                              const V&                     centerPos,
                              double                       radius);
  
  //! Destructor
 ~uqWignerVectorRealizerClass();
 //@}
  
  //! @name Realization-related methods
  //@{
  //! TODO: Draws a realization.
  /*! \todo: implement and explain me!*/
  void realization(V& nextValues) const;
  //@}
 
private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
  V*     m_centerPos;
  double m_radius;
};
// Constructor -------------------------------------
template<class V, class M>
uqWignerVectorRealizerClass<V,M>::uqWignerVectorRealizerClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& unifiedImageSet,
  const V&                     centerPos,
  double                       radius)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_centerPos(new V(centerPos)),
  m_radius   (radius)    
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqWignerVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_radius <= 0.,
                      m_env.worldRank(),
                      "uqWignerVectorRealizerClass<V,M>::constructor()",
                      "invalid radius");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqWignerVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqWignerVectorRealizerClass<V,M>::~uqWignerVectorRealizerClass()
{
  delete m_centerPos;
}
// -------------------------------------------------
// TODO: implement me, please!!!!
template<class V, class M>
void
uqWignerVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqWignerVectorRealizerClass<V,M>::realization()",
                      "not implemented yet");
  
  nextValues.cwSet(0.);
  return;
}

//*****************************************************
// LogNormal class [R-10]
//*****************************************************
/*! 
 * \class uqLogNormalVectorRealizerClass
 * \brief A class for handling sampling from a Log-Normal probability density distribution.
 *
 * This class handles sampling from a Log-Normal probability density distribution, of 
 * mean and variance given.*/

template<class V, class M>
class uqLogNormalVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the image set, a  vector of 
   * mean values, \c lawExpVector, and a lower triangular matrix resulting from Cholesky 
   * decomposition of the covariance matrix, \c lowerCholLawCovMatrix.  */ 
  uqLogNormalVectorRealizerClass(const char*                  prefix,
                                 const uqVectorSetClass<V,M>& unifiedImageSet,
                                 const V&                     lawExpVector, // vector of mean values
                                 const M&                     lowerCholLawCovMatrix); // lower triangular matrix resulting from Cholesky decomposition of the covariance matrix

  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the image set, a vector of
   * mean values, \c lawExpVector, and a set of two matrices and one vector resulting from the 
   * Single Value Decomposition of the covariance matrix, \c matU, \c vecSsqrt and \c matVt.*/ 
  uqLogNormalVectorRealizerClass(const char*                  prefix,
                                 const uqVectorSetClass<V,M>& unifiedImageSet,
                                 const V&                     lawExpVector, // vector of mean values
                                 const M&                     matU,
                                 const V&                     vecSsqrt,
                                 const M&                     matVt);
  //! Destructor
  ~uqLogNormalVectorRealizerClass();
  //@}
  
  //! @name Realization-related methods
  //@{
  //! Access to the vector of mean values and private attribute:  m_unifiedLawExpVector. 
  const V&   unifiedLawExpVector ()              const;
  
  //! Access to the vector of variance values and private attribute:  m_unifiedLawVarVector. 
  const V&   unifiedLawVarVector ()              const;
  
  //! Draws a realization.
  /*! This function draws a realization of a LogNormal distribution and saves it in \c nextValues.*/
  void realization         (V& nextValues) const;
  //@}
    
private:
  V* m_unifiedLawExpVector;
  V* m_unifiedLawVarVector;
  M* m_lowerCholLawCovMatrix;
  M* m_matU;
  V* m_vecSsqrt;
  M* m_matVt;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;
};
// Constructor -------------------------------------
template<class V, class M>
uqLogNormalVectorRealizerClass<V,M>::uqLogNormalVectorRealizerClass(const char* prefix,
                                                                    const uqVectorSetClass<V,M>& unifiedImageSet,
                                                                    const V& lawExpVector,
                                                                    const M& lowerCholLawCovMatrix)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(new M(lowerCholLawCovMatrix)),
  m_matU                 (NULL),
  m_vecSsqrt             (NULL),
  m_matVt                (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalVectorRealizerClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalVectorRealizerClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------								  
template<class V, class M>
uqLogNormalVectorRealizerClass<V,M>::uqLogNormalVectorRealizerClass(const char* prefix,
                                                                    const uqVectorSetClass<V,M>& unifiedImageSet,
                                                                    const V& lawExpVector,
                                                                    const M& matU,
                                                                    const V& vecSsqrt,
                                                                    const M& matVt)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(NULL),
  m_matU                 (new M(matU)),
  m_vecSsqrt             (new V(vecSsqrt)),
  m_matVt                (new M(matVt))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalVectorRealizerClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalVectorRealizerClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------								  
template<class V, class M>
uqLogNormalVectorRealizerClass<V,M>::~uqLogNormalVectorRealizerClass()
{
  delete m_matVt;
  delete m_vecSsqrt;
  delete m_matU;
  delete m_lowerCholLawCovMatrix;
  delete m_unifiedLawVarVector;
  delete m_unifiedLawExpVector;
}
// Realization-related methods----------------------
template <class V, class M>
const V&
uqLogNormalVectorRealizerClass<V,M>::unifiedLawExpVector() const
{
  return *m_unifiedLawExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
uqLogNormalVectorRealizerClass<V,M>::unifiedLawVarVector() const
{
  return *m_unifiedLawVarVector;
}
//--------------------------------------------------
template<class V, class M>
void
uqLogNormalVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  V iidGaussianVector(m_unifiedImageSet.vectorSpace().zeroVector());

  bool outOfSupport = true;
  do {
    iidGaussianVector.cwSetGaussian(0.0, 1.0);

    if (m_lowerCholLawCovMatrix) {
      nextValues = (*m_unifiedLawExpVector) + (*m_lowerCholLawCovMatrix)*iidGaussianVector;
    }
    else if (m_matU && m_vecSsqrt && m_matVt) {
      nextValues = (*m_unifiedLawExpVector) + (*m_matU)*( (*m_vecSsqrt) * ((*m_matVt)*iidGaussianVector) );
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqLogNormalVectorRealizerClass<V,M>::realization()",
                          "inconsistent internal state");
    }

    for (unsigned int i = 0; i < nextValues.sizeLocal(); ++i) {
      nextValues[i] = std::exp(nextValues[i]);
    }

    outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));
  } while (outOfSupport); // prudenci 2011-Oct-04

  return;
}

//*****************************************************
// Concatenated class [R-11]
//*****************************************************
/*!
 * \class uqConcatenatedVectorRealizerClass
 * \brief A class for handling sampling from concatenated probability density distributions.
 * 
 * This class allows the user draw samples from concatenated probability density distributions (two 
 * or more distincts probability distributions has(ve) been concatened into one single vector RV). 
 * This class used, for instance, to draw realization of concatenate priors from two or more RVs, 
 * where one of them has a uniform distribution whereas the other one(s) has a Gaussian distribution. */

template<class V, class M>
class uqConcatenatedVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Concatenates two RVs: \c rv1 and \c rv2 into one vector RV, given a prefix and the image set of the vector RV.*/
  uqConcatenatedVectorRealizerClass(const char*                           prefix,
                                    const uqBaseVectorRealizerClass<V,M>& realizer1,
                                    const uqBaseVectorRealizerClass<V,M>& realizer2,
                                    const uqVectorSetClass<V,M>&          unifiedImageSet);
  //! Constructor
  /*! Concatenates a sequence of RVs, given by: <c> std::vector<const uqBaseVectorRVClass<V,M>* >& rvs </c>
   * into one singke vector RV, given a prefix and the image set of the resulting vector RV.*/
  uqConcatenatedVectorRealizerClass(const char*                                                prefix,
                                    const std::vector<const uqBaseVectorRealizerClass<V,M>* >& realizers,
                                    unsigned int                                               minPeriod,
                                    const uqVectorSetClass<V,M>&                               unifiedImageSet);
  //! Destructor
  ~uqConcatenatedVectorRealizerClass();
  //@}
  
  //! @name Realization-related methods
  //@{
  void realization(V& nextValues) const;
  //@}

private:
  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_unifiedImageSet;
  using uqBaseVectorRealizerClass<V,M>::m_subPeriod;

  std::vector<const uqBaseVectorRealizerClass<V,M>* > m_realizers;
};
// Constructor -------------------------------------
template<class V, class M>
uqConcatenatedVectorRealizerClass<V,M>::uqConcatenatedVectorRealizerClass(
  const char*                           prefix,
  const uqBaseVectorRealizerClass<V,M>& realizer1,
  const uqBaseVectorRealizerClass<V,M>& realizer2,
  const uqVectorSetClass<V,M>&          unifiedImageSet)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::min(realizer1.subPeriod(),realizer2.subPeriod()) ), // 2011/Oct/02
  m_realizers(2,(const uqBaseVectorRealizerClass<V,M>*) NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_realizers[0] = &realizer1;
  m_realizers[1] = &realizer2;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedVectorRealizerClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
uqConcatenatedVectorRealizerClass<V,M>::uqConcatenatedVectorRealizerClass(
  const char*                                                prefix,
  const std::vector<const uqBaseVectorRealizerClass<V,M>* >& realizers,
  unsigned int                                               minPeriod,
  const uqVectorSetClass<V,M>&                               unifiedImageSet)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gen").c_str(),unifiedImageSet,minPeriod),
  m_realizers(realizers.size(),(const uqBaseVectorRealizerClass<V,M>*) NULL)
{
  for (unsigned int i = 0; i < m_realizers.size(); ++i) {
    m_realizers[i] = realizers[i];
  }
}
// Destructor --------------------------------------
template<class V, class M>
uqConcatenatedVectorRealizerClass<V,M>::~uqConcatenatedVectorRealizerClass()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
uqConcatenatedVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  std::vector<V*> vecs(m_realizers.size(),(V*)NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_realizers[i]->unifiedImageSet().vectorSpace().zeroVector());
    //std::cout << "In uqConcatenatedVectorRealizerClass<V,M>::realization: v[i]->sizeLocal() = " << v[i]->sizeLocal() << std::endl;
    m_realizers[i]->realization(*(vecs[i]));
  }

  //std::cout << "In uqConcatenatedVectorRealizerClass<V,M>::realization: nextValues.sizeLocal() = " << nextValues.sizeLocal() << std::endl;
  std::vector<const V*> constVecs(m_realizers.size(),(V*)NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    constVecs[i] = vecs[i];
  }
  nextValues.cwSetConcatenated(constVecs);
  //std::cout << "In uqConcatenatedVectorRealizerClass<V,M>::realization: succeeded" << std::endl;

  for (unsigned int i = 0; i < vecs.size(); ++i) {
    delete vecs[i];
  }

  return;
}
#endif // __UQ_REALIZER_H__
