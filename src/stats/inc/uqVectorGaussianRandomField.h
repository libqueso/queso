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

#ifndef __UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H__
#define __UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H__

#include <uqMatrixCovarianceFunction.h>
#include <uqVectorFunction.h>
#include <uqVectorRV.h>

namespace QUESO {

/*! 
 * \file uqVectorGaussianRandomField.h
 * \brief A class for handling Gaussian random fields (GRF).
 * 
 * \class VectorGaussianRandomField
 * \brief A class for handling vector Gaussian random fields (GRF).
 *
 * This class implements a vector Gaussian random field (GRF); i.e. a random field involving 
 * vector Gaussian probability density functions (PDFs) of the variables. */

template <class P_V, class P_M, class Q_V, class Q_M>
class VectorGaussianRandomField
{
 public:
   //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! Constructs a new object, given a prefix, an index set, and both a mean and a 
   * covariance function. This method deletes the previous saved positions. */ 
  VectorGaussianRandomField(const char*                                                 prefix,
                                   const VectorSet<P_V,P_M>&                            indexSet,
                                   const VectorSet<Q_V,Q_M>&                            imageSetPerIndex,
                                   const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&           meanFunction, 
                                   const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& covarianceFunction);
  
  //! TODO: Copy constructor.
  /*! \todo: implement me!*/
  VectorGaussianRandomField(const VectorGaussianRandomField&                     obj);
 
  //! Destructor.
  ~VectorGaussianRandomField();
  //@}
 
  //! @name Set methods
  //@{
  //! TODO: Assignment operator; it copies \c rhs to \c this.
  /*! \todo: implement me!*/
  VectorGaussianRandomField& operator=(const VectorGaussianRandomField& rhs);
  //@}

  //! @name Math methods
  //@{
  //! Index set; access to protected attribute m_indexSet.
  const VectorSet<P_V,P_M>&                   indexSet          () const;
  
  //! Mean function; access to protected attribute m_meanFunction.
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&  meanFunction      () const;
  
  //! Covariance function; access to protected attribute m_covarianceFunction.
  const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& covarianceFunction() const;
  
  //! Function that samples from a Gaussian PDF.
  /*! Given the field positions, this method performs a number of tests, calculates the mean vector,
   * the covariance matrix and then it samples from a Gaussian random vector as 
   * many positions as required.*/
  void                                              sampleFunction(const std::vector<P_V*>& fieldPositions, Q_V& sampleValues);
  //@}
protected:
  //! Copy method.
  void                                              copy          (const VectorGaussianRandomField& src);

  //! Environment.
  const BaseEnvironment&                     m_env;
  
  //! Prefix.
  std::string                                       m_prefix;
 
  //! Index set.	
  const VectorSet<P_V,P_M>&                  m_indexSet;
    
  //! Image set of the RV, per index.
  const VectorSet<Q_V,Q_M>&                  m_imageSetPerIndex;
  
  //! Mean function.
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& m_meanFunction;
  
  //! Covariance function.
  const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& m_covarianceFunction;

  //! Saved positions.
  std::vector<P_V*>                                 m_savedPositions;
  
  //! Image set of the RV.
  VectorSpace<Q_V,Q_M>*                      m_savedRvImageSpace;
  
  //! Vector of the mean value of the RV.
  Q_V*                                              m_savedRvLawExpVector;
  
  //! Covariance matrix of the RV.
  Q_M*                                              m_savedRvLawCovMatrix;
  
  //! My RV.
  GaussianVectorRV<Q_V,Q_M>*                 m_savedRv;
};

// Default constructor -----------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::VectorGaussianRandomField(
  const char*                                                 prefix,
  const VectorSet<P_V,P_M>&                            indexSet,
  const VectorSet<Q_V,Q_M>&                            imageSetPerIndex,
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&           meanFunction, 
  const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& covarianceFunction)
  :
  m_env                (indexSet.env()),
  m_prefix             ((std::string)(prefix)+"grf_"),
  m_indexSet           (indexSet),
  m_imageSetPerIndex   (imageSetPerIndex),
  m_meanFunction       (meanFunction),
  m_covarianceFunction (covarianceFunction),
  m_savedRvImageSpace  (NULL),
  m_savedRvLawExpVector(NULL),
  m_savedRvLawCovMatrix(NULL),
  m_savedRv            (NULL)
{
  m_savedPositions.clear();
}
// Destructor ---------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::~VectorGaussianRandomField()
{
}
// Math methods -------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
const VectorSet<P_V,P_M>&
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::indexSet() const
{
  return m_indexSet;
}
// --------------------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::meanFunction() const
{
  return m_meanFunction;
}
// --------------------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>&
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::covarianceFunction() const
{
  return m_covarianceFunction;
}
// --------------------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
void
VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction(const std::vector<P_V*>& fieldPositions, Q_V& sampleValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "Entering VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(( sampleValues.sizeLocal() % fieldPositions.size() ) != 0,
                      m_env.fullRank(),
                      "VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()",
                      "input data is not multiple of each other");

  unsigned int numberOfImageValuesPerIndex = sampleValues.sizeLocal()/fieldPositions.size();
  
  UQ_FATAL_TEST_MACRO(numberOfImageValuesPerIndex != m_imageSetPerIndex.vectorSpace().dimLocal(),
                      m_env.fullRank(),
                      "VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()",
                      "invalid input data dimension");

  if ((m_savedPositions.size() == 0   ) &&
      (m_savedRvImageSpace     == NULL) &&
      (m_savedRvLawExpVector   == NULL) &&
      (m_savedRvLawCovMatrix   == NULL) &&
      (m_savedRv               == NULL)) {
    // Ok
  }
  else if ((m_savedPositions.size() != 0   ) &&
           (m_savedRvImageSpace     != NULL) &&
           (m_savedRvLawExpVector   != NULL) &&
           (m_savedRvLawCovMatrix   != NULL) &&
           (m_savedRv               != NULL)) {
    // Ok
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.fullRank(),
                        "VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()",
                        "invalid combination of pointer values");
  }

  unsigned int numberOfPositions = fieldPositions.size();
  bool instantiate = true;
  if (m_savedPositions.size() == numberOfPositions) {
    bool allPositionsAreEqual = true;
    for (unsigned int i = 0; i < numberOfPositions; ++i) {
      UQ_FATAL_TEST_MACRO(m_savedPositions[i] == NULL,
                          m_env.fullRank(),
                          "VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()",
                          "m_savedPositions[i] should not be NULL");
      if ((m_savedPositions[i]->sizeLocal() == fieldPositions[i]->sizeLocal()) &&
          (*(m_savedPositions[i])           == *(fieldPositions[i])          )) {
        // Ok
      }
      else {
        allPositionsAreEqual = false;
        break;
      }
    } // for i
    instantiate = !allPositionsAreEqual;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                            << ": numberOfPositions = " << numberOfPositions
                            << ", instantiate = "       << instantiate
                            << std::endl;
  }

  if (instantiate) {
    delete m_savedRv;
    delete m_savedRvLawCovMatrix;
    delete m_savedRvLawExpVector;
    delete m_savedRvImageSpace;
    for (unsigned int i = 0; i < m_savedPositions.size(); ++i) {
      delete m_savedPositions[i];
    }
    m_savedPositions.clear();

    // Set m_savedPositions
    m_savedPositions.resize(numberOfPositions,NULL);
    for (unsigned int i = 0; i < m_savedPositions.size(); ++i) {
      m_savedPositions[i] = new P_V(*(fieldPositions[i]));
    }

    // Set m_savedRvImageSpace
    m_savedRvImageSpace = new VectorSpace<Q_V,Q_M>(m_env, "grf_", numberOfPositions*numberOfImageValuesPerIndex, NULL);

    // Set m_savedRvLawExpVector
    Q_V tmpVec(m_imageSetPerIndex.vectorSpace().zeroVector());
    m_savedRvLawExpVector = new Q_V(m_savedRvImageSpace->zeroVector());
    for (unsigned int i = 0; i < numberOfPositions; ++i) {
      m_meanFunction.compute(*(fieldPositions[i]),NULL,tmpVec,NULL,NULL,NULL);
      for (unsigned int j = 0; j < numberOfImageValuesPerIndex; ++j) {
        (*m_savedRvLawExpVector)[i*numberOfImageValuesPerIndex + j] = tmpVec[j];
      }
    }

    // Set m_savedRvLawCovMatrix
    Q_M tmpMat(m_imageSetPerIndex.vectorSpace().zeroVector());
    m_savedRvLawCovMatrix = new Q_M(m_savedRvImageSpace->zeroVector());
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                              << ": m_savedRvLawCovMatrix order = " << m_savedRvLawCovMatrix->numCols()
                              << ", numberOfPositions = "           << numberOfPositions
                              << ", tmpMat order = "                << tmpMat.numCols()
                              << ", numberOfImageValuesPerIndex = " << numberOfImageValuesPerIndex
                              << std::endl;
    }
    for (unsigned int i = 0; i < numberOfPositions; ++i) {
      for (unsigned int j = 0; j < numberOfPositions; ++j) {
        m_covarianceFunction.covMatrix(*(fieldPositions[i]),*(fieldPositions[j]),tmpMat);
#if 1
        Q_M testMat(tmpMat);
        if (testMat.chol() != 0) {
          *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                                  << ": i = " << i
                                  << ", j = " << j
                                  << ", *(fieldPositions[i]) = " << *(fieldPositions[i])
                                  << ", *(fieldPositions[j]) = " << *(fieldPositions[j])
                                  << ", tmpMat = "               << tmpMat
                                  << ", testMat = "              << testMat
                                  << ", tmpMat is not positive definite"
                                  << std::endl;
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()",
                              "tmpMat is not positive definite");
        }
#endif
        for (unsigned int k1 = 0; k1 < numberOfImageValuesPerIndex; ++k1) {
          for (unsigned int k2 = 0; k2 < numberOfImageValuesPerIndex; ++k2) {
            unsigned int tmpI = i*numberOfImageValuesPerIndex + k1;
            unsigned int tmpJ = j*numberOfImageValuesPerIndex + k2;
            (*m_savedRvLawCovMatrix)(tmpI,tmpJ) = tmpMat(k1,k2);
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
              *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                                      << ": i = " << i
                                      << ", j = " << j
                                      << ", k1 = " << k1
                                      << ", k2 = " << k2
                                      << ", tmpI = " << tmpI
                                      << ", tmpJ = " << tmpJ
                                      << ", *(fieldPositions[i]) = " << *(fieldPositions[i])
                                      << ", *(fieldPositions[j]) = " << *(fieldPositions[j])
                                      << ", (*m_savedRvLawCovMatrix)(tmpI,tmpJ) = " << (*m_savedRvLawCovMatrix)(tmpI,tmpJ)
                                      << std::endl;
            }
          }
        }
      }
    }

    // Set m_savedRv
    m_savedRv = new GaussianVectorRV<Q_V,Q_M>("grf_",
                                                     *m_savedRvImageSpace,
                                                     *m_savedRvLawExpVector,
                                                     *m_savedRvLawCovMatrix);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                              << ": just instantiated Gaussian RV"
                              << "\n *m_savedRvLawExpVector = " << *m_savedRvLawExpVector
                              << "\n *m_savedRvLawCovMatrix = " << *m_savedRvLawCovMatrix
                              << std::endl;
      for (unsigned int i = 0; i < numberOfPositions; ++i) {
        *m_env.subDisplayFile() << " *(m_savedPositions[" << i
                                << "]) = "                << *(m_savedPositions[i])
                                << std::endl;
      }
    }
  } // if (instantiate)

  // Generate sample function
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                            << ": about to realize sample values"
                            << std::endl;
  }
  m_savedRv->realizer().realization(sampleValues);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                            << ": just realized sample values"
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "Leaving VectorGaussianRandomField<P_V,P_M,Q_V,Q_M>::sampleFunction()"
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO

#endif // __UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H__ 
