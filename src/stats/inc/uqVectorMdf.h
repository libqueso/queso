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

#ifndef __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__
#define __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__

#include <uqArrayOfOneDGrids.h>
#include <uqArrayOfOneDTables.h>
#include <uqEnvironment.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Classes to accommodate a marginal density function
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
/*! \file uqVectorMdf.h
 * \brief Classes to accommodate a marginal density function of a vector RV.
 * 
 * \class BaseVectorMdfClass
 * \brief A templated (base) class for handling MDFs of vector functions.
 * 
 * To obtain the marginal distribution over a subset of multivariate (vector) RVs, one only 
 * needs to drop the irrelevant variables (the variables that one wants to marginalize out). 
 * If \b X is a Random Vector which contains the continuous random variables \f$ X_1, X_2, ..., X_n \f$. 
 * Then each \f$ X_i \f$ is a continuous random variable with its own Probability Distribution called 
 * the marginal distribution of \f$ X_i \f$.
 * This class handles MDFs of a vector RV.*/

template<class V, class M>
class BaseVectorMdfClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and domain set of the MDF.*/
  BaseVectorMdfClass(const char*                  prefix,
                                const VectorSetClass<V,M>& domainSet);
  //! Virtual destructor.
  virtual ~BaseVectorMdfClass();
  //@}

  //! @name Mathematical methods
  //@{ 
  //! Returns the domain set; access to protected attribute m_domainSet.
  const   VectorSetClass<V,M>& domainSet() const;
  
  //! Finds the value of the vector MDF at each element of \c paramValue, and saves it in \c mdfVec. See template specialization.  
  virtual void                   values   (const V& paramValues,
                                                 V& mdfVec)  const = 0;
  //@}
  
  //! @name I/O methods
  //@{ 
  //! Prints the vector MDF. See template specialization. 
  virtual void                   print    (std::ostream& os) const = 0;
  //@}

protected:

  const   BaseEnvironmentClass& m_env;
          std::string             m_prefix;
  const   VectorSetClass<V,M>&  m_domainSet;
};
// Default constructor -----------------------------
template<class V, class M>
BaseVectorMdfClass<V,M>::BaseVectorMdfClass(
  const char*                  prefix,
  const VectorSetClass<V,M>& domainSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"Mdf_"),
  m_domainSet(domainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering BaseVectorMdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving BaseVectorMdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseVectorMdfClass<V,M>::~BaseVectorMdfClass()
{
}
// Math methods--------------------------------------
template<class V, class M>
const VectorSetClass<V,M>&
BaseVectorMdfClass<V,M>::domainSet() const
{
  return m_domainSet;
}

//*****************************************************
// Generic marginal density function
//*****************************************************
/*!\class GenericVectorMdfClass
 * \brief A class for handling generic MDFs of vector functions. */

template<class V, class M>
class GenericVectorMdfClass : public BaseVectorMdfClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, and a routine 
   * (acting as a math function). */  
  GenericVectorMdfClass(const char*                    prefix,
                          const VectorSetClass<V,M>& domainSet,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
                          const void* routineDataPtr);
  //! Destructor
 ~GenericVectorMdfClass();
  //@}

    //! @name Mathematical method
  //@{
  //! Finds the values of the vector MDF at each element of \c paramValues, by calling \c m_routinePtr, and saves it at \c mdfValues.
  void values(const V& paramValues, V& mdfVec) const;
  //@}
  
  //! @name I/O method
  //@{ 
  //! TODO: Prints the vector MDF. 
  /*! \todo: implement me!*/
  void print (std::ostream& os)                const;
  //@}

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec);
  const void* m_routineDataPtr;

  using BaseVectorMdfClass<V,M>::m_env;
  using BaseVectorMdfClass<V,M>::m_prefix;
  using BaseVectorMdfClass<V,M>::m_domainSet;
};
// Default constructor -----------------------------
template<class V, class M>
GenericVectorMdfClass<V,M>::GenericVectorMdfClass(
  const char*                    prefix,
  const VectorSetClass<V,M>& domainSet,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
  const void* routineDataPtr)
  :
  BaseVectorMdfClass<V,M>(prefix,domainSet),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}
// Destructor ---------------------------------------
template<class V, class M>
GenericVectorMdfClass<V,M>::~GenericVectorMdfClass()
{
}
// Math methods--------------------------------------
template<class V, class M>
void
GenericVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  m_routinePtr(paramValues, m_routineDataPtr, mdfVec);
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
GenericVectorMdfClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// Gaussian marginal density function class
//*****************************************************
/*! 
 * \class GaussianVectorMdfClass
 * \brief TODO: A class for handling Gaussian MDFs.
 *
 * This class \b will implement a Gaussian vector marginal density function function (MDF).
 * \todo: Implement me! */

template<class V, class M>
class GaussianVectorMdfClass : public BaseVectorMdfClass<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   * Instantiates an object of the class given a prefix, the domain set, and the domain mean 
   * and expected values. */
  GaussianVectorMdfClass(const char*                    prefix,
                           const VectorSetClass<V,M>& domainSet,
                           const V&                       domainExpectedValues,
                           const V&                       domainVarianceValues);
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   * Instantiates an object of the class given a prefix, the domain set, and the domain mean 
   * and covariance matrix. */
  GaussianVectorMdfClass(const char*                    prefix,
                           const VectorSetClass<V,M>& domainSet,
                           const V&                       domainExpectedValues,
                           const M&                       covMatrix);
  //! Destructor
  ~GaussianVectorMdfClass();
  //@}

  //! @name Mathematical method
  //@{
  //! TODO: Returns the values of the vector CDF at each element of \c paramValues. 
  /*! \todo: implement me!*/
  void values(const V& paramValues, V& mdfVec) const;
  //@}
  
  //! @name I/O method
  //@{ 
  //! TODO: Prints the vector CDF. 
  /*! \todo: implement me!*/
  void print (std::ostream& os)                const;
  //@}

protected:
  const M*                         m_covMatrix;

  using BaseVectorMdfClass<V,M>::m_env;
  using BaseVectorMdfClass<V,M>::m_prefix;
  using BaseVectorMdfClass<V,M>::m_domainSet;

  void commonConstructor();
};
// Constructor -------------------------------------
template<class V,class M>
GaussianVectorMdfClass<V,M>::GaussianVectorMdfClass(
  const char*                    prefix,
  const VectorSetClass<V,M>& domainSet,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  BaseVectorMdfClass<V,M>(prefix,domainSet),
  m_covMatrix              (m_domainSet.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorMdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorMdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Constructor -------------------------------------
template<class V,class M>
GaussianVectorMdfClass<V,M>::GaussianVectorMdfClass(
  const char*                    prefix,
  const VectorSetClass<V,M>& domainSet,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  BaseVectorMdfClass<V,M>(prefix,domainSet),
  m_covMatrix              (new M(covMatrix))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorMdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorMdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor -------------------------------------
template<class V,class M>
GaussianVectorMdfClass<V,M>::~GaussianVectorMdfClass()
{
  delete m_covMatrix;
}
// Math method -------------------------------------
template<class V, class M>
void
GaussianVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "GaussianVectorMdfClass<V,M>::mdfVec()",
                      "incomplete code");
  return;
}
// I/O method --------------------------------------
template <class V, class M>
void
GaussianVectorMdfClass<V,M>::print(std::ostream& os) const
{
  return;
}
// Protected member function------------------------
template<class V,class M>
void
GaussianVectorMdfClass<V,M>::commonConstructor()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "GaussianVectorMdfClass<V,M>::commonConstructor()",
                      "incomplete code");
  return;
}

//*****************************************************
// Sampled marginal density function class
//*****************************************************
/*! 
 * \class SampledVectorMdfClass
 * \brief A class for handling sampled vector MDFs.
 *
 * This class implements a sampled vector marginal density function (MDF), given 
 * the grid points where it will be sampled and it returns its values.*/

template<class V, class M>
class SampledVectorMdfClass : public BaseVectorMdfClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix and the grid points 
   * where it will be sampled/evaluated.*/
  SampledVectorMdfClass(const char*                          prefix,
                          const ArrayOfOneDGridsClass <V,M>& oneDGrids,
                          const ArrayOfOneDTablesClass<V,M>& mdfValues);
  //! Destructor
  ~SampledVectorMdfClass();
  //@}
  
  //! @name Mathematical methods
  //@{
  //! TODO: Returns the values of the vector MDF at each element of \c paramValues. 
  /*! \todo: implement me!*/
  void values(const V& paramValues, V& mdfVec) const;
  //@}
  
  //! @name I/O methods
  //@{ 
  //! Prints the vector MDF (values of the grid points and of the MDF at such grid points).
  void print (std::ostream& os)                const;
  //@}

protected:
  using BaseVectorMdfClass<V,M>::m_env;
  using BaseVectorMdfClass<V,M>::m_prefix;
  using BaseVectorMdfClass<V,M>::m_domainSet;

  const ArrayOfOneDGridsClass <V,M>& m_oneDGrids;
  const ArrayOfOneDTablesClass<V,M>& m_mdfValues;
};
// Default constructor -----------------------------
template<class V,class M>
SampledVectorMdfClass<V,M>::SampledVectorMdfClass(
  const char*                          prefix,
  const ArrayOfOneDGridsClass <V,M>& oneDGrids,
  const ArrayOfOneDTablesClass<V,M>& mdfValues)
  :
  BaseVectorMdfClass<V,M>(prefix,oneDGrids.rowSpace()),
  m_oneDGrids(oneDGrids),
  m_mdfValues(mdfValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SampledVectorMdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SampledVectorMdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
SampledVectorMdfClass<V,M>::~SampledVectorMdfClass()
{
}
// Math methods--------------------------------------
template<class V, class M>
void
SampledVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "SampledVectorMdfClass<V,M>::mdfVec()",
                      "incomplete code");
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
SampledVectorMdfClass<V,M>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_oneDGrids;

  // Print *mdf* values *at* grid points
  os << m_mdfValues;

  return;
}

}  // End namespace QUESO

#endif // __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__
