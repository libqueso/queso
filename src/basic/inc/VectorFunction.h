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

#ifndef __UQ_VECTOR_FUNCTION_H__
#define __UQ_VECTOR_FUNCTION_H__

#include <queso/Environment.h>
#include <queso/DistArray.h>
#include <queso/VectorSet.h>

namespace QUESO {

//*****************************************************
// Base class
//*****************************************************


/*! \file uqVectorFunction.h
 * \brief Set of classes for handling vector functions.
 * 
 * \class BaseVectorFunction
 * \brief A templated (base) class for handling vector functions.
 *
 * This class allows the mathematical definition of a vector function such as:
 * \f$ \mathbf{q}: B \subset R^n \rightarrow R^m \f$. It requires the specification 
 * of the domain \f$ B \f$, which is a subset of the vector space (set) \f$ R^n \f$,
 * (and have already been introduced by the class VectorSet) and  of the 
 * image set \f$ R^m \f$.*/


template<class P_V,class P_M,class Q_V,class Q_M>
class BaseVectorFunction {
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default Constructor
  /*! Instantiates an object of the class, i.e. a vector function, given a prefix, its domain and image.*/
  BaseVectorFunction(const char*                      prefix,
			    const VectorSet<P_V,P_M>& domainSet,
			    const VectorSet<Q_V,Q_M>& imageSet);
  //! Destructor	   
  virtual ~BaseVectorFunction();
  //@}
  
  //! @name Mathematical methods.
  //@{  
  //! Access to the protected attribute \c m_domainSet: domain set of the vector function. It is an instance of the class VectorSet.  
  /*! This is one example of the advantage of how QUESO represents mathematical 
   * identities in a straightforward manner. */
  const VectorSet<P_V,P_M>& domainSet() const;
  
  //! Access to the protected attribute \c m_imageSet: image set of the vector function/ It is an instance of the class VectorSet.   
  const VectorSet<Q_V,Q_M>& imageSet () const;
  
  //! Computes the image vector. See template specialization.
  virtual void  compute  (const P_V&              domainVector,
			  const P_V*              domainDirection,
			  Q_V&                    imageVector,
			  DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
			  DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'	
			  DistArray<P_V*>* hessianEffects) const = 0;
  //@}
protected:
  const BaseEnvironment&    m_env;
        std::string                m_prefix;
	
  //! Domain set of the vector function.	
  const VectorSet<P_V,P_M>& m_domainSet;
  
  //! Image set of the vector function.
  const VectorSet<Q_V,Q_M>& m_imageSet;
};
// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::BaseVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"func_"),
  m_domainSet(domainSet),
  m_imageSet (imageSet)
{
}
// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::~BaseVectorFunction()
{
}
// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
const VectorSet<P_V,P_M>&
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_domainSet;
}
// --------------------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
const VectorSet<Q_V,Q_M>&
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::imageSet() const
{
  return m_imageSet;
}

//*****************************************************
// Generic class
//*****************************************************

/*!\class GenericVectorFunction
 * \brief A class for handling generic vector functions.
 *
 * This class allows the mathematical definition of a vector function such as:
 * \f$ \mathbf{q}: B \subset R^n \rightarrow R^m \f$. It is derived from 
 * BaseVectorFunction.
 */

template<class P_V,class P_M,class Q_V,class Q_M>
class GenericVectorFunction : public BaseVectorFunction<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Default constructor. 
  /*! Instantiates an object of \c this class given its prefix, domain and a pointer to a routine. 
   This routine plays the role of a vector-valued function, and it is useful, for instance, to 
   calculate the likelihood (and its image set).*/
  GenericVectorFunction(const char*                      prefix,
                               const VectorSet<P_V,P_M>& domainSet,
                               const VectorSet<Q_V,Q_M>& imageSet,
                               void (*routinePtr)(const P_V&                    domainVector,
                                                  const P_V*                    domainDirection,
                                                  const void*                   functionDataPtr,
                                                        Q_V&                    imageVector,
                                                        DistArray<P_V*>* gradVectors,
                                                        DistArray<P_M*>* hessianMatrices,
                                                        DistArray<P_V*>* hessianEffects),
                               const void* functionDataPtr);
  //! Virtual destructor.
  virtual ~GenericVectorFunction();

  
  //! @name Mathematical method
  //@{ 
  //! Calls the protected  member \c *m_routinePtr(), in order to calculate the image vector \c imageVector. 
  void compute  (const P_V&                    domainVector,
                 const P_V*                    domainDirection,
                       Q_V&                    imageVector,
                       DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                       DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                       DistArray<P_V*>* hessianEffects) const;
  //@}
		       
protected:
  //! Routine defining a vector-valued function.
  /*! The presence of the parameters \c gradVectors, \c hessianMatrices and \c hessianEffects
   * allows the user to calculate gradient vectors, Hessian matrices and Hessian effects; which 
   * can hold important information about her/his statistical application. Used, for instance to 
   * define the likelihood.  */
  void (*m_routinePtr)(const P_V&                    domainVector,
                       const P_V*                    domainDirection,
                       const void*                   functionDataPtr,
                             Q_V&                    imageVector,
                             DistArray<P_V*>* gradVectors,
                             DistArray<P_M*>* hessianMatrices,
                             DistArray<P_V*>* hessianEffects);
  const void* m_routineDataPtr;

  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_domainSet;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_imageSet;
};
// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::GenericVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  void (*routinePtr)(const P_V&                    domainVector,
                     const P_V*                    domainDirection,
                     const void*                   functionDataPtr,
                           Q_V&                    imageVector,
                           DistArray<P_V*>* gradVectors,
                           DistArray<P_M*>* hessianMatrices,
                           DistArray<P_V*>* hessianEffects),
  const void* functionDataPtr)
  :
  BaseVectorFunction<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}
// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::~GenericVectorFunction()
{
}
// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
void
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::compute(
  const P_V&                    domainVector,
  const P_V*                    domainDirection,
        Q_V&                    imageVector,
        DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        DistArray<P_V*>* hessianEffects) const
{
  //UQ_FATAL_TEST_MACRO(false,
  //                    domainVector.env().worldRank(),
  //                    "GenericVectorFunction<P_V,P_M,Q_V,Q_M>::compute()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(domainVector, domainDirection, m_routineDataPtr, imageVector, gradVectors, hessianMatrices, hessianEffects);

  return;
}

//*****************************************************
// Constant class
//*****************************************************

/*!\class ConstantVectorFunction
 * \brief A class for handling vector functions which image is constant.
 *
 * This class allows the mathematical definition of a vector-valued function which image 
 * set is constant vector. */

template<class P_V,class P_M,class Q_V,class Q_M>
class ConstantVectorFunction : public BaseVectorFunction<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default Constructor
  /*! Instantiates an object of the class, i.e. a vector function, given a prefix, its domain and constant image.*/
  ConstantVectorFunction(const char*                      prefix,
                                const VectorSet<P_V,P_M>& domainSet,
                                const VectorSet<Q_V,Q_M>& imageSet,
                                const Q_V&                       constantImageVector);
  //! Destructor
  virtual ~ConstantVectorFunction();
  //@}
  
  //! @name Mathematical method
  //@{ 
  //! Calculates the image vector: assigns to the protected attribute \c m_constantImageVector the value of the constant vector \c imageVector. 
  void compute  (const P_V&                    domainVector,
                 const P_V*                    domainDirection,
                       Q_V&                    imageVector,
                       DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                       DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                       DistArray<P_V*>* hessianEffects) const;
  //@}
protected:
  const Q_V* m_constantImageVector;

  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_domainSet;
  using BaseVectorFunction<P_V,P_M,Q_V,Q_M>::m_imageSet;
};
// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::ConstantVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  const Q_V&                       constantImageVector)
  :
  BaseVectorFunction<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
  m_constantImageVector(NULL)
{
  m_constantImageVector = new Q_V(constantImageVector);
}
// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::~ConstantVectorFunction()
{
  delete m_constantImageVector;
}
// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
void
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::compute(
  const P_V&                    domainVector,
  const P_V*                    domainDirection,
        Q_V&                    imageVector,
        DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        DistArray<P_V*>* hessianEffects) const
{
  UQ_FATAL_TEST_MACRO(m_constantImageVector == NULL,
                      domainVector.env().worldRank(),
                      "ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::compute()",
                      "m_constantImageVector is NULL");

  imageVector = *m_constantImageVector;

  return;
}

}  // End namespace QUESO

#endif // __UQ_VECTOR_FUNCTION_H__
