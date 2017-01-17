//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef UQ_BASE_VECTOR_FUNCTION_H
#define UQ_BASE_VECTOR_FUNCTION_H

#include <queso/Environment.h>
#include <queso/DistArray.h>
#include <queso/VectorSet.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorFunction.h
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

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_BASE_VECTOR_FUNCTION_H
