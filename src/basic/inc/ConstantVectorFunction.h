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

#ifndef UQ_CONSTANT_VECTOR_FUNCTION_H
#define UQ_CONSTANT_VECTOR_FUNCTION_H

#include <queso/Environment.h>
#include <queso/DistArray.h>
#include <queso/VectorSet.h>
#include <queso/VectorFunction.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Constant class
//*****************************************************

/*!\class ConstantVectorFunction
 * \brief A class for handling vector functions which image is constant.
 *
 * This class allows the mathematical definition of a vector-valued function which image
 * set is constant vector. */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector ,class Q_M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_CONSTANT_VECTOR_FUNCTION_H
