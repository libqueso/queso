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

#ifndef UQ_GENERIC_SCALAR_FUNCTION_H
#define UQ_GENERIC_SCALAR_FUNCTION_H

#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/Defines.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\class GenericScalarFunction
 * \brief A class for handling generic scalar functions.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ f: B \subset R \rightarrow R \f$. It is derived from BaseScalarFunction. */

template <class V = GslVector, class M = GslMatrix>
class GenericScalarFunction : public BaseScalarFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of \c this class given its prefix, domain set and a pointer to a routine.
   This routine plays the role of a scalar math function, and it is useful, for instance, to calculate
   the likelihood (and its image set).*/
  GenericScalarFunction(const char*                  prefix,
                               const VectorSet<V,M>& domainSet,
                               double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
                               const void* routinesDataPtr,
                               bool routineIsForLn);
  //! Virtual destructor
  virtual ~GenericScalarFunction();

  //! @name Mathematical method
  //@{
  //! Calculates the actual value of this scalar function.
  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Calculates the logarithm of value of this scalar function.
  /*! It is used in routines that calculate the likelihood and expect the logarithm of value.*/
  double lnValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;

   //! Routine defining a scalar function.
   /*! The presence of the parameters \c gradVectors, \c hessianMatrices and \c hessianEffects
   * allows the user to calculate gradient vectors, Hessian matrices and Hessian effects; which
   * can hold important information about her/his statistical application. Used, for instance to
   * define the likelihood.  */
  double (*m_valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect);
  const void* m_routinesDataPtr;
  bool m_routineIsForLn;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_SCALAR_FUNCTION_H
