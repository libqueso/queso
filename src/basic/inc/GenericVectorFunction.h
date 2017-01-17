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

#ifndef UQ_GENERIC_VECTOR_FUNCTION_H
#define UQ_GENERIC_VECTOR_FUNCTION_H

#include <queso/Environment.h>
#include <queso/DistArray.h>
#include <queso/VectorSet.h>
#include <queso/VectorFunction.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\class GenericVectorFunction
 * \brief A class for handling generic vector functions.
 *
 * This class allows the mathematical definition of a vector function such as:
 * \f$ \mathbf{q}: B \subset R^n \rightarrow R^m \f$. It is derived from
 * BaseVectorFunction.
 */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_GENERIC_VECTOR_FUNCTION_H
