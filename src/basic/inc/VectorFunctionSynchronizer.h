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

#ifndef UQ_VECTOR_FUNCTION_SYNCHRONIZER_H
#define UQ_VECTOR_FUNCTION_SYNCHRONIZER_H

#include <queso/VectorFunction.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorFunctionSynchronizer.h
 * \brief Class for synchronizing the calls of vector-valued functions
 *
 * \class VectorFunctionSynchronizer
 * \brief A templated class for synchronizing the calls of vector-valued functions.
 *
 * This class creates a synchronization point among processes which call vector-valued
 * functions. This means that all processes must reach a point in their code before they
 * can all begin executing again. */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class VectorFunctionSynchronizer
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  VectorFunctionSynchronizer(const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& inputFunction,
                                    const P_V&                                        auxPVec,
                                    const Q_V&                                        auxQVec);
  //! Destructor
 ~VectorFunctionSynchronizer();
  //@}

  //! @name Mathematical methods
  //@{
  //! Access to the domain set of the vector-valued function which will be synchronized.
  const VectorSet<P_V,P_M>& domainSet() const;
  //@}

  //! @name Sync method
  //@{
  //! Calls the vector-valued function which will be synchronized.
  /*! This procedure  forms a barrier, and no processes in the communicator can pass the
   * barrier until all of them call the function. */
  void callFunction(const P_V*                    vecValues,
                    const P_V*                    vecDirection,
                          Q_V*                    imageVector,
                          DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                          DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                          DistArray<P_V*>* hessianEffects) const;
  //@}
private:
  const BaseEnvironment&                     m_env;
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& m_vectorFunction;
  const P_V&                                        m_auxPVec;
  const Q_V&                                        m_auxQVec;
};

}  // End namespace QUESO

#endif  // UQ_VECTOR_FUNCTION_SYNCHRONIZER_H
