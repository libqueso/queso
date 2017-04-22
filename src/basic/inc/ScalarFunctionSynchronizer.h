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

#ifndef UQ_SCALAR_FUNCTION_SYNCHRONIZER_H
#define UQ_SCALAR_FUNCTION_SYNCHRONIZER_H

#include <queso/Environment.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class V, class M>
class BayesianJointPdf;

/*! \file ScalarFunctionSynchronizer.h
 * \brief Class for synchronizing the calls of scalar functions
 *
 * \class ScalarFunctionSynchronizer
 * \brief A templated class for synchronizing the calls of scalar functions (BaseScalarFunction and derived classes).
 *
 * This class creates a synchronization point among processes which call scalar functions.
 * This means that all processes must reach a point in their code before they can all begin
 * executing again. */

template <class V = GslVector, class M = GslMatrix>
class ScalarFunctionSynchronizer
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  ScalarFunctionSynchronizer(const BaseScalarFunction<V,M>& inputFunction,
                                    const V&                              auxVec);

  //! Destructor.
 ~ScalarFunctionSynchronizer();
  //@}

  //! @name Mathematical methods
  //@{
  //! Access to the domain set of the scalar function which will be synchronized.
  const VectorSet<V,M>& domainSet() const;
  //@}

  //! @name Sync method
  //@{
  //! Calls the scalar function which will be synchronized.
  /*! This procedure  forms a barrier, and no processes in the communicator can pass the
   * barrier until all of them call the function. */
  double callFunction(const V* vecValues,
                      const V* vecDirection,
                            V* gradVector,
                            M* hessianMatrix,
                            V* hessianEffect,
                            double* extraOutput1,
                            double* extraOutput2) const;

  double callFunction(const V* vecValues,
                      double* extraOutput1,
                      double* extraOutput2) const;
  //@}
private:
  const BaseEnvironment&         m_env;
  const BaseScalarFunction<V,M>& m_scalarFunction;
  const BayesianJointPdf<V,M>*   m_bayesianJointPdfPtr;
  const V&                              m_auxVec;
};

}  // End namespace QUESO

#endif // UQ_SCALAR_FUNCTION_SYNCHRONIZER_H
