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

#include "adaptive_metropolis.h"

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianJointPdf.h>
#include <queso/TKFactoryRandomWalk.h>

namespace QUESO {

template<class V, class M>
MyTransitionKernel<V, M>::MyTransitionKernel(
  const char*                    prefix,
  const VectorSpace<V, M>& vectorSpace, // FIX ME: vectorSubset ???
  const std::vector<double>&     scales,
  const M&                       covMatrix)
  :
  ScaledCovMatrixTKGroup<V, M>(prefix, vectorSpace, scales, covMatrix),
  m_vectorSpace(vectorSpace),
  m_counter(0),
  m_dirtyCovMatrix(false)
{
}

template<class V, class M>
MyTransitionKernel<V, M>::~MyTransitionKernel()
{
}

template <class V, class M>
void
MyTransitionKernel<V, M>::updateTK()
{
  if (m_counter == 15) {
    GslMatrix new_matrix(m_vectorSpace.zeroVector());

    new_matrix(0, 0) = 1.0;

    this->m_dirtyCovMatrix = true;

    // Update all the RVs with new cov matrix (including all the ones used in
    // Delayed Rejection)
    this->updateLawCovMatrix(new_matrix);
  }

  // Just keeping track of a counter so I can trigger a 'dirty' matrix
  m_counter++;
}

template <class V, class M>
bool
MyTransitionKernel<V, M>::covMatrixIsDirty()
{
  return m_dirtyCovMatrix;
}

template <class V, class M>
void
MyTransitionKernel<V, M>::cleanCovMatrix()
{
  this->m_dirtyCovMatrix = false;
}

// Explicit instantiation of the template
template class MyTransitionKernel<GslVector, GslMatrix>;

// Register this TK with the appropriate factory
TKFactoryRandomWalk<MyTransitionKernel<GslVector, GslMatrix> > tk_factory_mytk("my_tk");

}  // End namespace QUESO
