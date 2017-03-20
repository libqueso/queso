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

#ifndef QUESO_CUSTOM_TK_HESS_H
#define QUESO_CUSTOM_TK_HESS_H

#include <queso/ScaledCovMatrixTKGroup.h>
#include <queso/VectorRV.h>
#include <queso/ScalarFunctionSynchronizer.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class V = GslVector, class M = GslMatrix>
class MyTransitionKernel : public ScaledCovMatrixTKGroup<V, M> {
public:
  MyTransitionKernel(const char * prefix,
                     const VectorSpace<V, M> & vectorSpace,
                     const std::vector<double> & scales,
                     const M & covMatrix);

  virtual ~MyTransitionKernel();
  virtual void updateTK();
  virtual bool covMatrixIsDirty();
  virtual void cleanCovMatrix();

  const VectorSpace<V, M> & m_vectorSpace;

private:
  unsigned int m_counter;
  bool m_dirtyCovMatrix;
};

}  // End namespace QUESO

#endif // QUESO_CUSTOM_TK_HESS_H
