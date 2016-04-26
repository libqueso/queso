//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TransitionKernelFactory.h>
#include <queso/ScaledCovMatrixTKGroup.h>

namespace QUESO
{

template <>
std::map<std::string, Factory<BaseTKGroup<GslVector, GslMatrix> > *> &
Factory<BaseTKGroup<GslVector, GslMatrix> >::factory_map()
{
  static std::map<std::string, Factory<BaseTKGroup<GslVector, GslMatrix> > *> _factory_map;

  return _factory_map;
}

SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type
TransitionKernelFactory::build_tk(
    const MhOptionsValues & options,
    const VectorSpace<GslVector, GslMatrix> & v,
    const std::vector<double> & dr_scales,
    const ScalarFunctionSynchronizer<GslVector, GslMatrix> & pdf_synchronizer,
    const GslMatrix & initial_cov_matrix)
{
  SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk(
      new ScaledCovMatrixTKGroup<GslVector, GslMatrix>(
        "", v, dr_scales, initial_cov_matrix));

  return new_tk;
}

template<>
const VectorSpace<GslVector, GslMatrix> * FactoryWithVectorSpace<BaseTKGroup<GslVector, GslMatrix> >::m_vectorSpace = NULL;

const std::vector<double> * TransitionKernelFactory::m_dr_scales = NULL;

const ScalarFunctionSynchronizer<GslVector, GslMatrix> * TransitionKernelFactory::m_pdf_synchronizer = NULL;

const GslMatrix * TransitionKernelFactory::m_initial_cov_matrix = NULL;

const MhOptionsValues * TransitionKernelFactory::m_options = NULL;

} // namespace QUESO
