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

#include <queso/TKFactoryStochasticNewton.h>
#include <queso/HessianCovMatricesTKGroup.h>

namespace QUESO
{

template <class DerivedTK>
SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type
TKFactoryStochasticNewton<DerivedTK>::build_tk()
{
  SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk;

  new_tk.reset(new DerivedTK(this->m_options->m_prefix.c_str(),
                             *(this->m_vectorSpace),
                             *(this->m_dr_scales),
                             *(this->m_pdf_synchronizer)));

  return new_tk;
}

// Instantiate all the transition kernel factories
TKFactoryStochasticNewton<HessianCovMatricesTKGroup<GslVector, GslMatrix> > tk_factory_stochastic_newton("stochastic_newton");

} // namespace QUESO
