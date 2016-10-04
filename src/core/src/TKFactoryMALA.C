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

#include <queso/TKFactoryMALA.h>
#include <queso/MetropolisAdjustedLangevinTK.h>
#include <queso/BayesianJointPdf.h>

namespace QUESO
{

template <class DerivedTK>
SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type
TKFactoryMALA<DerivedTK>::build_tk()
{
  SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk;

  // Assume the problem is Bayesian
  const BayesianJointPdf<GslVector, GslMatrix> * target_bayesian_pdf =
    dynamic_cast<const BayesianJointPdf<GslVector, GslMatrix> *>(
        this->m_target_pdf);

  new_tk.reset(new DerivedTK(this->m_options->m_prefix.c_str(),
                             *target_bayesian_pdf,
                             *(this->m_dr_scales),
                             *(this->m_initial_cov_matrix)));

  return new_tk;
}

// Instantiate all the transition kernel factories
TKFactoryMALA<MetropolisAdjustedLangevinTK<GslVector, GslMatrix> > tk_factory_mala("mala");

} // namespace QUESO
