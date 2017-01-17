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

#ifndef QUESO_TK_FACTORY_LOGIT_RANDOM_WALK_H
#define QUESO_TK_FACTORY_LOGIT_RANDOM_WALK_H

#include <queso/TransitionKernelFactory.h>
#include <queso/TKGroup.h>

namespace QUESO
{

/**
 * TKFactoryLogitRandomWalk class defintion.  Implements a factory for
 * random-walk type transition kernels along with a variable transformation to
 * allow more efficient sampling for problems with parameters bounds in
 * high-dimensional state spaces.
 */
template <class DerivedTK>
class TKFactoryLogitRandomWalk : public TransitionKernelFactory
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  TKFactoryLogitRandomWalk(const std::string & name)
    : TransitionKernelFactory(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~TKFactoryLogitRandomWalk() {}

protected:
  virtual SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type build_tk()
  {
    SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk;

    // Cast the domain to a box.  Might this cast fail?
    const BoxSubset<GslVector, GslMatrix> & boxSubset =
      dynamic_cast<const BoxSubset<GslVector, GslMatrix> & >(
          this->m_target_pdf->domainSet());

    new_tk.reset(new DerivedTK(this->m_options->m_prefix.c_str(),
                               boxSubset,
                               *(this->m_dr_scales),
                               *(this->m_initial_cov_matrix)));

    return new_tk;
  }

};

} // namespace QUESO

#endif  // QUESO_TK_FACTORY_LOGIT_RANDOM_WALK_H
