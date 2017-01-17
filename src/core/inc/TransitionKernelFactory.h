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

#ifndef QUESO_TK_FACTORY_H
#define QUESO_TK_FACTORY_H

#include <queso/Factory.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/TKGroup.h>
#include <queso/MetropolisHastingsSGOptions.h>

namespace QUESO
{

/**
 * TransitionKernelFactory class defintion.
 */
class TransitionKernelFactory : public Factory<BaseTKGroup<GslVector, GslMatrix> >
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  TransitionKernelFactory(const std::string & name)
    : Factory<BaseTKGroup<GslVector, GslMatrix> >(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~TransitionKernelFactory() {}

  /**
   * Static method to set the vector space the transition kernel is defined on
   */
  static void set_vectorspace(const VectorSpace<GslVector, GslMatrix> & v)
  {
    m_vectorSpace = &v;
  }

  /**
   * Static method to set the vector of scale factor to scale a proposal
   * covariance matrix by for the purpose of delayed rejection
   */
  static void set_dr_scales(const std::vector<double> & scales)
  {
    m_dr_scales = &scales;
  }

  /**
   * Static method to set the pdf synchronizer.  Used by Stochastic Newton?
   */
  static void set_pdf_synchronizer(const ScalarFunctionSynchronizer<GslVector, GslMatrix> & synchronizer)
  {
    m_pdf_synchronizer = &synchronizer;
  }

  /**
   * Static method to set the initial proposal covariance matrix
   */
  static void set_initial_cov_matrix(GslMatrix & cov_matrix)
  {
    m_initial_cov_matrix = &cov_matrix;
  }

  /**
   * Static method to set the options object.  Useful for passing input file
   * options to transition kernels.
   */
  static void set_options(const MhOptionsValues & options)
  {
    m_options = &options;
  }

  /**
   * Static method to set the pdf we wish to draw samples from
   */
  static void set_target_pdf(const BaseJointPdf<GslVector, GslMatrix> & target_pdf)
  {
    m_target_pdf = &target_pdf;
  }

protected:
  virtual SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type build_tk() = 0;

  static const VectorSpace<GslVector, GslMatrix> * m_vectorSpace;
  static const std::vector<double> * m_dr_scales;
  static const ScalarFunctionSynchronizer<GslVector, GslMatrix> * m_pdf_synchronizer;
  static GslMatrix * m_initial_cov_matrix;
  static const MhOptionsValues * m_options;
  static const BaseJointPdf<GslVector, GslMatrix> * m_target_pdf;

private:
  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type create();
};

inline
SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type TransitionKernelFactory::create()
{
  queso_require_msg(m_vectorSpace, "ERROR: must call set_vectorspace() before building tk!");
  queso_require_msg(m_dr_scales, "ERROR: must call set_dr_scales() before building tk!");
  queso_require_msg(m_pdf_synchronizer, "ERROR: must call set_pdf_synchronizer() before building tk!");
  queso_require_msg(m_initial_cov_matrix, "ERROR: must call set_initial_cov_matrix() before building tk!");
  queso_require_msg(m_options, "ERROR: must call set_options() before building tk!");
  queso_require_msg(m_target_pdf, "ERROR: must call set_target_pdf() before building tk!");

  SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk = this->build_tk();

  queso_assert(new_tk);

  return new_tk;
}

} // namespace QUESO

#endif // QUESO_TK_FACTORY_H
